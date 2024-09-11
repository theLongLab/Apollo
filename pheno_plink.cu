#include "pheno_plink.cuh"

pheno_plink::pheno_plink(string vcf_File, string output_File,string output_VCF)
{
    cout << "Intializing Phenotype generation\n";
    functions_library functions = functions_library();

    fstream vcf_File_original;
    vcf_File_original.open(vcf_File, ios::in);

    fstream vcf_File_rename;
    vcf_File_rename.open(output_VCF, ios::out);

    string line;

    if (vcf_File_original.is_open())
    {
        cout << "Reading VCF: " << vcf_File << endl;
        while (getline(vcf_File_original, line))
        {
            if (line.substr(0, 2) != "##")
            {
                break;
                // vcf_File_rename << line << endl;
            }
            else
            {
                vcf_File_rename << line << endl;
            }
        }
    }
    else
    {
        cout << "ERROR: UNABLE TO OPEN VCF FILE: " << vcf_File << endl;
        exit(-1);
    }

    cout << "Completed reading VCF\n";
    functions.split(headers, line, '\t');
    cout << to_string(headers.size() - 9) << " samples found\n";

    for (int col = 0; col < 9; col++)
    {
        vcf_File_rename << headers[col] << "\t";
    }

    int count = 0;

    for (int sample = 9; sample < headers.size(); sample = sample + 2)
    {
        headers[sample] = to_string(count) + "vNORMAL";
        headers[sample + 1] = to_string(count) + "vCANCER";
        count++;
    }

    for (int col = 9; col < headers.size(); col++)
    {
        vcf_File_rename << headers[col];
        if (col + 1 != headers.size())
        {
            vcf_File_rename << "\t";
        }
        else
        {
            vcf_File_rename << "\n";
        }
    }

    while (getline(vcf_File_original, line))
    {
        vcf_File_rename << line << endl;
    }

    vcf_File_rename.close();
    vcf_File_original.close();

    this->output_File = output_File;
}

void pheno_plink::ingress()
{
    cout << "Generating phenotpye file: " << output_File << "\n";
    functions_library functions = functions_library();

    functions.create_File(this->output_File, "FID\tIID\tPhenotype");

    fstream pheno_File;
    pheno_File.open(output_File, ios::app);
    if (pheno_File.is_open())
    {
        vector<string> line_Data_1;
        vector<string> line_Data_2;
        int family_ID = 0;
        for (int sample = 9; sample < headers.size(); sample = sample + 2)
        {
            functions.split(line_Data_1, headers[sample], 'v');
            functions.split(line_Data_2, headers[sample + 1], 'v');

            int line_1 = 0;
            if (line_Data_1[1] != "NORMAL")
            {
                line_1 = 1;
            }

            int line_2 = 0;
            if (line_Data_2[1] != "NORMAL")
            {
                line_2 = 1;
            }

            pheno_File << headers[sample] << "\t" << headers[sample] << "\t" << line_1 << endl;
            pheno_File << headers[sample + 1] << "\t" << headers[sample + 1] << "\t" << line_2 << endl;
            family_ID++;
        }
    }
    else
    {
        cout << "ERROR: UNABLE TO OPEN PHENOTPYE FILE: " << output_File << endl;
        exit(-1);
    }
    pheno_File.close();

    cout << "Done generating phenotype file\n";
}