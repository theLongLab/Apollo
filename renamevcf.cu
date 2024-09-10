#include "renamevcf.cuh"

renamevcf::renamevcf(string vcf_Folder, string sample_sheet_vcf)
{
    cout << "Starting TCGA VCF file and creating PLINK phenotype file\n\n";

    vcf_Folder_name = filesystem::path(vcf_Folder).filename().string();

    this->vcf_Folder = vcf_Folder;

    cout << "Processing vcf Folder: " << vcf_Folder_name << endl;

    cout << "Finding folder in sample sheet file\n";

    fstream sample_Sheet;
    sample_Sheet.open(sample_sheet_vcf, ios::in);

    functions_library functions = functions_library();

    string sample_Type = "";

    if (sample_Sheet.is_open())
    {
        string line;
        vector<string> line_Data;

        getline(sample_Sheet, line);

        char delim = '\t';

        functions.split(line_Data, line, delim);

        if (line_Data.size() == 1)
        {
            delim = ',';
        }

        while (getline(sample_Sheet, line))
        {
            functions.split(line_Data, line, delim);
            if (line_Data[0] == vcf_Folder_name)
            {
                cout << "VCF file found\n";
                sample_Type = line_Data[7];
                vcf_File_name_only = line_Data[1].substr(0, line_Data[1].find_last_of(".gz") - 2);
                vcf_File = vcf_Folder + "/" + line_Data[1].substr(0, line_Data[1].find_last_of(".gz") - 2);
                cout << "vcf File location: " << vcf_File << endl;
                break;
            }
        }
    }
    else
    {
        cout << "ERROR: UNABLE TO OPEN SAMPLE SHEET FILE: " << sample_sheet_vcf << endl;
        exit(-1);
    }

    if (sample_Type != "")
    {

        vector<string> sample_Types;
        functions.split(sample_Types, sample_Type, ',');

        if (sample_Types[0] == " Metastatic" || sample_Types[1] == " Metastatic" || sample_Types[0] == "Metastatic" || sample_Types[1] == "Metastatic")
        {
            tumor_Column_Name = vcf_Folder_name + "_METASTATIC_TUMOR";
        }
        else if (sample_Types[0] == "Primary Tumor" || sample_Types[1] == "Primary Tumor"||sample_Types[0] == " Primary Tumor" || sample_Types[1] == " Primary Tumor")
        {
            tumor_Column_Name = vcf_Folder_name + "_PRIMARY_TUMOR";
        }
        else
        {
            cout << "No tumor data found for the VCF\n";
            exit(-1);
        }
    }
    else
    {
        cout << "\nERROR: NO VCF folder match found\n";
        exit(-1);
    }
}

void renamevcf::ingress()
{
    if (filesystem::exists(vcf_Folder + "/renamed_" + vcf_File_name_only))
    {
        cout << "Already exists\n";
    }
    else
    {
        cout << "\nReading and renaming vcf files\n";
        fstream vcf_File_rename;
        fstream vcf_File_original;

        vcf_File_original.open(vcf_File, ios::in);

        functions_library functions = functions_library();

        if (vcf_File_original.is_open())
        {
            vcf_File_rename.open(vcf_Folder + "/renamed_" + vcf_File_name_only, ios::out);

            string line;

            while (getline(vcf_File_original, line))
            {
                if (line.substr(0, 2) == "##")
                {
                    vcf_File_rename << line << endl;
                }
                else
                {
                    break;
                }
            }

            vector<string> line_Data;
            functions.split(line_Data, line, '\t');

            for (int col = 0; col < 9; col++)
            {
                vcf_File_rename << line_Data[col] << "\t";
            }

            if (line_Data[9] == "NORMAL")
            {
                vcf_File_rename << vcf_Folder_name << "_NORMAL\t" << tumor_Column_Name;
            }
            else
            {
                vcf_File_rename << tumor_Column_Name << "\t" << vcf_Folder_name << "_NORMAL";
            }

            vcf_File_rename << endl;

            while (getline(vcf_File_original, line))
            {
                vcf_File_rename << line << endl;
            }

            cout << "Renaming VCF columns completed: ";

            vcf_File_rename.close();
            vcf_File_original.close();

            cout << vcf_Folder << "/renamed_" << vcf_File_name_only << endl;
        }
        else
        {
            cout << "ERROR: UNABLE TO OPEN FILE: " << vcf_File << endl;
        }
    }
}