#include "test.cuh"

test::test()
{
    // seed
    unsigned seed = chrono::steady_clock::now().time_since_epoch().count();
    default_random_engine e(seed);
    cout << "Seed: " << seed << endl;
    // dist_Nomral(mean,std)
    normal_distribution<double> dist_Normal(5, 2);
    // dist_Gamma(alpha,beta)
    gamma_distribution<double> dist_Gamma(0.5, 2);

    // display values
    for (size_t i = 0; i < 1000; i++)
    {
        // cout << dist_Gamma(e) << endl;
    }
}