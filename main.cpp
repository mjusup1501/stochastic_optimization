#include <iostream>
#include <vector>
#include <random>

using namespace std;

vector<double> wiener_process(int T){
    vector<double> wiener(T, 0);
    normal_distribution<double> rnorm(0, 1);
    // random device class instance, source of 'true' randomness for initializing random seed
    std::random_device rd;
// Mersenne twister PRNG, initialized with seed from previous random device instance
    std::mt19937 generator(rd());

    for (int t = 0; t < T; ++t){
        double Wt = sqrt(t + 1) * rnorm(generator);
        wiener[t] = Wt;
    }

    return wiener;
}



vector<double> GBM(int T, double S0, double mu, double sigma){
    vector<double> gbm(T + 1, S0);
    vector<double> wiener = wiener_process(T);
    for (int t = 1; t <= T; ++t){
        double drift = (mu - pow(sigma, 2) / 2) * t;
        double vol = sigma * wiener[t - 1];
        double St = S0 * exp(drift + vol);
        gbm[t] = St;
    }

    return gbm;
}

double mean(vector<double> v){
    return accumulate(v.begin(), v.end(), 0) / double(v.size());
}

int print_vector(vector<double> v){
    for (auto elem : v) cout << elem << " ";
    cout << endl;
    return 0;
}

int main() {
    int sample_size = 10000;
    int T = 10;
    double S0 = 40;
    double mu = 0;
    double sigma = 0.25;

    /*
    vector<double> wiener = wiener_process(T);
    for (auto w : wiener) cout << w << endl;
    cout << mean(wiener);
    */

    vector<vector<double>> sample;
    for (int i = 0; i < sample_size; ++i){
        vector<double> gbm = GBM(T, S0, mu, sigma);
        sample.push_back(gbm);
    }

    // Expecting mean to be S0 because ES_t = S0 * exp(mu * t)
    vector<double> means(T + 1, S0);
    for (int t = 0; t <= T; ++t){
        double mean_t = 0;
        for (int i = 0; i < sample_size; ++i){
            mean_t += sample[i][t];
        }
        mean_t /= sample_size;
        means[t] = mean_t;
    }

    for (auto m : means){
        cout << m << endl;
    }

    return 0;
}
