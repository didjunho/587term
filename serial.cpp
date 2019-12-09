#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>

using namespace std;

//efficient estimation formula for z-score --> probability
inline double z_norm(const double & z) {
    if(z > 3.6 || z < -3.6) { return 1.0; }
    double z_loc = z;
    if(z < 0) {
        z_loc *= -1;
        return 1 - (0.46375418 + 0.065687194*tanh(1.280022196 - 0.720528073*z_loc) - 0.602383931*tanh(0.033142223 - 0.682842425*z_loc));
    }
    return (0.46375418 + 0.065687194*tanh(1.280022196 - 0.720528073*z) - 0.602383931*tanh(0.033142223 - 0.682842425*z));
}

//pricing algorithm for black scholes on call options
inline double blackScholesCall(const double & S, const double & X, const double & T,
                               const double & r, const double & d, const double & v) {
    double d1, d2;
    d1 = ((log(S/X) + (r - d + pow(v,2)/2) * T) / v) / sqrt(T);
    d2 = d1 -  v * sqrt(T);
    return (exp(-d * T) * S * z_norm(d1)) - (X * exp(-r * T) * z_norm(d2));
}

//calculating implied volatility for an option
inline double impliedVolatility(const double & S, const double & X, const double & T,
                                const double & r, const double & d, const double & p) {
    double epsilon_abs = 0.0000001;
    double epsilon_step = 0.0000001;
    double low_vol = 0.001;
    double high_vol = 1; 
    int iterations = 0;
    double mid_vol = 0;
    double prev_mid_vol = -1;

    while((high_vol - low_vol >= epsilon_step) || (abs(blackScholesCall(S,X,T,r,d,low_vol)-p) >= epsilon_abs && abs(blackScholesCall(S, X, T, r, d, high_vol)-p) >= epsilon_abs)) {
        mid_vol = (low_vol + high_vol)/2;
        //identified implied volatility
        if(abs(blackScholesCall(S,X,T,r,d,mid_vol)-p) <= epsilon_abs) {
            break;
        }
        //reduce upper bound
        else if ((blackScholesCall(S,X,T,r,d,low_vol)-p) * (blackScholesCall(S,X,T,r,d,mid_vol)-p) < 0) {
            high_vol = mid_vol;
        }
        //increase lower bound
        else {
            low_vol = mid_vol;
        }
        ++iterations;
        if(mid_vol == prev_mid_vol) { break; }
        prev_mid_vol = mid_vol;
    }
    cout << "Iterations to converge: " << iterations << endl;
    return low_vol;
}

int main(int argc, char ** argv) {
    //ios_base::sync_with_stdio(false);
    //S = stock price
    //X = strike price
    //T = time to maturity in years
    //r = risk free rate
    //d = dividend yield
    //p = option price

    //READ IN DATA
    //start timing
    //one by one calculate each one, append to a list
    //end timing
    //generate surfaces
    vector<double> X;
    for(int i = 34; i < 38; ++i) {
        X.push_back(i);
    }
    vector<double> T;
    for(int i = 2; i < 3; ++i) {
        T.push_back(i);
    }
    double S = 30;
    double r = 0.05;
    double d = 0.02;
    double p = 3;
    cout << impliedVolatility(S, 35, 2, r, d, p) << endl;
    for(int i = 0; i < X.size(); ++i) {
        for(int j = 0; j < T.size(); ++j) {
            cout << "X: " << X[i] << " T: " << T[j] << " iv: " << impliedVolatility(S, X[i], T[i], r, d, p) << endl;
        }
    }
    return 0;
} 