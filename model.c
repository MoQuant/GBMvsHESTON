#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double dWT(){
    int num = 20;
    double dw = (rand() % (2*num + 1)) - num;
    return dw/100.0;
}

double Heston(double kappa, double theta, double sigma, double S, double mu, double v, double dt, int N, int P){
    srand(time(NULL));
    double price = 0;
    for(int p = 0; p < P; ++p){
        double S0 = S;
        double v0 = v;
        for(int i = 0; i < N; ++i){
            double dw = dWT();
            v0 += kappa*(theta - v0)*dt + sigma*pow(v0, 0.5)*dw;
            S0 += mu*S0*dt + pow(v0, 0.5)*S0*dw;
        }
        price += S0;
    }
    return price / (double) P;
}

double GBM(double S, double mu, double v, double dt, int N, int P){
    srand(time(NULL));
    double price = 0;
    for(int p = 0; p < P; ++p){
        double S0 = S;
        for(int i = 0; i < N; ++i){
            S0 += mu*S0*dt + pow(v, 0.5)*S0*dWT();
        }
        price += S0;
    }
    return price / (double) P;
}

/*
int main()
{
    srand(time(NULL));

    double S = 100;
    double mu = 0.03;
    double vt = 0.05;
    double v = 0.05;
    double t = 30.0/365.0;
    int N = 1500;
    int P = 100;
    double dt = t / (double) N; 

    double kappa = 7.3;
    double theta = 0.00025;
    double sigma = 0.00007;

    for(int k = 0; k < 100; ++k){
        printf("%f | %f\n", Heston(kappa, theta, sigma, S, mu, vt, dt, N, P), GBM(S, mu, v, dt, N, P));
    }

    return 0;
}
*/