#include <iostream>
#include <armadillo>
#include <cmath>

using namespace std;
using namespace arma;

double driftfun(double y)
{
    return(y/2.0);
}

double difffun(double y)
{
    return(sqrt(1 + y*y));
}

int main(void)
{
    double h = 0.002;
    double s = 0.75;
    double k = pow(h,s);
    double yM = datum::pi/k;
    int bigm = ceil(yM/k);
    unsigned int veclen = 2*bigm+1;
    double init = 0;

    vec phatn = zeros<vec>(veclen);
    vec phatnp1 = zeros<vec>(veclen);
    
    // pdf after one time step
    double mydiff;
    double temp;
    for (int i=-bigm;i<=bigm;i++)
    {
        mydiff = pow(difffun(init),2);
        phatn(i+bigm) = exp(-pow(i*k-init-driftfun(init)*h,2)/(2.0*mydiff*h))/sqrt(2.0*datum::pi*mydiff*h);
    }

    // iterate
    double T = 1.0;
    int bign = ceil(T/h);
    double thresh = 1e-16;
    for (int n=1;n<bign;n++)
    {
        cout << "n = " << n << "\n";
        cout.flush();

        #pragma omp parallel for
        for (int i=-bigm;i<=bigm;i++)
        {
            bool keepgoing = true;
            int j = i;
            double tally = 0.0;
            while (keepgoing)
            {
                double thisdiff = pow(difffun(j*k),2);
                double ker = k*exp(-pow(i*k - j*k - driftfun(j*k)*h,2)/(2.0*thisdiff*h))/sqrt(2.0*datum::pi*thisdiff*h);
                if (ker < thresh) keepgoing = false;
                tally += ker*phatn(j+bigm);
                j++;
                if (j > bigm) keepgoing = false;
            }
            if (i > -bigm)
            {
                keepgoing = true;
                j = i-1;
            }
            while (keepgoing)
            {
                double thisdiff = pow(difffun(j*k),2);
                double ker = k*exp(-pow(i*k - j*k - driftfun(j*k)*h,2)/(2.0*thisdiff*h))/sqrt(2.0*datum::pi*thisdiff*h);
                if (ker < thresh) keepgoing = false;
                tally += ker*phatn(j+bigm);
                j--;
                if (j < -bigm) keepgoing = false;
            }
            phatnp1(i+bigm) = tally;
        }
        #pragma omp barrier
        phatn = phatnp1;
    }
    cout << "\n\n" << "********" << "\n\n" << phatn << "\n";
}


