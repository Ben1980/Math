#include "integration.h"
#include <iostream>

double NumLib::TrapezoidalIntegration(double x1, double x2, size_t N, const std::function<double (double)> &f) {
    if(f) {
        if(N < 1 ) N = 1;
        const double width = (x2-x1)/N;

        double trapezoidal_integral = 0;
        for(size_t step = 0; step < N; step++) {
            const double x_i = x1 + step*width;
            const double x_i1 = x1 + (step+1)*width;

            trapezoidal_integral += 0.5*(x_i1-x_i)*(f(x_i) + f(x_i1));
        }

        return trapezoidal_integral;
    }
    
    return 0;
}

double NumLib::SimpsonIntegration(double x1, double x2, size_t N, const std::function<double (double)> &f) {
    if(f) {
        if(N < 1 ) N = 1;
        const double width = (x2-x1)/N;

        double simpson_integral = 0;
        for(size_t step = 0; step < N; step++) {
            const double x_i = x1 + step*width;
            const double x_i1 = x1 + (step+1)*width;

            simpson_integral += (x_i1-x_i)/6.0*(f(x_i) + 4.0*f(0.5*(x_i+x_i1)) + f(x_i1));
        }

        return simpson_integral;
    }
    
    return 0;
}
