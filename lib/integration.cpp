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
