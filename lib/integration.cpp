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

std::vector<std::vector<double>> NumLib::RombergIntegration(double x1, double x2, size_t N, const std::function<double (double)> &f) {
    if(N < 1 ) N = 1;
    std::vector<std::vector<double>> romberg_integral(N, std::vector<double>(N));
    
    if(f) {

        //R(0,0) Start with trapezoidal integration with N = 1
        romberg_integral.front().front() = NumLib::TrapezoidalIntegration(x1, x2, 1, f);

        double h = x2-x1;
        for(size_t step = 1; step < N; step++) {
            h *= 0.5;

            //R(step, 0) Improve trapezoidal integration with decreasing h
            double trapezoidal_integration = 0;
            size_t stepEnd = pow(2, step - 1);
            for(size_t tzStep = 1; tzStep <= stepEnd; tzStep++) {
                const double deltaX = (2*tzStep - 1)*h;
                trapezoidal_integration += f(x1 + deltaX);
            }
            romberg_integral[step].front() = 0.5*romberg_integral[step - 1].front() + trapezoidal_integration*h;

            //R(m,n) Romberg integration with R(m,1) -> Simpson rule, R(m,2) -> Boole's rule
            for(size_t rbStep = 1; rbStep <= step; rbStep++) {
                const double k = pow(4, rbStep);
                romberg_integral[step][rbStep] = (k*romberg_integral[step][rbStep-1] - romberg_integral[step-1][rbStep-1])/(k-1);
            }
        }
    }

    return romberg_integral;
}

double NumLib::GaussLegendreIntegration::operator () (double x1, double x2, size_t N, const std::function<double (double)> &f) const {

}