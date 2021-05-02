#include "integration.h"
#include <iostream>

using namespace NumLib;

double TrapezoidalIntegration(double x1, double x2, size_t N, const std::function<double (double)> &f) {
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

double SimpsonIntegration(double x1, double x2, size_t N, const std::function<double (double)> &f) {
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

std::vector<std::vector<double>> RombergIntegration(double x1, double x2, size_t N, const std::function<double (double)> &f) {
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

const double GaussLegendreIntegration::LegendrePolynomial::EPSILON = 1e-15;

const std::vector<double> & GaussLegendreIntegration::LegendrePolynomial::getWeight() const {
    return mWeight;
}

const std::vector<double> & GaussLegendreIntegration::LegendrePolynomial::getRoot() const {
    return mRoot;
}

void GaussLegendreIntegration::LegendrePolynomial::calculateWeightAndRoot() {
    for(int step = 0; step <= mNumberOfIterations; step++) {
        double root = cos(M_PI * (step-0.25)/(mNumberOfIterations+0.5));
        Result result = calculatePolynomialValueAndDerivative(root);

        double newtonRaphsonRatio;
        do {
            newtonRaphsonRatio = result.value/result.derivative;
            root -= newtonRaphsonRatio;
            result = calculatePolynomialValueAndDerivative(root);
        } while (fabs(newtonRaphsonRatio) > EPSILON);

        mRoot[step] = root;
        mWeight[step] = 2.0/((1-root*root)*result.derivative*result.derivative);
    }
}

GaussLegendreIntegration::LegendrePolynomial::Result GaussLegendreIntegration::LegendrePolynomial::calculatePolynomialValueAndDerivative(double x) {
    Result result(x, 0);

    double value_minus_1 = 1;
    const double f = 1/(x*x-1);
    for(int step = 2; step <= mNumberOfIterations; step++) {
        const double value = ((2*step-1)*x*result.value-(step-1)*value_minus_1)/step;
        result.derivative = step*f*(x*value - result.value);

        value_minus_1 = result.value;
        result.value = value;
    }

    return result;
}

double GaussLegendreIntegration::operator () (double x1, double x2, size_t N, const std::function<double (double)> &f) const {
    const LegendrePolynomial legendrePolynomial(x1, x2, N);
    const std::vector<double> & weight = legendrePolynomial.getWeight();
    const std::vector<double> & root = legendrePolynomial.getRoot();

    const double width = 0.5*(x2-x1);
    const double mean = 0.5*(x1+x2);

    double gaussLegendre = 0;
    for(int step = 1; step <= N; step++) {
        gaussLegendre += weight[step]*f(width * root[step] + mean);
    }

    return gaussLegendre * width;
}



// NumLib::GaussLegendreIntegration::LegendrePolynomial(double lowerBound, double upperBound, size_t numberOfIterations)
//     : mLowerBound(lowerBound), mUpperBound(upperBound), mNumberOfIterations(numberOfIterations), mWeight(numberOfIterations+1), mRoot(numberOfIterations+1) {
//     calculateWeightAndRoot();
// }