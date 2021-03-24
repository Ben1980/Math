#ifndef INTEGRATION_H
#define INTEGRATION_H

#include <doctest/doctest.h>
#include <functional>

namespace NumLib {
  /**
   * @brief A numerical integration method using the trapezoidal approach
   * 
   * TrapezoidalIntegration is using the trapezoidal approach to acomplish numerical integrations as described at [thoughts-on-coding.com](https://thoughts-on-coding.com/2019/04/17/numerical-methods-in-cpp-part-1-newton-cotes-integration/)
   * and calculates the integral of function \f$f\f$ in range \f$[a-b]\f$.
   * 
   * \f[
   * w=\frac{b-a}{N}
   * \f]
   * \f[
   * x_i=a+i \cdot w
   * \f]
   * \f[
   * x_{i+1}=a+\left( i+1 \right) \cdot w
   * \f]
   * \f[
   * I=\textstyle \sum_{0}^{N}\frac{x_{i+1}-x_i}{2}\left( f(x_i)+f(x_{i+1}) \right)
   * \f]
   * 
   * @image html trapezoid.png
   * 
   * @test Test TrapezoidalIntegration, valid case: Testing method with function \f$I=\int_0^{\pi/2} \frac{5}{e^\pi-2}\exp(2x)\cos(x)dx=1.0\f$
   * @test Test TrapezoidalIntegration, x1=x2=0: Testing method with zero range
   * @test Test TrapezoidalIntegration, N=0: Testing method with zero iteration steps
   * @test Test TrapezoidalIntegration, f=nullptr: Testing method with error in function f
   * 
   * @param x1 begining of Integration range \f$x_1=a\f$
   * @param x2 end of Integration range \f$x_2=b\f$
   * @param N number of iteration steps
   * @param f function \f$f\f$ to integrate
   * 
   * @return Integral of function \f$f\f$ in range \f$[a-b]\f$
   */
  double TrapezoidalIntegration(double x1, double x2, size_t N, const std::function<double (double)> &f);
};

#include "testUtils.h"

TEST_CASE ("Test TrapezoidalIntegration, valid case") {
  const double PI = 3.14159265358979323846;
  const double E = 2.71828182845904523536;

  auto f = [E, PI](double x) { 
    return 5.0/(std::pow(E, PI) - 2.0) * exp(2.0*x) * cos(x); 
  };

  const double result = NumLib::TrapezoidalIntegration(0, 0.5*PI, 100, f);
  CHECK(TestUtils::IsValid(result, 1.0, 1e-3));
}

TEST_CASE ("Test TrapezoidalIntegration, x1=x2=0") {
  const double PI = 3.14159265358979323846;
  const double E = 2.71828182845904523536;

  auto f = [E, PI](double x) { 
    return 5.0/(std::pow(E, PI) - 2.0) * exp(2.0*x) * cos(x); 
  };

  const double result = NumLib::TrapezoidalIntegration(0, 0, 100, f);
  CHECK(TestUtils::IsValid(result, 0.0, 1e-3));
}

TEST_CASE ("Test TrapezoidalIntegration, N=0") {
  const double PI = 3.14159265358979323846;
  const double E = 2.71828182845904523536;

  auto f = [E, PI](double x) { 
    return 5.0/(std::pow(E, PI) - 2.0) * exp(2.0*x) * cos(x); 
  };

  const double result = NumLib::TrapezoidalIntegration(0, 0.5*PI, 0, f);
  CHECK(TestUtils::IsValid(result, 0.18575, 1e-3));
}

TEST_CASE ("Test TrapezoidalIntegration, f=nullptr") {
  const double PI = 3.14159265358979323846;

  const double result = NumLib::TrapezoidalIntegration(0, 0.5*PI, 0, nullptr);
  CHECK(TestUtils::IsValid(result, 0.0, 1e-3));
}


#endif
