#ifndef INTEGRATION_H
#define INTEGRATION_H

#include <doctest/doctest.h>
#include <functional>
#include <vector>

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

  /**
   * @brief A numerical integration method using the simpson approach
   * 
   * SimpsonIntegration is using the simpson approach to acomplish numerical integrations as described at [thoughts-on-coding.com](https://thoughts-on-coding.com/2019/04/17/numerical-methods-in-cpp-part-1-newton-cotes-integration/)
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
   * I=\textstyle \sum_{0}^{N}\frac{x_{i+1}-x_i}{6}\left( f(x_i)+4f\left( \frac{x_1+x_2}{2} \right)+f(x_{i+1}) \right)
   * \f]
   * 
   * @image html simpson.png
   * 
   * @test Test SimpsonIntegration, valid case: Testing method with function \f$I=\int_0^{\pi/2} \frac{5}{e^\pi-2}\exp(2x)\cos(x)dx=1.0\f$
   * @test Test SimpsonIntegration, x1=x2=0: Testing method with zero range
   * @test Test SimpsonIntegration, N=0: Testing method with zero iteration steps
   * @test Test SimpsonIntegration, f=nullptr: Testing method with error in function f
   * 
   * @param x1 begining of Integration range \f$x_1=a\f$
   * @param x2 end of Integration range \f$x_2=b\f$
   * @param N number of iteration steps
   * @param f function \f$f\f$ to integrate
   * 
   * @return Integral of function \f$f\f$ in range \f$[a-b]\f$
   */
  double SimpsonIntegration(double x1, double x2, size_t N, const std::function<double (double)> &f);

  /**
   * @brief A numerical integration method using the romberg approach
   * 
   * RombergIntegration is using the simpson approach to acomplish numerical integrations as described at [thoughts-on-coding.com](https://thoughts-on-coding.com/2019/04/17/numerical-methods-in-cpp-part-1-newton-cotes-integration/)
   * and calculates the integral of function \f$f\f$ in range \f$[a-b]\f$.
   * 
   * \f[
   * R\left( 0,0 \right)=h_1\left( f\left( a \right) + f\left( b \right)\right)
   * \f]
   * \f[
   * R\left( n,0 \right)=\frac{1}{2}R\left( n-1,0 \right)+h_n\sum_{k=1}^{2^{\pi-1}}f\left( a+\left( 2k-1 \right)h_n \right)
   * \f]
   * \f[
   * R\left( n,m \right)=R\left( n,m-1 \right)+\frac{R\left( n,m-1 \right)-R\left( n-1,m-1 \right)}{4^m-1}
   * \f]
   * 
   * @image html romberg.png
   * 
   * @test Test RombergIntegration, valid case: Testing method with function \f$I=\int_0^{\pi/2} \frac{5}{e^\pi-2}\exp(2x)\cos(x)dx=1.0\f$
   * @test Test RombergIntegration, x1=x2=0: Testing method with zero range
   * @test Test RombergIntegration, N=0: Testing method with zero iteration steps
   * @test Test RombergIntegration, f=nullptr: Testing method with error in function f
   * 
   * @param x1 begining of Integration range \f$x_1=a\f$
   * @param x2 end of Integration range \f$x_2=b\f$
   * @param N number of iteration steps
   * @param f function \f$f\f$ to integrate
   * 
   * @return Integral of function \f$f\f$ in range \f$[a-b]\f$
   */
  std::vector<std::vector<double>> RombergIntegration(double x1, double x2, size_t N, const std::function<double (double)> &f);

  class GaussLegendreIntegration {
    public:
      double operator () (double x1, double x2, size_t N, const std::function<double (double)> &f) const;

    private:
      class LegendrePolynomial {
        public:
            LegendrePolynomial(double lowerBound, double upperBound, size_t numberOfIterations);

            const std::vector<double> & getWeight() const;

            const std::vector<double> & getRoot() const;

        private:
            const static double EPSILON;

            struct Result {
                double value;
                double derivative;

                Result() : value(0), derivative(0) {}
                Result(double val, double deriv) : value(val), derivative(deriv) {}
            };

            void calculateWeightAndRoot();

            Result calculatePolynomialValueAndDerivative(double x);

            const double mLowerBound;
            const double mUpperBound;
            const int mNumberOfIterations;
            std::vector<double> mWeight;
            std::vector<double> mRoot;
        };
  };
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

TEST_CASE ("Test SimpsonIntegration, valid case") {
  const double PI = 3.14159265358979323846;
  const double E = 2.71828182845904523536;

  auto f = [E, PI](double x) { 
    return 5.0/(std::pow(E, PI) - 2.0) * exp(2.0*x) * cos(x); 
  };

  const double result = NumLib::SimpsonIntegration(0, 0.5*PI, 100, f);
  CHECK(TestUtils::IsValid(result, 1.0, 1e-3));
}

TEST_CASE ("Test SimpsonIntegration, x1=x2=0") {
  const double PI = 3.14159265358979323846;
  const double E = 2.71828182845904523536;

  auto f = [E, PI](double x) { 
    return 5.0/(std::pow(E, PI) - 2.0) * exp(2.0*x) * cos(x); 
  };

  const double result = NumLib::SimpsonIntegration(0, 0, 100, f);
  CHECK(TestUtils::IsValid(result, 0.0, 1e-3));
}

TEST_CASE ("Test SimpsonIntegration, N=0") {
  const double PI = 3.14159265358979323846;
  const double E = 2.71828182845904523536;

  auto f = [E, PI](double x) { 
    return 5.0/(std::pow(E, PI) - 2.0) * exp(2.0*x) * cos(x); 
  };

  const double result = NumLib::SimpsonIntegration(0, 0.5*PI, 0, f);
  CHECK(TestUtils::IsValid(result, 0.905, 1e-3));
}

TEST_CASE ("Test SimpsonIntegration, f=nullptr") {
  const double PI = 3.14159265358979323846;

  const double result = NumLib::SimpsonIntegration(0, 0.5*PI, 0, nullptr);
  CHECK(TestUtils::IsValid(result, 0.0, 1e-3));
}

TEST_CASE ("Test RombergIntegration, valid case") {
  const double PI = 3.14159265358979323846;
  const double E = 2.71828182845904523536;

  auto f = [E, PI](double x) { 
    return 5.0/(std::pow(E, PI) - 2.0) * exp(2.0*x) * cos(x); 
  };

  const double result = NumLib::RombergIntegration(0, 0.5*PI, 4, f).back().back();
  CHECK(TestUtils::IsValid(result, 1.0, 1e-5));
}

TEST_CASE ("Test RombergIntegration, x1=x2=0") {
  const double PI = 3.14159265358979323846;
  const double E = 2.71828182845904523536;

  auto f = [E, PI](double x) { 
    return 5.0/(std::pow(E, PI) - 2.0) * exp(2.0*x) * cos(x); 
  };

  const double result = NumLib::RombergIntegration(0, 0, 4, f).back().back();
  CHECK(TestUtils::IsValid(result, 0.0, 1e-5));
}

TEST_CASE ("Test RombergIntegration, N=0") {
  const double PI = 3.14159265358979323846;
  const double E = 2.71828182845904523536;

  auto f = [E, PI](double x) { 
    return 5.0/(std::pow(E, PI) - 2.0) * exp(2.0*x) * cos(x); 
  };

  const double result = NumLib::RombergIntegration(0, 0.5*PI, 0, f).back().back();
  CHECK(TestUtils::IsValid(result, 0.18575, 1e-3));
}

TEST_CASE ("Test RombergIntegration, f=nullptr") {
  const double PI = 3.14159265358979323846;

  const double result = NumLib::RombergIntegration(0, 0.5*PI, 0, nullptr).back().back();
  CHECK(TestUtils::IsValid(result, 0.0, 1e-3));
}

#endif
