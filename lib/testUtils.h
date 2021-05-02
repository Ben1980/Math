#ifndef TESTUTILS_H
#define TESTUTILS_H

#include <cmath>
#include <limits>

namespace TestUtils {
    bool IsValid(double result, double expected, double epsilon);
};

#endif