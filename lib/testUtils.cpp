#include "testUtils.h"

bool TestUtils::IsValid(double result, double expected, double epsilon) {
    static const double CLOSEST_TO_ZERO = std::numeric_limits<double>::min();

    if(fabs(result) >= CLOSEST_TO_ZERO && fabs(expected) >= CLOSEST_TO_ZERO) {
        return fabs(result/expected - 1) <= epsilon;
    }
    if(fabs(result) < CLOSEST_TO_ZERO && fabs(expected) < CLOSEST_TO_ZERO) {
        return true;
    }

    return false;
}