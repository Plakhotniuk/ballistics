#ifndef PROJECT_BALLISTICS_UTILITY_HPP
#define PROJECT_BALLISTICS_UTILITY_HPP

#include "../types/BasicTypes.hpp"
#include <vector>
#include <array>
#include <cmath>

namespace Ballistics
{
    /** Function for calculating the partial sum of a power series
     * SUM = a0 + a1 * x + a2 * x^2 + a3 * x^3 + ...
     * @param coeffs coefficients of the power series
     * @param var variable of the power series
     * @return partial sum of a power series
     */
    template<typename T>
    scalar calcPartialSumPowerSeries(const std::vector<T> &coeffs,
                                     const T var)
    {
        T res = 0;
        T variableInPowerI = 1;
        for (int i = 0; i < coeffs.size(); ++i)
        {
            res += coeffs[i] * variableInPowerI;
            variableInPowerI *= var;
        }
        return res;
    }
}

#endif //PROJECT_BALLISTICS_UTILITY_HPP
