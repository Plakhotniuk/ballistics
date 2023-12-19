#ifndef BALLISTICS_ATMOSPHEREMODEL_HPP
#define BALLISTICS_ATMOSPHEREMODEL_HPP

#include "../time/Time.hpp"
#include "../types/BasicTypes.hpp"
#include "../constants/AtmosphereConstants.hpp"
#include "Utility.hpp"
#include <cmath>

namespace Ballistics::Atmosphere
{
    class AtmosphereModel
    {
    public:
        /** Method for calculating density of atmosphere in model GOST R 25645.166-2004
         * (link: https://docs.cntd.ru/document/1200036026)
         * @param r coordinate vector in the system ITRS
         * @param t time in scale UTC
         * @param h height above the Earth ellipsoid in kilometres
         * @param F10_7 average daily index of solar activity
         * @param F81 averaged F10_7 index of solar activity for three solar revolutions
         * @param F0 index F81 rounded to the nearest number multiple of 25
         * @param kpp modified values of 3-hour geomagnetic disturbance indexes
         * @param alpha right ascension
         * @param delta declination
         * @param S Greenwich mean sidereal time
         * @param d number of days from the beginning of the year
         * @return density of atmosphere
         */
        static scalar calcDensity(const Eigen::Vector3d &r,
                                  const Time::Time<Time::Scale::UTC> &t,
                                  const scalar h,
                                  const scalar F10_7,
                                  const scalar F81,
                                  const integer F0,
                                  const scalar kpp,
                                  const scalar alpha,
                                  const scalar delta,
                                  const scalar S,
                                  const scalar d);
    };
}

#endif //BALLISTICS_ATMOSPHEREMODEL_HPP