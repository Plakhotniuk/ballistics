#ifndef PROJECT_BALLISTICS_SHADOW_HPP
#define PROJECT_BALLISTICS_SHADOW_HPP

#include "ballistics/types/BasicTypes.hpp"

namespace Ballistics::Shadow
{
    /** Function for calculating the shadow function
     * @param r position vector of the satellite in the geocentric coordinate system
     * @param rs position vector of the Sun in the geocentric coordinate system
     * @param solarRadius radius of the Sun
     * @param earthRadius radius of the Earth
     * @return value of a shadow function
     */
    [[nodiscard]] scalar calcShadowFunction(const Vector3d &r,
                                            const Vector3d &rs,
                                            const scalar solarRadius,
                                            const scalar earthRadius)
    {
        const scalar normR = r.norm();
        const scalar normDeltaRsR = (rs - r).norm();

        const scalar thetaE = asin(earthRadius / normR);
        const scalar thetaS = asin(solarRadius / normDeltaRsR);
        const scalar thetaES = acos((-r.transpose().dot(rs - r)) / (normR * normDeltaRsR));

        scalar S = 0;

        if (thetaE - thetaS > thetaES)
        {
            S = M_PI * thetaS * thetaS;
        } else
        {
            if (thetaS - thetaE > thetaES)
            {
                S = M_PI * thetaE * thetaE;
            } else
            {
                if (abs(thetaE - thetaS) < thetaES && thetaES < thetaE + thetaS)
                {
                    const scalar thetaSSqr = thetaS * thetaS;
                    const scalar thetaESSqr = thetaES * thetaES;
                    const scalar thetaESqr = thetaE * thetaE;

                    const scalar angleCAF = acos((thetaSSqr + thetaESSqr - thetaESqr) / (2 * thetaS * thetaES));
                    const scalar angleCBD = acos((thetaESqr + thetaESSqr - thetaSSqr) / (2 * thetaE * thetaES));

                    const scalar areaAFC = (1. / 2.) * angleCAF * thetaSSqr;
                    const scalar areaAEC = (1. / 2.) * (thetaS * sin(angleCAF)) * (thetaS * cos(angleCAF));
                    const scalar areaBDC = (1. / 2.) * angleCBD * thetaESqr;
                    const scalar areaBEC = (1. / 2.) * (thetaE * sin(angleCBD)) * (thetaE * cos(angleCBD));

                    S = 2 * (areaAFC - areaAEC) + 2 * (areaBDC - areaBEC);
                }
            }
        }

        return 1 - S / (M_PI * thetaS * thetaS);
    }
}

#endif //PROJECT_BALLISTICS_SHADOW_HPP
