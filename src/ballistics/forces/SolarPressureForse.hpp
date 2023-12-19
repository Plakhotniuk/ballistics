#ifndef PROJECT_BALLISTICS_SOLARPRESSUREFORSE_HPP
#define PROJECT_BALLISTICS_SOLARPRESSUREFORSE_HPP

#include "ballistics/types/BasicTypes.hpp"
#include "ballistics/constants/FundamentalConstants.hpp"

namespace Ballistics::Forces
{

     struct SolarPressureForce{

         [[nodiscard]] static Vector3d calcForce(const Vector3d& r, const scalar F, const scalar S,
                                                 const scalar TSI, const scalar AU) noexcept
         {
             /** Function for calculating the force of solar pressure
              * @param r position vector of the satellite in the heliocentric coordinate system
              * @param F value of the shadow function
              * @param S effective area
              * @param TSI total solar irradiance
              * @param AU astronomical unit in metres
              * @return force of solar pressure
              */

             const scalar rNormSqr = r.squaredNorm();
             const scalar rNorm = std::sqrt(rNormSqr);
             return TSI * (AU * AU) / (rNorm * rNormSqr * c) * F * S * r;
         }
     };

}

#endif //PROJECT_BALLISTICS_SOLARPRESSUREFORSE_HPP
