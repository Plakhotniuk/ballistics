#ifndef PROJECT_BALLISTICS_ELLIPSOID_HPP
#define PROJECT_BALLISTICS_ELLIPSOID_HPP

#include "ballistics/types/BasicTypes.hpp"

namespace Ballistics::Ellipsoid
{
    class Ellipsoid
    {
    private:
        scalar a_;  // large semi-axis of the ellipsoid
        scalar f_;  // flattening of the ellipsoid
        scalar e2_;
        scalar e2m_;
        scalar e2a_;
        scalar e4a_;
        scalar maxRad_;
        scalar polarRadius_;
    public:
        /** Constructor of Ellipsoid class
         * @param a the large semi-axis of the ellipsoid
         * @param f flattening of the ellipsoid
         * @return Ellipsoid class object
         */
        constexpr inline Ellipsoid(scalar a,
                                   scalar f)
                : a_(a),
                  f_(f),
                  e2_(f * (2 - f)),
                  e2m_(1 - e2_),
                  e2a_(e2_),
                  e4a_(e2_ * e2_),
                  maxRad_(2 * a / std::numeric_limits<scalar>::epsilon()),
                  polarRadius_(a * (1 - f)) {}

        /** Method for calculating the height above the ellipsoid by geocentric coordinates
         * @param position geocentric coordinates in meters
         * @return height above the ellipsoid in meters
         */
        [[nodiscard]] scalar calcHeight(const Vector3d &position) const noexcept;

        [[nodiscard]] scalar getLargeSemiAxis() const noexcept { return a_; }

        [[nodiscard]] scalar getFlattening() const noexcept { return f_; }
    };

    static constexpr Ellipsoid WGS84Ellipsoid(6378137.,
                                              1. / (static_cast<scalar>(298257223563LL) / 1000000000.));
}

#endif //PROJECT_BALLISTICS_ELLIPSOID_HPP