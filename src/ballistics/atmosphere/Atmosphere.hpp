#ifndef BALLISTICS_ATMOSPHERE_HPP
#define BALLISTICS_ATMOSPHERE_HPP

#include "AtmosphereModel.hpp"
#include "ballistics/ellipsoid/Ellipsoid.hpp"
#include "ballistics/ephemeris/Ephemeris.hpp"
#include "../../../third_party/include/rapidcsv.h"
#include "ballistics/time/TimeConverter.hpp"
#include <cmath>

namespace Ballistics::Atmosphere
{
    struct Indexes
    {
        scalar F10_7;
        scalar F81;
        scalar kpp;
    };

    class Atmosphere
    {
    private:
        std::string csvLocation_;
        Ephemeris::Ephemeris ephemeris_;
        Ballistics::Time::TimeConverter converter_;
        Ellipsoid::Ellipsoid EarthEllipsoid_;

        std::vector<scalar> jdCol_;
        std::vector<scalar> F10_7Col_;
        std::vector<scalar> F81Col_;
        std::array<std::vector<scalar>, 8> kppArr_;

        /** Method for parsing scv file with solar weather and initializing private variables
         */
        void initialize();

        /** Method for converting kp indexes to the correct format
         * @param kp 3-hour geomagnetic disturbance index
         * @return kp index in correct format
         */
        [[nodiscard]] scalar calcKp(const integer kp) const noexcept;

        /** Method for interpolating indexes F10_7, F81 and kpp
         * @param jd moment of time in scale UTC
         * @return indexes F10_7, F81 and kpp
         */
        [[nodiscard]] Indexes interpolate(const scalar jd) const noexcept;

        /** Method for calculating the number of the eighth part of the day at a given time
         * @param jd moment of time in scale UTC
         * @return number of the eighth part of the day
         */
        [[nodiscard]] integer calcNumEighthPartOfDay(const scalar jd) const;

        /** Method for calculating index F0 by the index F81
         * index F0 is index F81 rounded to the nearest number multiple of 25
         * @param F81 averaged F10_7 index of solar activity for three solar revolutions
         * @return index F0
         */
        [[nodiscard]] integer calcF0(const scalar F81) const noexcept;

        /** Method for calculating Greenwich mean sidereal time
         * @param timeUTC moment of time in scale UTC
         * @return Greenwich mean sidereal time
         */
        [[nodiscard]] scalar calcGMST(Time::Time<Time::Scale::UTC> timeUTC) const;

        /** Method for calculating position of the Sun in the system ITRS
         * @param t moment of time in scale UTC
         * @return coordinates of Sun
         */
        [[nodiscard]] Vector3d calcSunPosition(const Time::Time<Time::Scale::UTC> &t) const;

        /** Method for calculating right ascension
         * @param pos coordinates of the Sun in the system ITRS
         * @return right ascension
         */
        [[nodiscard]] scalar calcRightAscension(const Vector3d &pos) const noexcept;

        /** Method for calculating declination
         * @param pos coordinates of the Sun in the system ITRS
         * @return declination
         */
        [[nodiscard]] scalar calcDeclination(const Vector3d &pos) const noexcept;

        /** Method for calculating height above the Earth ellipsoid in metres
         * @param position coordinates of the point in the system ITRS
         * @return height above the Earth ellipsoid
         */
        [[nodiscard]] scalar calcHeight(const Vector3d &position) const noexcept;
    public:
        /** Constructor of Atmosphere class
         * @param csvLocation location of scv file with solar weather data
         * @param a the large semi-axis of the Earth ellipsoid
         * (default value of a is Earth equatorial radius in system WGS 84)
         * @param f flattening of the Earth ellipsoid
         * (default value of f is Earth flattening in system WGS 84)
         * @return Atmosphere class object
         */
        explicit Atmosphere(const std::string &csvLocation,
                            const scalar a = Ellipsoid::WGS84Ellipsoid.getLargeSemiAxis(),
                            const scalar f = Ellipsoid::WGS84Ellipsoid.getFlattening())
                            : csvLocation_(csvLocation),
                              ephemeris_("../../resources/ephemeris/de405.bin"),
                              converter_("../../resources/rotation/earth_rotation.csv"),
                              EarthEllipsoid_(a, f)
        { initialize(); }

        /** Method for calculating density of atmosphere in model GOST R 25645.166-2004
         * (link: https://docs.cntd.ru/document/1200036026)
         * @param r coordinates of the point in the system ITRS
         * @param t moment of time in scale UTC
         * @return density of atmosphere in given point
         */
        [[nodiscard]] scalar getDensity(const Vector3d &r,
                                        const Time::Time<Time::Scale::UTC> &t) const;
    };
}

#endif //BALLISTICS_ATMOSPHERE_HPP