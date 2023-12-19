//
// Created by Арсений Плахотнюк on 24.11.2022.
//

#ifndef PROJECT_BALLISTICS_EPHEMERIS_HPP
#define PROJECT_BALLISTICS_EPHEMERIS_HPP

#include "../../../third_party/calceph/calceph.h"
#include "../exception/ProjectException.hpp"
#include "../time/Time.hpp"
#include <string>

namespace Ballistics::Ephemeris
{
    enum Order
    {
        ONLY_POS         = 0,
        ONLY_POS_AND_VEL = 1
    };

    enum Body
    {
        MercuryBarycenter     = 1,
        VenusBarycenter       = 2,
        Earth                 = 3,
        MarsBarycenter        = 4,
        JupiterBarycenter     = 5,
        SaturnBarycenter      = 6,
        UranusBarycenter      = 7,
        NeptuneBarycenter     = 8,
        PlutoBarycenter       = 9,
        Moon                  = 10,
        Sun                   = 11,
        SolarSystemBarycenter = 12,
        EarthMoonBarycenter   = 13
    };


    class Ephemeris
    {
        t_calcephbin *eph_;
    public:
        /** Constructor of a Ephemeris
         * @param fileLocation location of the binary file with information about ephemeris
         */
        explicit Ephemeris(const std::string fileLocation = "../../resources/ephemeris/de405.bin")
        {
            eph_ = calceph_open(fileLocation.c_str());
            if (eph_ == NULL)
            {
                std::stringstream buff;
                buff << "Error opening the file: " << fileLocation << ". File: " << __FILE__ << ". Line: " << __LINE__;
                throw ProjectException(buff.str());
            }
        }

        /** Destructor of a Ephemeris
         */
        ~Ephemeris() { calceph_close(eph_); }

        /** Method for calculating the coordinates of the target body relative to the central one
         * at a certain point in time
         * @param jdDate date in jd format (possible scales TCB or TDB only)
         * @param target the body whose coordinates we want to get
         * @param center the body around which the target body rotates
         * @return coordinates of the target body in kilometres
         */
        template<Time::Scale T>
        [[nodiscard]] Vector3d calcCoordinates(const Time::Time<T> &jdDate,
                                               const integer target,
                                               const integer center)
        const
        {
            if (T != Time::Scale::TCB && T != Time::Scale::TDB)
            {
                std::stringstream buff;
                buff << "The jdDay argument can only be in TDB or TCB scales." << ". File: " << __FILE__ << ". Line: " << __LINE__;
                throw ProjectException(buff.str());
            }

            integer e = calceph_gettimescale(eph_);
            if (e == 0)
            {
                std::stringstream buff;
                buff << "An error occurred in the function calceph_gettimescale" << ". File: " << __FILE__ << ". Line: " << __LINE__;
                throw ProjectException(buff.str());
            } else
            {
                if (e == 1 && T != Time::Scale::TDB)
                {
                    std::stringstream buff;
                    buff << "Received timescale does not match with the timescale of the ephemeris file" << ". File: " << __FILE__ << ". Line: " << __LINE__;
                    throw ProjectException(buff.str());
                }

                if (e == 2 && T != Time::Scale::TCB)
                {
                    std::stringstream buff;
                    buff << "Received timescale does not match with the timescale of the ephemeris file" << ". File: " << __FILE__ << ". Line: " << __LINE__;
                    throw ProjectException(buff.str());
                }
            }

            e = calceph_prefetch(eph_);
            if (e == 0)
            {
                std::stringstream buff;
                buff << "An error occurred in the function calceph_prefetch" << ". File: " << __FILE__ << ". Line: " << __LINE__;
                throw ProjectException(buff.str());
            }
            scalar positions[3];
            e = calceph_compute_order(eph_,
                                      jdDate.get_jdDay(),
                                      jdDate.get_jdPart(),
                                      target,
                                      center,
                                      CALCEPH_UNIT_KM + CALCEPH_UNIT_SEC,
                                      ONLY_POS,
                                      positions);
            if (e == 0)
            {
                std::stringstream buff;
                buff << "An error occurred in the function calceph_compute_order" << ". File: " << __FILE__ << ". Line: " << __LINE__;
                throw ProjectException(buff.str());
            }

            Eigen::Vector3<scalar> vector = {positions[0], positions[1], positions[2]};
            return vector;
        }

        /** Method for calculating the velocities of the target body relative to the central one
        * at a certain point in time
        * @param jdDate date in jd format (possible scales TCB or TDB only)
        * @param target the body whose velocities we want to get
        * @param center the body around which the target body rotates
        * @return velocities of the target body in kilometres per second
        */
        template<Time::Scale T>
        [[nodiscard]] Vector3d calcVelocities(const Time::Time<T> &timeTDB,
                                              const integer target,
                                              const integer center)
        const
        {
            if (T != Time::Scale::TCB && T != Time::Scale::TDB)
            {
                std::stringstream buff;
                buff << "The jdDay argument can only be in TDB or TCB scales." << ". File: " << __FILE__ << ". Line: " << __LINE__;
                throw ProjectException(buff.str());
            }

            integer e = calceph_gettimescale(eph_);
            if (e == 0)
            {
                std::stringstream buff;
                buff << "An error occurred in the function calceph_gettimescale" << ". File: " << __FILE__ << ". Line: " << __LINE__;
                throw ProjectException(buff.str());
            } else
            {
                if (e == 1 && T != Time::Scale::TDB)
                {
                    std::stringstream buff;
                    buff << "Received timescale does not match with the timescale of the ephemeris file" << ". File: " << __FILE__ << ". Line: " << __LINE__;
                    throw ProjectException(buff.str());
                }

                if (e == 2 && T != Time::Scale::TCB)
                {
                    std::stringstream buff;
                    buff << "Received timescale does not match with the timescale of the ephemeris file" << ". File: " << __FILE__ << ". Line: " << __LINE__;
                    throw ProjectException(buff.str());
                }
            }

            e = calceph_prefetch(eph_);
            if (e == 0)
            {
                std::stringstream buff;
                buff << "An error occurred in the function calceph_prefetch" << ". File: " << __FILE__ << ". Line: " << __LINE__;
                throw ProjectException(buff.str());
            }

            scalar posAndVel[6];
            e = calceph_compute_order(eph_,
                                      timeTDB.get_jdDay(),
                                      timeTDB.get_jdPart(),
                                      target,
                                      center,
                                      CALCEPH_UNIT_KM + CALCEPH_UNIT_SEC,
                                      ONLY_POS_AND_VEL,
                                      posAndVel);
            if (e == 0)
            {
                std::stringstream buff;
                buff << "An error occurred in the function calceph_compute_order" << ". File: " << __FILE__ << ". Line: " << __LINE__;
                throw ProjectException(buff.str());
            }
            scalar velocities[3];
            for (int i = 0; i < 3; ++i)
            {
                velocities[i] = posAndVel[i + 3];
            }

            Eigen::Vector3<scalar> vector = {velocities[0], velocities[1], velocities[2]};
            return vector;
        }
    };
}

#endif //PROJECT_BALLISTICS_EPHEMERIS_HPP