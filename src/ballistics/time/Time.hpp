//
// Created by Арсений Плахотнюк on 15.09.2022.
//
#ifndef TIME_HPP
#define TIME_HPP
#include "ballistics/types/BasicTypes.hpp"
#include "third_party/sofa/sofa.h"
#include "cmath"
#include <utility>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <third_party/expected/Expected.h>
#include <third_party/g3log/g3log.hpp>
#include "Error.hpp"
#include "third_party/fast-cpp-csv-parser/csv.h"
#include "limits"
#include "ballistics/exception/ProjectException.hpp"
#include "ballistics/constants/TimeConstants.hpp"

namespace Ballistics::Time {
    enum class Scale {
        UT1 = 1,
        UTC = 2,
        TAI = 3,
        TT = 4,
        TCG = 5,
        TCB = 6,
        TDB = 7,

    };

    template<Scale scale>
    class Time {
        integer jdDay_;
        scalar jdPart_;

    public:

        /** The default constructor of a Time class
        */
        inline constexpr Time() noexcept = default;

        void correct_format(){
            while(jdPart_ < 0){
                jdPart_ += 1.;
                jdDay_ -= 1;
            }
            while(jdPart_ >= 1){
                jdPart_ -= 1.;
                jdDay_ += 1;
            }
        }

        inline constexpr  Time(const scalar day, const scalar part) {
            jdDay_ = static_cast<integer>(day);
            jdPart_ = day - std::floor(day) + part;
            while(jdPart_ >= 1.)
            {
                jdPart_ -= 1.;
                jdDay_ += 1;
            }
            while(jdPart_ < 0.)
            {
                jdPart_ += 1.;
                jdDay_ -= 1;
            }

        };

        inline constexpr  Time(const integer jdDay, const scalar jdPart) : jdDay_(jdDay), jdPart_(jdPart)
        {
            correct_format();
        }


        explicit Time(scalar jdDate) {
            double day;
            double part = modf(jdDate, &day);
            jdDay_ = static_cast<integer>(day);
            jdPart_ = static_cast<scalar>(part);
            correct_format();
        }

        [[nodiscard]] scalar toScalar() const noexcept{
            return jdDay_ + jdPart_;
        }

        /** Operator for adding seconds to Time class object
        * @param s number of seconds
        * @return new moment of time
        */
        [[nodiscard]] constexpr Time operator+(const scalar s) const noexcept
        {
            integer jdDay = jdDay_;
            scalar jdPart = jdPart_;
            return {jdDay, jdPart + s / SECONDS_IN_DAY};
        }

        /** Operator for subtracting seconds from Time class object
        * @param s number of seconds
        * @return new moment of time
        */
        [[nodiscard]] constexpr Time operator-(const scalar s) const noexcept
        {
            integer jdDay = jdDay_;
            scalar jdPart = jdPart_;
            return {jdDay, jdPart - s / SECONDS_IN_DAY};
        }

        /** Operator for subtracting one moment of time from another
        * @param r subtracted point in time
        * @return the difference between two points in time in seconds
        */
        [[nodiscard]] constexpr scalar operator-(const Time &r) const noexcept
        {
            integer jdDay = jdDay_;
            scalar jdPart = jdPart_;
            return ((jdDay - r.jdDay_) + (jdPart - r.jdPart_)) * SECONDS_IN_DAY;
        }

        [[nodiscard]] integer get_jdDay() const { return jdDay_; }

        [[nodiscard]] scalar get_jdPart() const { return jdPart_; }

        [[nodiscard]] scalar from_jd_to_mjd() const noexcept {
            return (jdDay_ - 2400000.5) + jdPart_;
        }

        [[nodiscard]] scalar from_mjd_to_jd() const noexcept {
            return (jdDay_ + 2400000.5) + jdPart_;
        }

        void add_n_sec(integer n) noexcept {
            double days;
            modf(n / 86400., &days);
            if (jdPart_ + n / 86400. > 1.) {
                jdDay_ += static_cast<integer>(days);
                jdPart_ = jdPart_ + n / 86400. - days;
            } else
                jdPart_ += n / 86400.;
            correct_format();
        }

        void subtract_n_sec(integer n) noexcept {
            double days;
            modf(n / 86400., &days);
            if (jdPart_ - n / 86400. < 0.) {
                jdDay_ -= static_cast<integer>(days);
                jdPart_ = jdPart_ - n / 86400.0 + days;
            } else
                jdPart_ -= n / 86400.0;
            correct_format();
        }

        Time operator<(const Time<scale> &time2) {
            if (jdDay_ < time2.jdDay_) { return true; }
            else if (jdDay_ == time2.jdDay_) { return jdPart_ < time2.jdPart_; }
            else { return false; }
        }

        bool operator==(const Time<scale> &time2) {
            if (jdDay_ == time2.jdDay_ && std::fabs(jdPart_ - time2.jdPart_) < 1.e-12){
                return true;
            }
            else if(jdDay_ - time2.jdDay_ == 1 &&
            time2.jdPart_ - 1 <= std::numeric_limits<double>::epsilon() &&
            jdPart_ <= std::numeric_limits<double>::epsilon()){
                return true;
            }
            else if(time2.jdDay_ - jdDay_ == 1 &&
                    jdPart_ - 1 <= std::numeric_limits<double>::epsilon() &&
                    time2.jdPart_ <= std::numeric_limits<double>::epsilon()){
                return true;
            }
            return false;
        }


        Time operator>(const Time<scale> &time2) {
            if (jdDay_ > time2.jdDay_) { return true; }
            else if (jdDay_ == time2.jdDay_) { return jdPart_ > time2.jdPart_; }
            else { return false; }
        }

    };

    /** Function for calculating number of days since the beginning of the year
     * @param t moment of time in scale UTC
     * @return number of days since the beginning of the year
     */
    [[nodiscard]] static  scalar calcDaysSinceJan1(const Time<Scale::UTC> &t)
    {
        int year, month, day;
        double fd;
        const int i = iauJd2cal(t.get_jdDay(), t.get_jdPart(), &year, &month, &day, &fd);
        if (i != 0)
        {
            std::stringstream buff;
            buff << "The return value of the iauJd2cal() function must be equal to 0. The received value: " << i << ". File: " << __FILE__ << ". Line: " << __LINE__;
            throw ProjectException(buff.str());
        }

        scalar t1, t2;
        const int j = iauDtf2d("UTC", year - 1, 12, 31, 0, 0, 0., &t1, &t2);
        if (j != 0)
        {
            std::stringstream buff;
            buff << "The return value of the iauDtf2d() function must be equal to 0. The received value: " << j << ". File: " << __FILE__ << ". Line: " << __LINE__;
            throw ProjectException(buff.str());
        }

        return static_cast<scalar>(t.get_jdDay()) + t.get_jdPart() - t1 - t2;
    }

}
#endif //TIME_HPP
