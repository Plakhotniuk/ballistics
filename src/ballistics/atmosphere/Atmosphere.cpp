#include "Atmosphere.hpp"

namespace Ballistics::Atmosphere
{
    void Atmosphere::initialize()
    {
        const rapidcsv::Document doc(csvLocation_);
        // Reading column "DATE" from the file
        const std::vector<std::string> jdStringCol = doc.GetColumn<std::string>("DATE");
        // Converting the date from the calendar to jd format and filling jdCol
        for (int i = 0; i < jdStringCol.size(); ++i)
        {
            const std::string year = {jdStringCol[i][0], jdStringCol[i][1], jdStringCol[i][2], jdStringCol[i][3]};
            const std::string month = {jdStringCol[i][5], jdStringCol[i][6]};
            const std::string day = {jdStringCol[i][8], jdStringCol[i][9]};
            scalar Day, Part;
            const std::string Scale = "UTC";
            int j = iauDtf2d(Scale.c_str(), atoi(year.c_str()), atoi(month.c_str()), atoi(day.c_str()),
                             0, 0, 0., &Day, &Part);
            if (j != 0)
            {
                std::stringstream buff;
                buff << "The return value of the iauDtf2d() function must be equal to 0. The received value: " << j << ". File: " << __FILE__ << ". Line: " << __LINE__;
                throw ProjectException(buff.str());
            }

            jdCol_.push_back(Day + Part);
        }
        // Reading columns from the file
        F10_7Col_ = doc.GetColumn<double>("F10.7_OBS");
        std::array<std::vector<scalar>, 8> kpArr;
        kpArr[0] = doc.GetColumn<scalar>("KP1");
        kpArr[1] = doc.GetColumn<scalar>("KP2");
        kpArr[2] = doc.GetColumn<scalar>("KP3");
        kpArr[3] = doc.GetColumn<scalar>("KP4");
        kpArr[4] = doc.GetColumn<scalar>("KP5");
        kpArr[5] = doc.GetColumn<scalar>("KP6");
        kpArr[6] = doc.GetColumn<scalar>("KP7");
        kpArr[7] = doc.GetColumn<scalar>("KP8");

        // Filling first 80 element of F81Col_ with zeros
        for (int i = 0; i <= 79; ++i)
        {
            F81Col_.push_back(0);
        }

        // Calculating F81[i] (where i starts with 80) and filling F81Col_
        for (int i = 80; i < F10_7Col_.size(); ++i)
        {
            scalar W = 0;
            scalar numerator = 0;
            scalar denominator = 0;
            for (int j = -80; j <= 0; ++j)
            {
                W = 1 + (0.5 * j) / 80;
                numerator += F10_7Col_[i + j] * W;
                denominator += W;
            }
            F81Col_.push_back(numerator / denominator);
        }

        const integer m = kpArr.size();       // Size of vector of arrays
        const integer n = kpArr[0].size();    // Size of each array
        // Converting kp to the right format
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                kpArr[i][j] = calcKp(static_cast<integer>(kpArr[i][j]));
            }
        }

        // Filling in initial value (kpp[0] = kp[0])
        kppArr_[0].push_back(kpArr[0][0]);

        // Calculating of kpp and filling in kppArr_
        for (int j = 0; j < n; ++j)
        {
            for (int i = 0; i < m; ++i)
            {
                scalar delta;
                if (i == 0)
                {
                    if (j == 0)
                    {
                        continue;
                    }
                    delta = kpArr[i][j] - kppArr_[m - 1][j - 1];
                } else
                {
                    delta = kpArr[i][j] - kppArr_[i - 1][j];
                }
                scalar r;
                delta > 0 ? r = 0.3 : r = 0.7;
                kppArr_[i].push_back(kpArr[i][j] - r * delta);
            }
        }
    }

    [[nodiscard]] scalar Atmosphere::calcKp(const integer kp) const noexcept
    {
        scalar val = 0.;
        if (kp % 10 == 0)
        {
            val += static_cast<scalar>(kp / 10) + 0.0000;
        } else
        {
            if (kp % 10 == 3)
            {
                val += static_cast<scalar>(kp / 10) + 0.3333;
            } else
            {
                if (kp % 10 == 7)
                {
                    val += static_cast<scalar>(kp / 10) + 0.6667;
                }
            }
        }
        return val;
    }

    [[nodiscard]] Indexes Atmosphere::interpolate(const scalar jd) const noexcept
    {
        Indexes indexes{};
        const scalar jd3hTimeLag = jd - timeObsKpThreeHours;
        const scalar jd107TimeLag = jd - F107TimeLag - timeObsF10_7;

        for(int i = 0; i < jdCol_.size(); ++i)
        {
            if(jdCol_[i] <= jd107TimeLag && jdCol_[i + 1] > jd107TimeLag)
            {
                indexes.F10_7 = ((jd107TimeLag - jdCol_[i]) / (jdCol_[i + 1] - jdCol_[i])) *
                                (F10_7Col_[i + 1] - F10_7Col_[i]) + F10_7Col_[i];
                indexes.F81 = ((jd107TimeLag - jdCol_[i]) / (jdCol_[i + 1] - jdCol_[i])) *
                              (F81Col_[i + 1] - F81Col_[i]) + F81Col_[i];
            }

            if(jdCol_[i] <= jd3hTimeLag && jdCol_[i + 1] > jd3hTimeLag)
            {
                // k is an array number in vector kppArr_
                const integer k = calcNumEighthPartOfDay(jd3hTimeLag);
                if (k == 7)
                {
                    indexes.kpp = ((jd3hTimeLag - jdCol_[i]) / (jdCol_[i + 1] - jdCol_[i])) *
                                  (kppArr_[0][i + 1] - kppArr_[k][i]) + kppArr_[k][i];
                } else
                {
                    indexes.kpp = ((jd3hTimeLag - jdCol_[i]) / (jdCol_[i + 1] - jdCol_[i])) *
                                  (kppArr_[k + 1][i] - kppArr_[k][i]) + kppArr_[k][i];
                }
            }
        }
        return indexes;
    }

    [[nodiscard]] integer Atmosphere::calcNumEighthPartOfDay(const scalar jd) const
    {
        const Time::Time<Time::Scale::UTC> timeUTC(jd, 0.);
        int year, month, day;
        double fd;
        const int i = iauJd2cal(timeUTC.get_jdDay(), timeUTC.get_jdPart(), &year, &month, &day, &fd);
        if (i != 0)
        {
            std::stringstream buff;
            buff << "The return value of the iauJd2cal() function must be equal to 0. The received value: " << i << ". File: " << __FILE__ << ". Line: " << __LINE__;
            throw ProjectException(buff.str());
        }

        scalar t1, t2;
        const int j = iauDtf2d("UTC", year, month, day, 0, 0, 0., &t1, &t2);
        if (j != 0)
        {
            std::stringstream buff;
            buff << "The return value of the iauDtf2d() function must be equal to 0. The received value: " << j << ". File: " << __FILE__ << ". Line: " << __LINE__;
            throw ProjectException(buff.str());
        }

        const scalar hours = (jd - t1 - t2) * HOURS_IN_DAY;
        return floor(hours / 3);
    }

    [[nodiscard]] integer Atmosphere::calcF0(const scalar F81) const noexcept
    {
        const integer F81Int = round(F81);
        const integer lowBorder = (F81Int / 25) * 25;
        const integer upperBorder = lowBorder + 25;
        if (abs(F81Int - lowBorder) >= abs(F81Int - upperBorder))
        {
            return upperBorder;
        } else
        {
            return lowBorder;
        }
    }

    [[nodiscard]] scalar Atmosphere::calcGMST(Time::Time<Time::Scale::UTC> timeUTC) const
    {
        int year, month, day;
        double fd;
        const int i = iauJd2cal(timeUTC.get_jdDay(), timeUTC.get_jdPart(), &year, &month, &day, &fd);
        if (i != 0)
        {
            std::stringstream buff;
            buff << "The return value of the iauJd2cal() function must be equal to 0. The received value: " << i << ". File: " << __FILE__ << ". Line: " << __LINE__;
            throw ProjectException(buff.str());
        }

        scalar t1, t2;
        const int j = iauDtf2d("UTC", year, month, day, 0, 0, 0., &t1, &t2);
        if (j != 0)
        {
            std::stringstream buff;
            buff << "The return value of the iauDtf2d() function must be equal to 0. The received value: " << j << ". File: " << __FILE__ << ". Line: " << __LINE__;
            throw ProjectException(buff.str());
        }
        Ballistics::Time::Time<Ballistics::Time::Scale::UTC> time(t1, t2);
        const Expected<Ballistics::Time::Time<Ballistics::Time::Scale::UT1>, Error> timeUT1 = converter_.convert<Ballistics::Time::Scale::UT1, Ballistics::Time::Scale::UTC>(time);
        const Expected<Ballistics::Time::Time<Ballistics::Time::Scale::TT>, Error> timeTT = converter_.convert<Time::Scale::TT, Time::Scale::UTC>(time);
        return iauGmst06(timeUT1.value().get_jdDay(), timeUT1.value().get_jdPart(), timeTT.value().get_jdDay(), timeTT.value().get_jdPart());
    }

    [[nodiscard]] Vector3d Atmosphere::calcSunPosition(const Time::Time<Time::Scale::UTC> &t) const
    {
        const Expected<Ballistics::Time::Time<Ballistics::Time::Scale::TDB>, Error> timeTDB = converter_.convert<Time::Scale::TDB, Time::Scale::UTC>(t);
        return ephemeris_.calcCoordinates(timeTDB.value(), 11, 3);
    }

    [[nodiscard]] scalar Atmosphere::calcRightAscension(const Vector3d &pos) const noexcept
    {
        scalar rightAscension = atan(pos.y() / pos.x());
        if(rightAscension < 0.)
        {
            rightAscension += 2 * M_PI;
        }
        return rightAscension;
    }

    [[nodiscard]] scalar Atmosphere::calcDeclination(const Vector3d &pos) const noexcept
    {
        const scalar declination = asin(pos.z() / sqrt(pos.x() * pos.x() + pos.y() * pos.y() + pos.z() * pos.z()));
        return declination;
    }

    [[nodiscard]] scalar Atmosphere::calcHeight(const Vector3d &position) const noexcept
    {
        return EarthEllipsoid_.calcHeight(position);
    }

    [[nodiscard]] scalar Atmosphere::getDensity(const Vector3d &r,
                                  const Time::Time<Time::Scale::UTC> &t) const
    {
        const Indexes indexes = interpolate(static_cast<scalar>(t.get_jdDay()) + t.get_jdPart());
        const Vector3d position = calcSunPosition(t);
        return AtmosphereModel::calcDensity(r,
                                            t,
                                            calcHeight(r) / 1000,
                                            indexes.F10_7,
                                            indexes.F81,
                                            calcF0(indexes.F81),
                                            indexes.kpp,
                                            calcRightAscension(position),
                                            calcDeclination(position),
                                            calcGMST(t),
                                            calcDaysSinceJan1(t));
    }

}
