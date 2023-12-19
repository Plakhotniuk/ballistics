//
// Created by Арсений Плахотнюк on 06.10.2022.
//

#ifndef BALLISTICS_CONVERTER_H
#define BALLISTICS_CONVERTER_H

#include "ballistics/time/Time.hpp"
#include "third_party/sofa/sofa.h"
#include "third_party/sofa/sofam.h"

using namespace Ballistics::Time;

namespace Ballistics::Frame {

    enum Frame{
        ITRS = 0,
        GCRS = 1
    };

    struct Pole
    {
        double xp;
        double yp;
    };


    class FrameConverter {
    private:
        std::vector<int> mjd_;
        std::vector<double> x_arcsec_;
        std::vector<double> y_arcsec_;
    public:
        inline FrameConverter(std::string const & csvFileLocation){
            io::CSVReader<3> data(csvFileLocation);
            int mjd;
            double x, y;
            data.read_header(io::ignore_extra_column, "mjd", "x arcsec", "y arcsec");
            while (data.read_row(mjd, x, y)) {
                mjd_.emplace_back(mjd);
                x_arcsec_.emplace_back(x);
                y_arcsec_.emplace_back(y);
            }
        }

        [[nodiscard]] Pole interpolate(const double mjd) const noexcept
        {
            double xp = 0.;
            double yp = 0.;
            for(int i = 0; i < mjd_.size(); ++i)
            {
                if(mjd_[i] <= mjd && mjd_[i + 1] > mjd)
                {
                    xp = ((mjd - mjd_[i]) / (mjd_[i + 1] - mjd_[i])) * (x_arcsec_[i + 1] - x_arcsec_[i]) + x_arcsec_[i];
                    yp = ((mjd - mjd_[i]) / (mjd_[i + 1] - mjd_[i])) * (y_arcsec_[i + 1] - y_arcsec_[i]) + y_arcsec_[i];
                }
            }
            return {xp * DAS2R, yp * DAS2R};
        }

        [[nodiscard]] Pole getPole(const scalar mjd) const noexcept
        {
            return interpolate(mjd);
        }

        template<Frame From, Frame To>
        [[nodiscard]] Eigen::Quaterniond calcQuaternion(const Time::Time<Time::Scale::TT> &TT,
                                                        const Time::Time<Scale::UT1> &UT1,
                                                        const Time::Time<Time::Scale::UTC> &UTC) const noexcept;
    };

    template<>
    inline Eigen::Quaterniond FrameConverter::calcQuaternion<Frame::GCRS, Frame::ITRS>(const Time::Time<Time::Scale::TT> &TT,
                                                                                const Time::Time<Time::Scale::UT1> &UT1,
                                                                                const Time::Time<Time::Scale::UTC> &UTC)
    const noexcept
    {
        double x, y, s;
        iauXys06a(TT.get_jdDay(), TT.get_jdPart(), &x, &y, &s);

        double rc2i[3][3];
        iauC2ixys(x, y, s, rc2i);

        const double era = iauEra00(UT1.get_jdDay(), UT1.get_jdPart());
        double rc2ti[3][3];
        iauCr(rc2i, rc2ti);
        iauRz(era, rc2ti);

        const double sp = iauSp00(TT.get_jdDay(), TT.get_jdPart());
        double rpom[3][3];
        Pole pole = getPole(UTC.from_jd_to_mjd());
        iauPom00(pole.xp, pole.yp, sp, rpom);
        double rc2it[3][3];
        iauRxr(rpom, rc2ti, rc2it);

        Eigen::Matrix<double, 3, 3> M {{rc2it[0][0], rc2it[0][1], rc2it[0][2]},
                                  {rc2it[1][0], rc2it[1][1], rc2it[1][2]},
                                  {rc2it[2][0], rc2it[2][1], rc2it[2][2]}};

        return Eigen::Quaterniond(M);
    }

    template<>
    inline Eigen::Quaterniond FrameConverter::calcQuaternion<Frame::ITRS, Frame::GCRS>(const Time::Time<Time::Scale::TT> &TT,
                                                                                const Time::Time<Scale::UT1> &UT1,
                                                                                const Time::Time<Time::Scale::UTC> &UTC)
    const noexcept
    {
        return calcQuaternion<Frame::GCRS, Frame::ITRS>(TT, UT1, UTC).conjugate();
    }

}
#endif //BALLISTICS_CONVERTER_H
