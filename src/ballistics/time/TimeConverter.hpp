//
// Created by Арсений Плахотнюк on 09.10.2022.
//

#ifndef BALLISTICS_CONVERTER_HPP
#define BALLISTICS_CONVERTER_HPP
#include "Time.hpp"
#include "ballistics/types/BasicTypes.hpp"
#include "sofa/sofa.h"
#include "cmath"
#include <utility>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include "expected/Expected.h"
#include "g3log/g3log.hpp"
#include "Error.hpp"
#include "fast-cpp-csv-parser/csv.h"

namespace Ballistics::Time {
    class TimeConverter {
        std::vector<scalar> dut_;
        integer dut_start_mjd_;
    public:
        [[nodiscard]] unsigned long get_dut_size()
        {
            return dut_.size();
        }
        void bulletin_parser(std::string const &filename) {
#ifdef INITIALIZE_LOG
            LOG(INFO) << "Вызов функции парсинга CSV файла";
#endif
            io::CSVReader<2> data(filename);
            integer a;
            scalar b;
            std::vector<integer> mjd;
            std::vector<scalar> dUT;
            data.read_header(io::ignore_extra_column, "mjd", "UT1-UTC s");
            while (data.read_row(a, b)) {
                mjd.emplace_back(a);
                dUT.emplace_back(b);
            }
            dut_ = dUT;
            dut_start_mjd_ = mjd[0];
        }

        inline TimeConverter(std::string const &filename) { bulletin_parser(filename); };

        static scalar extrapolate_dut(integer day, scalar part) noexcept {
            // received all parameters in MJD
            scalar T = iauEpb(day + 2400000, part + 0.5);
            scalar UT2_UT1 = 0.022 * sin(2 * M_PI * T) - 0.012 * cos(2 * M_PI * T)
                             - 0.006 * sin(4 * M_PI * T) + 0.007 * cos(4 * M_PI * T);
            return -0.0417 + 0.00035 * (day + part - 59852) - (UT2_UT1);
        }

        [[nodiscard]] Expected<scalar, Error> get_dut(integer day, scalar part) const noexcept {
#ifdef INITIALIZE_LOG
            LOG(INFO) << "Вызов функции для расчета DUT";
#endif
            // recieved all parametrs in JD
            // convertation to MJD
            scalar mjd_frac = part + 0.5;
            integer mjd_days = day - 2400001;
            while (mjd_frac > 1){
                mjd_frac -= 1;
                mjd_days += 1;
            }
            day = mjd_days;
            part = mjd_frac;

            if (dut_start_mjd_ <= day && day < dut_start_mjd_ + static_cast<integer>(dut_.size()) - 1) {
                return dut_[day - dut_start_mjd_] +
                       (dut_[day - dut_start_mjd_ + 1] - dut_[day - dut_start_mjd_]) * part;
            } else if (day == dut_start_mjd_ + static_cast<integer>(dut_.size()) - 1) {
                return dut_[day - dut_start_mjd_];
            } else if (day > dut_start_mjd_ + static_cast<integer>(dut_.size()) - 1) {
                return extrapolate_dut(day, part);
            }

#ifdef INITIALIZE_LOG
            LOG(WARNING) << "Некорректные данные для нахожения DUT! Код ошибки:" << Error::DUT_ERROR;
#endif
            return make_unexpected(Error::DUT_ERROR);
        }

        template<Scale To, Scale From>
        Expected<Time <To>, Error>
        convert(const Time <From> &time) const noexcept;

    };

    template<>
    inline Expected<Time < Scale::UTC>, Error>
    TimeConverter::convert(const Time <Scale::TAI> &time) const noexcept {

#ifdef INITIALIZE_LOG
        LOG(INFO) << "Вызов функции для конвертации TAI->UTC";
#endif
        double day, part;
        int status = iauTaiutc(time.get_jdDay(), time.get_jdPart(), &day, &part);
        if (status == -1) {
            LOG(WARNING) << "Вызвана функция для неприемлемой даты! Код ошибки:" << Error::UNACCEPTABLE_DATE;
            return make_unexpected(Error::UNACCEPTABLE_DATE);
        } else if (status == 1) {
            LOG(WARNING) << "Вызвана функция для даты, находящейся слишком далеко в будущем! "
                            "Результатам работы функции не стоит доверять! Код ошибки:" << Error::DUBIOUS_YEAR;
            return make_unexpected(Error::DUBIOUS_YEAR);
        }
        return Time<Scale::UTC>(static_cast<integer>(day), static_cast<scalar>(part));
    }

    template<>
    inline Expected<Time < Scale::TAI>, Error>
    TimeConverter::convert(const Time <Scale::TT> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации TT->TAI";
        double day, part;
        iauTttai(time.get_jdDay(), time.get_jdPart(), &day, &part);
        return Time<Scale::TAI>(static_cast<integer>(day), static_cast<scalar>(part));
    }


    template<>
    inline Expected<Time < Scale::UTC>, Error>
    TimeConverter::convert(const Time <Scale::UT1> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации UT1->UTC";
        double day = time.get_jdDay();
        double part = time.get_jdPart();
        int status = 0;
        for (int i = 0; i < 3; i++) {
            scalar dut = get_dut(static_cast<int>(day), part).value();
            status = iauUt1utc(time.get_jdDay(), time.get_jdPart(), dut, &day, &part);
        }
        if (status == -1) {
            LOG(WARNING) << "Вызвана функция для неприемлемой даты! Код ошибки:" << Error::UNACCEPTABLE_DATE;
            return make_unexpected(Error::UNACCEPTABLE_DATE);
        } else if (status == 1) {
            LOG(WARNING) << "Вызвана функция для даты, находящейся слишком далеко в будущем! "
                            "Результатам работы функции не стоит доверять! Код ошибки:" << Error::DUBIOUS_YEAR;
            return make_unexpected(Error::DUBIOUS_YEAR);
        }
        return Time<Scale::UTC>(static_cast<integer>(day), static_cast<scalar>(part));
    }

    template<>
    inline Expected<Time < Scale::UT1>, Error>
    TimeConverter::convert(const Time <Scale::UTC> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации UTC->UT1";
        double day, part;
        scalar dut = get_dut(time.get_jdDay(), time.get_jdPart()).value();
        int status = iauUtcut1(time.get_jdDay(), time.get_jdPart(), dut, &day, &part);
        if (status == -1) {
            LOG(WARNING) << "Вызвана функция для неприемлемой даты! Код ошибки:" << Error::UNACCEPTABLE_DATE;
            return make_unexpected(Error::UNACCEPTABLE_DATE);
        } else if (status == 1) {
            LOG(WARNING) << "Вызвана функция для даты, находящейся слишком далеко в будущем! "
                            "Результатам работы функции не стоит доверять! Код ошибки:" << Error::DUBIOUS_YEAR;
            return make_unexpected(Error::DUBIOUS_YEAR);
        }
        return Time<Scale::UT1>(static_cast<integer>(day), static_cast<scalar>(part));
    }

    template<>
    inline Expected<Time < Scale::UTC>, Error>
    TimeConverter::convert(const Time <Scale::TT> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации TT->UTC";
        return convert < Scale::UTC > (convert < Scale::TAI > (time).value());
    }

    template<>
    inline Expected<Time < Scale::UT1>, Error>
    TimeConverter::convert(const Time <Scale::TT> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации TT->UT1";
        return convert<Scale::UT1> (convert<Scale::UTC>(time).value());
    }

    template<>
    inline Expected<Time < Scale::TCG>, Error>
    TimeConverter::convert(const Time <Scale::TT> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации TT->TCG";
        double day, part;
        iauTttcg(time.get_jdDay(), time.get_jdPart(), &day, &part);
        return Time<Scale::TCG>(static_cast<integer>(day), static_cast<scalar>(part));
    }

    template<>
    inline Expected<Time < Scale::TDB>, Error>
    TimeConverter::convert(const Time <Scale::TT> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации TT->TDB";
        double day, part;
        Time<Scale::UT1> ut1 = convert < Scale::UT1 > (time).value();
        iauTttdb(time.get_jdDay(), time.get_jdPart(), iauDtdb(time.get_jdDay(), time.get_jdPart(),
                                                              ut1.get_jdPart(), 0, 0, 0), &day, &part);
        return Time<Scale::TDB>(static_cast<integer>(day), static_cast<scalar>(part));
    }

    template<>
    inline Expected<Time < Scale::TDB>, Error>
    TimeConverter::convert(const Time <Scale::TCB> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации TCB->TDB";
        double day, part;
        iauTcbtdb(time.get_jdDay(), time.get_jdPart(), &day, &part);
        return Time<Scale::TDB>(static_cast<integer>(day), static_cast<scalar>(part));
    }

    template<>
    inline Expected<Time < Scale::TCB>, Error>
    TimeConverter::convert(const Time <Scale::TDB> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации TDB->TCB";
        double day, part;
        iauTdbtcb(time.get_jdDay(), time.get_jdPart(), &day, &part);
        return Time<Scale::TCB>(static_cast<integer>(day), static_cast<scalar>(part));
    }

    template<>
    inline Expected<Time < Scale::TAI>, Error>
    TimeConverter::convert(const Time <Scale::UTC> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации UTC->TAI";
        double day, part;
        int status = iauUtctai(time.get_jdDay(), time.get_jdPart(), &day, &part);
        if (status == -1) {
            LOG(WARNING) << "Вызвана функция для неприемлемой даты! Код ошибки:" << Error::UNACCEPTABLE_DATE;
            return make_unexpected(Error::UNACCEPTABLE_DATE);
        } else if (status == 1) {
            LOG(WARNING) << "Вызвана функция для даты, находящейся слишком далеко в будущем! "
                            "Результатам работы функции не стоит доверять! Код ошибки:" << Error::DUBIOUS_YEAR;
            return make_unexpected(Error::DUBIOUS_YEAR);
        }
        return Time<Scale::TAI>(static_cast<integer>(day), static_cast<scalar>(part));
    }

    template<>
    inline Expected<Time < Scale::TT>, Error>
    TimeConverter::convert(const Time <Scale::TDB> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации TDB->TT";
        scalar day = 0., part = 0.;
        const scalar dtdb = iauDtdb(time.get_jdDay(), time.get_jdPart(), 0, 0, 0, 0);
        iauTdbtt(time.get_jdDay(), time.get_jdPart(), dtdb, &day, &part);

        return Time<Scale::TT>(static_cast<integer>(day), static_cast<scalar>(part));
    }

    template<>
    inline Expected<Time < Scale::TT>, Error>
    TimeConverter::convert(const Time <Scale::TAI> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации TAI->TT";
        double day, part;
        iauTaitt(time.get_jdDay(), time.get_jdPart(), &day, &part);
        return Time<Scale::TT>(static_cast<integer>(day), static_cast<scalar>(part));
    }

    template<>
    inline Expected<Time < Scale::TT>, Error>
    TimeConverter::convert(const Time <Scale::UTC> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации UTC->TT";
        return convert<Scale::TT>(convert<Scale::TAI>(time).value());
    }

    template<>
    inline Expected<Time < Scale::TT>, Error>
    TimeConverter::convert(const Time <Scale::UT1> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации UT1->TT";
        return convert<Scale::TT>(convert<Scale::UTC>(time).value());
    }

    template<>
    inline Expected<Time < Scale::TT>, Error>
    TimeConverter::convert(const Time <Scale::TT> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации TT->TT";
        return time;
    }

    template<>
    inline Expected<Time < Scale::UTC>, Error>
    TimeConverter::convert(const Time <Scale::UTC> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации UTC->UTC";
        return time;
    }

    template<>
    inline Expected<Time < Scale::TCG>, Error>
    TimeConverter::convert(const Time <Scale::UTC> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации UTC->TCG";
        return convert < Scale::TCG > (convert < Scale::TT > (time).value());
    }

    template<>
    inline Expected<Time < Scale::UT1>, Error>
    TimeConverter::convert(const Time <Scale::TAI> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации TAI->UT1";
        return convert < Scale::UT1 > (convert < Scale::UTC > (time).value());
    }

    template<>
    inline Expected<Time < Scale::TCG>, Error>
    TimeConverter::convert(const Time <Scale::TAI> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации TAI->TCG";
        return convert < Scale::TCG > (convert < Scale::TT > (time).value());
    }

    template<>
    inline Expected<Time < Scale::TAI>, Error>
    TimeConverter::convert(const Time <Scale::UT1> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации UT1->TAI";
        return convert < Scale::TAI > (convert < Scale::UTC > (time).value());
    }


    template<>
    inline Expected<Time < Scale::TAI>, Error>
    TimeConverter::convert(const Time <Scale::TAI> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации TAI->TAI";
        return time;
    }

    template<>
    inline Expected<Time < Scale::TDB>, Error>
    TimeConverter::convert(const Time <Scale::TAI> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации TAI->TDB";
        return convert < Scale::TDB > (convert < Scale::TT > (time).value());
    }

    template<>
    inline Expected<Time < Scale::UT1>, Error>
    TimeConverter::convert(const Time <Scale::UT1> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации UT1->UT1";
        return time;
    }

    template<>
    inline Expected<Time < Scale::TDB>, Error>
    TimeConverter::convert(const Time <Scale::UT1> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации UT1->TDB";
        return convert < Scale::TDB > (convert < Scale::TT > (time).value());
    }

    template<>
    inline Expected<Time < Scale::TCG>, Error>
    TimeConverter::convert(const Time <Scale::UT1> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации UT1->TCG";
        return convert < Scale::TCG > (convert < Scale::TT > (time).value());
    }

    template<>
    inline Expected<Time < Scale::TT>, Error>
    TimeConverter::convert(const Time <Scale::TCG> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации TCG->TT";
        double day, part;
        iauTcgtt(time.get_jdDay(), time.get_jdPart(), &day, &part);
        return Time<Scale::TT>(static_cast<integer>(day), static_cast<scalar>(part));
    }

    template<>
    inline Expected<Time < Scale::TCG>, Error>
    TimeConverter::convert(const Time <Scale::TCG> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации TCG->TCG";
        return time;
    }

    template<>
    inline Expected<Time < Scale::TAI>, Error>
    TimeConverter::convert(const Time <Scale::TCG> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации TCG->TAI";
        return convert < Scale::TAI > (convert < Scale::TT > (time).value());
    }

    template<>
    inline Expected<Time < Scale::UTC>, Error>
    TimeConverter::convert(const Time <Scale::TCG> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации TCG->UTC";
        return convert < Scale::UTC > (convert < Scale::TT > (time).value());
    }

    template<>
    inline Expected<Time < Scale::UT1>, Error>
    TimeConverter::convert(const Time <Scale::TCG> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации TCG->UT1";
        return convert < Scale::UT1 > (convert < Scale::TT > (time).value());
    }

    template<>
    inline Expected<Time < Scale::TDB>, Error>
    TimeConverter::convert(const Time <Scale::TCG> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации TCG->TDB";
        return convert < Scale::TDB > (convert < Scale::TT > (time).value());
    }

    template<>
    inline Expected<Time < Scale::TCB>, Error>
    TimeConverter::convert(const Time <Scale::TCG> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации TCG->TCB";
        return convert < Scale::TCB > (convert < Scale::TDB > (time).value());
    }

    template<>
    inline Expected<Time < Scale::TCB>, Error>
    TimeConverter::convert(const Time <Scale::TT> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации TT->TCB";
        return convert < Scale::TCB > (convert < Scale::TDB > (time).value());
    }

    template<>
    inline Expected<Time < Scale::TCB>, Error>
    TimeConverter::convert(const Time <Scale::TAI> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации TAI->TCB";
        return convert < Scale::TCB > (convert < Scale::TT > (time).value());
    }

    template<>
    inline Expected<Time < Scale::TCB>, Error>
    TimeConverter::convert(const Time <Scale::UTC> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации UTC->TCB";
        return convert < Scale::TCB > (convert < Scale::TT > (time).value());
    }

    template<>
    inline Expected<Time < Scale::TCB>, Error>
    TimeConverter::convert(const Time <Scale::UT1> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации UT1->TCB";
        return convert < Scale::TCB > (convert < Scale::TT > (time).value());
    }

    template<>
    inline Expected<Time < Scale::TCB>, Error>
    TimeConverter::convert(const Time <Scale::TCB> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации TCB->TCB";
        return time;
    }

    template<>
    inline Expected<Time < Scale::TDB>, Error>
    TimeConverter::convert(const Time <Scale::TDB> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации TDB->TDB";
        return time;
    }

    template<>
    inline Expected<Time < Scale::UTC>, Error>
    TimeConverter::convert(const Time <Scale::TDB> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации TDB->UTC";
        return convert < Scale::UTC > (convert < Scale::TT > (time).value());
    }

    template<>
    inline Expected<Time < Scale::TAI>, Error>
    TimeConverter::convert(const Time <Scale::TDB> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации TDB->TAI";
        return convert < Scale::TAI > (convert < Scale::TT > (time).value());
    }

    template<>
    inline Expected<Time < Scale::UT1>, Error>
    TimeConverter::convert(const Time <Scale::TDB> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации TDB->UT1";
        return convert < Scale::UT1 > (convert < Scale::TT > (time).value());
    }

    template<>
    inline Expected<Time < Scale::TT>, Error>
    TimeConverter::convert(const Time <Scale::TCB> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации TCB->TT";
        return convert < Scale::TT > (convert < Scale::TDB > (time).value());
    }

    template<>
    inline Expected<Time < Scale::TCG>, Error>
    TimeConverter::convert(const Time <Scale::TCB> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации TCB->TCG";
        return convert < Scale::TCG > (convert < Scale::TT > (time).value());
    }

    template<>
    inline Expected<Time < Scale::TCG>, Error>
    TimeConverter::convert(const Time <Scale::TDB> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации TDB->TCG";
        return convert < Scale::TCG > (convert < Scale::TT > (time).value());
    }

    template<>
    inline Expected<Time < Scale::UT1>, Error>
    TimeConverter::convert(const Time <Scale::TCB> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации TCB->UT1";
        return convert < Scale::UT1 > (convert < Scale::TT > (time).value());
    }

    template<>
    inline Expected<Time < Scale::UTC>, Error>
    TimeConverter::convert(const Time <Scale::TCB> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации TCB->UTC";
        return convert < Scale::UTC > (convert < Scale::TT > (time).value());
    }

    template<>
    inline Expected<Time < Scale::TAI>, Error>
    TimeConverter::convert(const Time <Scale::TCB> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации TCB->TAI";
        return convert < Scale::TAI > (convert < Scale::TT > (time).value());
    }

    template<>
    inline Expected<Time < Scale::TDB>, Error>
    TimeConverter::convert(const Time <Scale::UTC> &time) const noexcept {
        LOG(INFO) << "Вызов функции для конвертации UTC->TDB";
        return convert < Scale::TDB > (convert < Scale::TT > (time).value());
    }
}


#endif //BALLISTICS_CONVERTER_HPP
