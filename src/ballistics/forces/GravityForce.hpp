//
// Created by Арсений Плахотнюк on 10.11.2022.
//

#ifndef BALLISTICS_GRAVITYFORCE_HPP
#define BALLISTICS_GRAVITYFORCE_HPP
#include <vector>
#include <geographiclib/include/GeographicLib/GravityModel.hpp>
#include <ballistics/types/BasicTypes.hpp>
#include <ballistics/time/Time.hpp>
#include "ballistics/frame_converter/FrameConverter.hpp"
#include "ballistics/time/TimeConverter.hpp"

namespace Ballistics::Forces {

    [[nodiscard]] Vector3d CalcGravityForce(const Ballistics::scalar mass,
                                                   const Ballistics::Vector3d &rGCRS,
                                                   const Ballistics::quaterniond &q,
                                                   const Time::Time<Time::Scale::TT> &timeTT,
                                                   const Time::Time<Time::Scale::UT1> &timeUT1) const noexcept
    {
        const Vector3d rITRS = q._transformVector(rGCRS);
        Vector3d aITRS;
        model_.V(rITRS.x(), rITRS.y(), rITRS.z(), aITRS.x(), aITRS.y(), aITRS.z());
        const Vector3d aGCRS = q.conjugate()._transformVector(aITRS);
        return aGCRS * mass;
    }



    struct GravityForceParams{
        scalar mass_;
        Vector3d rGcrs_;
        Quaternion gcrs2Itrs_;
    };

    class GravityForce {
    private:
        const GeographicLib::GravityModel GM_;

    public:
        GravityForce(const std::string& modelName, const std::string& fileLocation,
                     const integer n = -1, const integer m = -1) : GM_(modelName, fileLocation, n, m) {
            /**
             * Конструктор класса гравитационной силы
             * @param modelName - название модели гравитационного потенциала
             * @param fileLocation - путь к файлу с данными модели
             * @param n - степень разложения гравитационного потенциала
             * @param m - порядок разложения потенциала
             * @return вектор силы в системе GCRS
             */
        };


        [[nodiscard]] Vector3d calcForce(const scalar mass,
                                         const Vector3d& rGcrs,
                                         const Quaternion& gcrs2Itrs) const
        {
            /**
             * Функция пересчитывает гравитационную силу с учетом вращения Земли
             * @param mass - масса объекта
             * @param rGcrs - вектор положения КА в GCRS
             * @param gcrs2Itrs - кватернион перехода из GCRS в ITRS
             * @return вектор силы в системе GCRS
             */

            const Vector3d r_itrs = gcrs2Itrs._transformVector(rGcrs);
            Vector3d g;
            GM_.V(r_itrs.x(), r_itrs.y(), r_itrs.z(), g.x(), g.y(), g.z());
            const Vector3d acceleration = gcrs2Itrs.conjugate()._transformVector(g);
            return g * mass;
        }

        [[nodiscard]] Vector3d calcForce(const GravityForceParams& params) const
        {
            /**
             * Функция пересчитывает гравитационную силу с учетом вращения Земли
             * @param params - структура, содержащая параметры для расчета силы
             * @return вектор силы в шкале GCRS
             */

            const Vector3d rItrs = params.gcrs2Itrs_._transformVector(params.rGcrs_);
            Vector3d g;
            GM_.V(rItrs.x(), rItrs.y(), rItrs.z(), g.x(), g.y(), g.z());
            const Vector3d acceleration = params.gcrs2Itrs_.conjugate()._transformVector(g);
            return g * params.mass_;
        };

        [[nodiscard]] Vector3d calcForce(const scalar mass,
                                         const Vector3d& rGcrs,
                                         const Time::Time<Time::Scale::TT>& timeTt,
                                         const Frame::FrameConverter& frameConverter,
                                         const Time::TimeConverter& timeConverter) const
        {
            /**
             * Функция пересчитывает гравитационную силу с учетом вращения Земли
             * @param GM - модель гравитационного поленциала Земли
             * @param mass - масса объекта
             * @param rGcrs - вектор положения объекта в GCRS
             * @param timeTt - время в шкале TT
             * @param frameConverter - класс конвертации систем координат
             * @param timeConverter - класс конвертации шкал времен
             * @return вектор силы в шкале GCRS
             */

            Time::Time<Time::Scale::UT1> timeUt1 = timeConverter.convert<Time::Scale::UT1>(timeTt).value();
            Time::Time<Time::Scale::UTC> timeUtc = timeConverter.convert<Time::Scale::UTC>(timeTt).value();

            Eigen::Quaterniond gcrs2Itrs = frameConverter.calcQuaternion<Frame::GCRS, Frame::ITRS>(timeTt, timeUt1, timeUtc);
            const Vector3d rItrs = gcrs2Itrs._transformVector(rGcrs);

            Vector3d aITRS;
            GM_.V(rItrs.x(), rItrs.y(), rItrs.z(), aITRS.x(), aITRS.y(), aITRS.z());
            const Vector3d aGCRS = gcrs2Itrs.conjugate()._transformVector(aITRS);
            return aGCRS * mass;
        };

        [[nodiscard]] Vector3d calcForce(const scalar mass,
                                         const Vector3d &rGCRS,
                                         const quaterniond &q,
                                         const Time::Time<Time::Scale::TT> &timeTT,
                                         const Time::Time<Time::Scale::UT1> &timeUT1) const noexcept{
            const Vector3d rITRS = q._transformVector(rGCRS);
            Vector3d aITRS;
            GM_.V(rITRS.x(), rITRS.y(), rITRS.z(), aITRS.x(), aITRS.y(), aITRS.z());
            const Vector3d aGCRS = q.conjugate()._transformVector(aITRS);
            return aGCRS * mass;

        };


    };
}
#endif //BALLISTICS_GRAVITYFORCE_HPP
