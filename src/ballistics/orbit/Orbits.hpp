//
// Created by Арсений Плахотнюк on 02.02.2023.
//

#ifndef BALLISTICS_ORBITS_HPP
#define BALLISTICS_ORBITS_HPP

#include "ballistics/types/BasicTypes.hpp"
#include "third_party/fast-cpp-csv-parser/csv.h"

namespace Ballistics::Orbits{

    struct PosVel{
        Vector3d r;
        Vector3d v;
    };

    struct KeplerTrue{
        scalar inclination; // i
        scalar ascendingNodeLongitude; // Omega
        scalar periapsisArgument; // omega
        scalar semiMajorAxis; // a
        scalar eccentricity; // e
        scalar trueAnomaly; // nu
    };

    struct KeplerMean{
        scalar inclination; // i
        scalar ascendingNodeLongitude; // Omega
        scalar periapsisArgument; // omega
        scalar semiMajorAxis; // a
        scalar eccentricity; // e
        scalar meanAnomaly; // M
    };


    template<typename To,typename From>
    [[nodiscard]] To convert(const From& val, const scalar& mu) noexcept;

    template<>
    inline KeplerTrue convert(const PosVel& rv, const scalar& mu) noexcept{
        // Вектор орбитального момента L
        const Vector3d L = rv.r.cross(rv.v);
        // Наклонение орбиты
        const scalar inclination = std::atan2(std::sqrt(L(0)*L(0) + L(1)*L(1)), L(2));

        const Vector3d z = Vector3d ::UnitZ();
        const Vector3d x = Vector3d::UnitX();
        const Vector3d z_cross_L =  z.cross(L);
        const scalar lNorm = z_cross_L.norm();

        // Направление на восходящий узел
        const Vector3d N = lNorm > 0 ? z_cross_L / lNorm : x;

        // Долгота восходящего узла
        const scalar ascending_node_longitude = std::atan2(N(1), N(0));

        // Векторный эксцентриситет
        const Vector3d eccentricity_vec = ((rv.v.dot(rv.v) - mu / rv.r.norm()) * rv.r - rv.r.dot(rv.v)*rv.v)/mu;

        // Эксцентриситет
        const scalar eccentricity = eccentricity_vec.norm();

        // Для исключения вырожение вводится вектор е_1
        const Vector3d e_1 = eccentricity_vec.norm() > 0 ? eccentricity_vec : N;

        // Вводим базис в плоскости орбиты d_1, d_2
        const Vector3d d_1 = N;
        const Vector3d d_2 = L.cross(N) / L.norm();

        // Аргумент перицентра
        const scalar periapsis_argument = std::atan2(d_2.dot(e_1), d_1.dot(e_1));

        // mu/r - v*v/2 = mu/2a; a = mu/(2mu/r - v*v)
        // Главная полуось орбиты a из интеграла энергии
        const scalar semi_major_axis = mu/(2*mu/rv.r.norm() - rv.v.dot(rv.v));

        const Vector3d e_2 = L.cross(e_1) / L.norm();
        // Истинная аномалия
        const scalar true_anomaly = std::atan2(e_2.dot(rv.r), e_1.dot(rv.r));

        return {inclination, ascending_node_longitude, periapsis_argument, semi_major_axis, eccentricity, true_anomaly};

    };

    template<>
    inline KeplerMean convert(const KeplerTrue& kt, const scalar& mu) noexcept{
        // Эксцентрическая аномалия
        const scalar E = std::atan2(std::sin(kt.trueAnomaly) * std::sqrt(1 - kt.eccentricity * kt.eccentricity), kt.eccentricity + std::cos(kt.trueAnomaly));
        // Средняя аномалия
        const scalar mean_anomaly = E - kt.eccentricity*std::sin(E);

        return {kt.inclination, kt.ascendingNodeLongitude, kt.periapsisArgument, kt.semiMajorAxis, kt.eccentricity, mean_anomaly};
    }

    template<>
    inline KeplerTrue convert(const KeplerMean& km, const scalar& mu) noexcept{
        const scalar tolerance = 1.e-15;
        // Метод Ньютона для поиска истинной аномалии по средней
        scalar E = km.meanAnomaly - km.eccentricity;

        std::function<scalar(scalar)> func = [km](scalar E) { return (E - km.eccentricity*std::sin(E) - km.meanAnomaly) / (1 - km.eccentricity * std::cos(E)); };

        while(std::abs(func(E)) > tolerance){
            E = E - func(E);
        }

        const scalar true_anomaly = std::atan2(std::sin(E)*std::sqrt(1 - km.eccentricity*km.eccentricity), std::cos(E) - km.eccentricity);

        return {km.inclination, km.ascendingNodeLongitude, km.periapsisArgument, km.semiMajorAxis, km.eccentricity, true_anomaly};
    }

    template<>
    inline PosVel convert(const KeplerTrue& kt, const scalar& mu) noexcept{
        // Calculations of trigonometric functions of Kepler's elements
        const scalar cosTrueAnomaly           = cos(kt.trueAnomaly);
        const scalar sinTrueAnomaly           = sin(kt.trueAnomaly);
        const scalar cosLongitOfAscendingNode = cos(kt.ascendingNodeLongitude);
        const scalar sinLongitOfAscendingNode = sin(kt.ascendingNodeLongitude);
        const scalar cosPericenterArg         = cos(kt.periapsisArgument);
        const scalar sinPericenterArg         = sin(kt.periapsisArgument);
        const scalar cosInclination           = cos(kt.inclination);
        const scalar sinInclination           = sin(kt.inclination);

        const scalar orbParam = kt.semiMajorAxis * (static_cast<scalar>(1) - kt.eccentricity *
                                                                 kt.eccentricity);
        // Some auxiliary multipliers
        const scalar multiplier1 = orbParam / (static_cast<scalar>(1) + kt.eccentricity * cosTrueAnomaly);
        const scalar multiplier2 = sqrt(mu / orbParam);

        const Vector3d positionTmp = {multiplier1 * cosTrueAnomaly,
                                      multiplier1 * sinTrueAnomaly,
                                      static_cast<scalar>(0)};
        const Vector3d velocityTmp = {-multiplier2 * sinTrueAnomaly,
                                      multiplier2 * (kt.eccentricity + cosTrueAnomaly),
                                      static_cast<scalar>(0)};
        const Matrix3x3 rotationMatrix {{cosLongitOfAscendingNode * cosPericenterArg - sinLongitOfAscendingNode * sinPericenterArg * cosInclination,
                                               -cosLongitOfAscendingNode * sinPericenterArg - sinLongitOfAscendingNode * cosPericenterArg * cosInclination,
                                               sinLongitOfAscendingNode * sinInclination},
                                       {sinLongitOfAscendingNode * cosPericenterArg + cosLongitOfAscendingNode * sinPericenterArg * cosInclination,
                                               -sinLongitOfAscendingNode * sinPericenterArg + cosLongitOfAscendingNode * cosPericenterArg * cosInclination,
                                               -cosLongitOfAscendingNode * sinInclination},
                                       {sinPericenterArg * sinInclination,
                                               cosPericenterArg * sinInclination,
                                               cosInclination}};
        return {rotationMatrix * positionTmp, rotationMatrix * velocityTmp};
    }

    template<>
    inline PosVel convert(const KeplerMean& km, const scalar& mu) noexcept{
        const KeplerTrue kt = convert<KeplerTrue>(km, mu);
        return convert<PosVel>(kt, mu);
    }

    template<>
    inline KeplerMean convert(const PosVel& rv, const scalar& mu) noexcept{
        const KeplerTrue kt = convert<KeplerTrue>(rv, mu);
        return convert<KeplerMean>(kt, mu);
    }


    [[nodiscard]] KeplerMean propagate(const KeplerMean& km, const scalar& mu, const scalar& DeltaT){

        const scalar mean_anomaly = std::sqrt(mu/(km.semiMajorAxis * km.semiMajorAxis * km.semiMajorAxis)) * DeltaT;

        return {km.inclination, km.ascendingNodeLongitude, km.periapsisArgument,
                km.semiMajorAxis, km.eccentricity, km.meanAnomaly + mean_anomaly};
    }

}

#endif //BALLISTICS_ORBITS_HPP
