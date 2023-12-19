//
// Created by Арсений Плахотнюк on 09.03.2023.
//

#ifndef BALLISTICS_RIGHTPARTOFDIFFEQ_HPP
#define BALLISTICS_RIGHTPARTOFDIFFEQ_HPP

#include "ballistics/time/Time.hpp"
#include "ballistics/forces/DragForce.hpp"
#include "ballistics/forces/SolarPressureForse.hpp"
#include "ballistics/forces/ThirdBodyGravity.hpp"
#include "ballistics/types/BasicTypes.hpp"
#include "ballistics/time/TimeConverter.hpp"
#include "ballistics/forces/GravityForce.hpp"
#include "ballistics/ephemeris/Ephemeris.hpp"
#include "ballistics/shadow/Shadow.hpp"
#include "ballistics/atmosphere/Atmosphere.hpp"


namespace Ballistics::Integrators {

    struct RightPartPointPotential {
        // Точечный потенциал

        // dr
        // -- = v
        // dt

        // dv      /  mu \
        // -- = - (  ---  ) * r
        // dt      \ r^3 /

        constexpr static index size = 6;

        struct Params{
            const scalar gravParam_;
        };

        using TimeState = scalar;

        struct State {
            Vector<scalar, size> u;
            TimeState t;
        };

        [[nodiscard]] static inline Vector<scalar, size> calc(const State &stateVector,
                                                              const Params &params) {
            const scalar rSquaredNorm =
                    Vector<scalar, size / 2>{stateVector.u[0], stateVector.u[1], stateVector.u[2]}.squaredNorm();
            const scalar coeff = (-1) * (params.gravParam_ / (std::sqrt(rSquaredNorm) * rSquaredNorm));
            return Vector<scalar, size>{stateVector.u[3],
                                        stateVector.u[4],
                                        stateVector.u[5],
                                        coeff * stateVector.u[0],
                                        coeff * stateVector.u[1],
                                        coeff * stateVector.u[2]};
        };
    };


    struct RightPartRealPotential{

        // dr/dt = v

        // dv/dt = F(t)/m

        constexpr static index size = 6;

        struct Params{
            const scalar mass_;
            const Forces::GravityForce& gravityForce_;
            const Frame::FrameConverter& frameConverter_;
            const Time::TimeConverter& timeConverter_;
        };

        using TimeState = Ballistics::Time::Time<Ballistics::Time::Scale::TT>;

        struct State {
            Vector<scalar, size> u;
            TimeState t;
        };

        [[nodiscard]] static inline Vector<scalar, size> calc(const State &stateVector,
                                                              const Params &params) {

            Vector3d u1 = stateVector.u.segment<3>(0);
            using TimeTt = Ballistics::Time::Time<Ballistics::Time::Scale::TT>;
            Vector3d force = 1 / params.mass_ * params.gravityForce_.calcForce(params.mass_,
                                                                               stateVector.u.segment<3>(0),
                                                                               stateVector.t,
                                                                               params.frameConverter_,
                                                                               params.timeConverter_);

            return Vector<scalar, size>{stateVector.u[3],
                                        stateVector.u[4],
                                        stateVector.u[5],
                                        force[0],
                                        force[1],
                                        force[2]};
        };
    };


    struct RightPartAllForces
    {
        // dr/dt = v

        // dv/dt = F(t)/m

        constexpr static index size = 6;

        struct Params {
            const scalar mass_;

            // Гравитационный потенциал Земли (GravityForce)
            const Forces::GravityForce& gravityForce_;
            const Frame::FrameConverter& frameConverter_;
            const Time::TimeConverter& timeConverter_;

            // Солнечное давление (SolarPressureForce)
            const scalar solarRadius_;
            const scalar earthRadius_;
            const scalar S_;
            const scalar TSI_;
            const scalar AU_;



            // Гравитация от других небесных тел (ThirdBodyGravity)
            const std::vector<integer> bodiesNumbers_;
            const std::vector<scalar> bodiesGravParams_;
            const Ephemeris::Ephemeris ephems_;

            // Сила трения об атмосферу
            const Atmosphere::Atmosphere atmosphere_;
            const scalar Cd_ = 2.2;
        };

        using TimeState = Ballistics::Time::Time<Ballistics::Time::Scale::TT>;

        struct State {
            Vector<scalar, size> u;
            TimeState t;
        };

        [[nodiscard]] static inline Vector<scalar, size> calc(const State &stateVector,
                                                              const Params &params) {

            const Vector3d satPosGCRS = stateVector.u.segment<3>(0);
            const Vector3d earthGravForce = 1 / params.mass_ * params.gravityForce_.calcForce(params.mass_,
                                                                                              satPosGCRS,
                                                                                              stateVector.t,
                                                                                              params.frameConverter_,
                                                                                              params.timeConverter_);

            const Vector3d drugForce = Forces::DragForce::calcForce(
                    stateVector.u.segment<3>(3),
                    params.atmosphere_.getDensity(satPosGCRS, params.timeConverter_.convert<Time::Scale::UTC>(stateVector.t).value()),
                    params.S_,
                    params.Cd_
                    );

            const Vector3d geoSunPos = 1000 * params.ephems_.calcCoordinates(
                    params.timeConverter_.convert<Time::Scale::TDB>(stateVector.t).value(),
                    Ephemeris::Body::Sun,
                    Ephemeris::Body::Earth);

            const Vector3d helioSatPos = stateVector.u.segment<3>(0) - geoSunPos;

            const scalar F = Shadow::calcShadowFunction(satPosGCRS,
                                                        geoSunPos,
                                                        params.solarRadius_,
                                                        params.earthRadius_);

            const Vector3d solarPressureForce = Forces::SolarPressureForce::calcForce(helioSatPos,
                                                                                      F,
                                                                                      params.S_,
                                                                                      params.TSI_,
                                                                                      params.AU_);
            // Calculating of the gravitational force from space bodies except Earth
            std::vector<Vector3d> BodyPositions;
            for (integer i = 0; i < params.bodiesNumbers_.size(); ++i) {
                BodyPositions.emplace_back(params.ephems_.calcCoordinates(
                        params.timeConverter_.convert<Time::Scale::TDB>(stateVector.t).value(),
                        params.bodiesNumbers_[i],
                        Ephemeris::Body::Earth) * 1000);
            }

            Vector3d thirdBodyGravityForce = Forces::ThirdBodyGravityForce::calcForce(params.mass_,
                                                                                      satPosGCRS,
                                                                                      BodyPositions,
                                                                                      params.bodiesGravParams_);

            Vector3d ResForce = earthGravForce + drugForce + thirdBodyGravityForce + solarPressureForce;

            return Vector<scalar, size>{stateVector.u[3],
                                        stateVector.u[4],
                                        stateVector.u[5],
                                        ResForce[0],
                                        ResForce[1],
                                        ResForce[2]};
        };


    };


    // du   u
    // -- = --
    // dt   t
    class RightPart1 {
    public:
        constexpr static index size = 1;

        struct Params {
            scalar gravParam;
        };

        using Time = scalar;

        struct State {
            vector<scalar, size> u;
            Time t;
        };

        [[nodiscard]] static inline vector<scalar, size> calc(const State &stateVector,
                                                              const Params &params) {
            return vector<scalar, size>{stateVector.u / stateVector.t};
        };
    };

    // du   1
    // -- = --
    // dt   u
    class RightPart2 {
    public:
        constexpr static index size = 1;

        struct Params {
            scalar gravParam;
        };

        using Time = scalar;

        struct State {
            vector<scalar, size> u;
            Time t;
        };

        [[nodiscard]] static inline vector<scalar, size> calc(const State &stateVector,
                                                              const Params &params) {
            return vector<scalar, size>{1. / stateVector.u[0]};
        };
    };

    // du
    // -- = sin(t) * u^2
    // dt
    class RightPart3 {
    public:
        constexpr static index size = 1;

        struct Params {
            scalar gravParam;
        };

        using Time = scalar;

        struct State {
            vector<scalar, size> u;
            Time t;
        };

        [[nodiscard]] static inline vector<scalar, size> calc(const State &stateVector,
                                                              const Params &params) {
            return vector<scalar, size>{sin(stateVector.t) * stateVector.u.dot(stateVector.u)};
        };
    };

    // du
    // -- = cos(t)
    // dt
    class RightPart4 {
    public:
        constexpr static index size = 1;

        struct Params {
            scalar gravParam;
        };

        using Time = scalar;

        struct State {
            vector<scalar, size> u;
            Time t;
        };

        [[nodiscard]] static inline vector<scalar, size> calc(const State &stateVector,
                                                              const Params &params) {
            return vector<scalar, size>{cos(stateVector.t)};
        };
    };

    // dr
    // -- = v
    // dt
    // dv      /  mu \
    // -- = - (  ---  ) * r
    // dt      \ r^3 /
    class RightPart5 {
    public:
        constexpr static index size = 6;

        struct Params {
            scalar gravParam;
        };

        using Time = scalar;

        struct State {
            vector<scalar, size> u;
            Time t;
        };

        [[nodiscard]] static inline vector<scalar, size> calc(const State &stateVector,
                                                              const Params &params) {
            const scalar rSquaredNorm =
                    vector<scalar, size / 2>{stateVector.u[0], stateVector.u[1], stateVector.u[2]}.squaredNorm();
            const scalar coeff = (-1) * (params.gravParam / (sqrt(rSquaredNorm) * rSquaredNorm));
            return vector<scalar, size>{stateVector.u[3],
                                        stateVector.u[4],
                                        stateVector.u[5],
                                        coeff * stateVector.u[0],
                                        coeff * stateVector.u[1],
                                        coeff * stateVector.u[2]};
        };
    };

    // Arenstorff orbit
    class RightPart6 {
    public:
        constexpr static index size = 4;

        struct Params {
            static constexpr scalar mu = 0.012277471;
            static constexpr scalar eta = 1 - mu;
        };

        using Time = scalar;

        struct State {
            vector<scalar, size> u;
            Time t;
        };

        [[nodiscard]] static inline vector<scalar, size> calc(const State &stateVector,
                                                              const Params &params) {
            const scalar x = stateVector.u[0];
            const scalar u = stateVector.u[1];
            const scalar y = stateVector.u[2];
            const scalar nu = stateVector.u[3];
            const scalar multiplier1 = (x + params.mu) * (x + params.mu) + y * y;
            const scalar A = sqrt(multiplier1 * multiplier1 * multiplier1);
            const scalar multiplier2 = (x - params.eta) * (x - params.eta) + y * y;
            const scalar B = sqrt(multiplier2 * multiplier2 * multiplier2);
            return vector<scalar, size>{u,
                                        x + 2 * nu - params.eta * ((x + params.mu) / A) - params.mu * ((x - params.eta) / B),
                                        nu,
                                        y - 2 * u - params.eta * (y / A) - params.mu * (y / B)};
        };
    };



}

#endif //BALLISTICS_RIGHTPARTOFDIFFEQ_HPP
