//
// Created by Арсений Плахотнюк on 07.06.2023.
//

#include "ballistics/integrators/explicitRungeKutta/ButcherTableERK4.hpp"
#include "ballistics/integrators/explicitRungeKutta/ExplicitRungeKutta.hpp"
#include "ballistics/integrators/RightPart/RightPartOfDiffEq.hpp"

#include "ballistics/forces/GravityForce.hpp"
#include "ballistics/forces/ThirdBodyGravity.hpp"
#include "ballistics/forces/SolarPressureForse.hpp"
#include "ballistics/forces/DragForce.hpp"

#include "ballistics/shadow/Shadow.hpp"
#include "ballistics/atmosphere/Atmosphere.hpp"

#include "ballistics/time/TimeConverter.hpp"
#include "ballistics/frame_converter/FrameConverter.hpp"
#include "ballistics/ephemeris/Ephemeris.hpp"
#include "ballistics/orbit/Orbits.hpp"

#include "ballistics/optimization/DifferentialCorrection.hpp"

#include <fstream>
#include <random>
#include <chrono>

const std::string FILE_PATH = __FILE__;
const int amountOfSymbolsBeforeRootDirectory = 39;
const std::string DIR_PATH = FILE_PATH.substr(0, FILE_PATH.size() - amountOfSymbolsBeforeRootDirectory);

using namespace Ballistics;
using namespace Ballistics::Time;

class quaternionCalculator {
    Frame::FrameConverter frameConverter_;
    Ballistics::Time::TimeConverter timeConverter_;
public:
    quaternionCalculator(const std::string &csvTimeLocation,
                         const std::string &csvFrameLocation) noexcept: frameConverter_(csvFrameLocation),
                                                                        timeConverter_(csvTimeLocation) {};

    [[nodiscard]] Eigen::Quaterniond calcQuaternion_(Ballistics::Time::Time<Ballistics::Time::Scale::TT> &timeTT) const {

        const auto timeUT1 = timeConverter_.convert<Ballistics::Time::Scale::UT1>(timeTT).value();
        const auto timeUTC = timeConverter_.convert<Ballistics::Time::Scale::UTC>(timeTT).value();

        return frameConverter_.calcQuaternion<Ballistics::Frame::Frame::GCRS, Ballistics::Frame::Frame::ITRS>(timeTT,
                                                                                                             timeUT1,
                                                                                                             timeUTC);
    }
};

class RightPart34 {
public:

    constexpr static Ballistics::index size = 6;

    struct Params {
        scalar mass;
        scalar solarRadius;
        scalar earthRadius;
        scalar S;
        scalar TSI;
        scalar AU;
        scalar Cd;
        std::vector<integer> bodiesNumbers;
        std::vector<scalar> bodiesGravParams;
        quaternionCalculator quatCalculator;
        Ballistics::Time::Time<Ballistics::Time::Scale::TT> &timeTT;
        Ballistics::Time::Time<Ballistics::Time::Scale::UT1> &timeUT1;
        Ballistics::Time::Time<Ballistics::Time::Scale::UTC> &timeUTC;
        Ballistics::Time::Time<Ballistics::Time::Scale::TDB> &timeTDB;
        Forces::GravityForce model;
        Ephemeris::Ephemeris ephems;
        Atmosphere::Atmosphere atmosphere;
    };

    using Time = Ballistics::Time::Time<Ballistics::Time::Scale::TT>;

    struct State {
        vector<scalar, size> u;
        Time t{};
    };

    [[nodiscard]] static inline vector<scalar, size> calc(const State &stateVector,
                                                          const Params &params) {
        const Vector3d geoSatPos = {stateVector.u[0], stateVector.u[1], stateVector.u[2]};
        // Calculating of the gravitational force from Earth
        const quaterniond GCRStoITRSquat = params.quatCalculator.calcQuaternion_(params.timeTT);
        const Vector3d gravityForceEarth = params.model.calcForce(params.mass,
                                                                  geoSatPos,
                                                                  GCRStoITRSquat,
                                                                  params.timeTT,
                                                                  params.timeUT1);
        // Calculating of the gravitational force from space bodies except Earth
        std::vector<Vector3d> BodyPositions;
        for (integer i = 0; i < params.bodiesNumbers.size(); ++i) {
            BodyPositions.emplace_back(params.ephems.calcCoordinates(params.timeTDB,
                                                                     params.bodiesNumbers[i],
                                                                     Ephemeris::Body::Earth) * 1000);
        }
        const Vector3d gravityForceThirdBodies = Ballistics::Forces::ThirdBodyGravityForce::calcForce(params.mass,
                                                                                                      geoSatPos,
                                                                                                      BodyPositions,
                                                                                                      params.bodiesGravParams);
        // Calculating of the solar pressure force
        const Vector3d geoSunPos = params.ephems.calcCoordinates(params.timeTDB,
                                                                 Ephemeris::Body::Sun,
                                                                 Ephemeris::Body::Earth) * 1000;
        const Vector3d helioSatPos = geoSatPos - geoSunPos;
        const scalar F = Shadow::calcShadowFunction(geoSatPos,
                                                    geoSunPos,
                                                    params.solarRadius,
                                                    params.earthRadius);
        const Vector3d solarPressureForce = Ballistics::Forces::SolarPressureForce::calcForce(helioSatPos,
                                                                           F,
                                                                           params.S,
                                                                           params.TSI,
                                                                           params.AU);
        // Calculating of the drag force
        const Vector3d ITRSSatPos = GCRStoITRSquat._transformVector(geoSatPos);
        const Vector3d dragForce = Ballistics::Forces::DragForce::calcForce(
                                                         {stateVector.u[3], stateVector.u[4], stateVector.u[5]},
                                                         params.atmosphere.getDensity(ITRSSatPos,
                                                                                      params.timeUTC),
                                                         params.S,
                                                         params.Cd);
        // Calculating sum of all forces
        const Vector3d force = gravityForceEarth + gravityForceThirdBodies + solarPressureForce + dragForce;
        return vector<scalar, size>{stateVector.u[3],
                                    stateVector.u[4],
                                    stateVector.u[5],
                                    force[0] / params.mass,
                                    force[1] / params.mass,
                                    force[2] / params.mass};
    };
};

int main() {
    // initial values
    const scalar gravParam = 3.986004415e14;
    const Ballistics::Orbits::KeplerMean meanElements = {0.,
                                            0.,
                                            6800e3,
                                            0.,
                                            0.,
                                            0.};
    const Ballistics::Orbits::PosVel CartesianElements = Ballistics::Orbits::convert<Ballistics::Orbits::PosVel>(meanElements, gravParam);

    // let's integrate this stuff
    const Ballistics::Time::TimeConverter timeConverter = Ballistics::Time::TimeConverter("../../resources/earth_rotation.csv");
    auto timeTT = Ballistics::Time::Time<Ballistics::Time::Scale::TT>(2459842, 0.);
    auto timeUT1 = timeConverter.convert<Ballistics::Time::Scale::UT1>(timeTT).value();
    auto timeUTC = timeConverter.convert<Ballistics::Time::Scale::UTC>(timeTT).value();
    auto timeTDB = timeConverter.convert<Ballistics::Time::Scale::TDB>(timeTT).value();
    const RightPart34::State initState = {vector<scalar, 6>{CartesianElements.r[0],
                                                            CartesianElements.r[1],
                                                            CartesianElements.r[2],
                                                            CartesianElements.v[0],
                                                            CartesianElements.v[1],
                                                            CartesianElements.v[2]}, timeTT};

    const quaternionCalculator quatCalculator = quaternionCalculator("../../resources/earth_rotation.csv",
                                                                     "../../resources/earth_rotation.csv");
    const RightPart34::Params params = {1.,
                                        695.700e+6,
                                        6.371e+6,
                                        0.1,
                                        1366,
                                        149597870700,
                                        2.2,
                                        {Ephemeris::Body::Sun, Ephemeris::Body::MarsBarycenter,
                                         Ephemeris::Body::Moon, Ephemeris::Body::JupiterBarycenter,
                                         Ephemeris::Body::SaturnBarycenter},
                                        {1.3271244e+20, 4.2828e+13,
                                         4.9028e+12, 1.26686534e+17,
                                         3.7931187e+16},
                                        quatCalculator,
                                        timeTT,
                                        timeUT1,
                                        timeUTC,
                                        timeTDB,
                                        {"egm2008", "../../resources/geographicLib", 16, 16},
                                        Ephemeris::Ephemeris(),
                                        Atmosphere::Atmosphere("../../resources/atmosphere/SW-Last5Years.csv")};
    const scalar intStep = 1;

    RightPart34::State tmpState = initState;
    std::vector<RightPart34::State> measurements = {initState};
    for (integer i = 0; i < 50; ++i) {
        const std::vector<typename RightPart34::State> tmpMeasurements =
                Ballistics::Integrators::ExplicitRK::ExplicitRungeKutta<Integrators::ExplicitRK::ButcherTableERK4>::
                calc<RightPart34>(tmpState,
                                  params,
                                  {intStep, tmpState.t + intStep});
        tmpState = tmpMeasurements.back();
        measurements.push_back(tmpState);
    }

    // Here we are adding a noice to the measurements
    unsigned int seed = std::chrono::steady_clock::now().time_since_epoch().count();
    std::default_random_engine gen(seed);

    double meanValueCoord = 0.;
    double standardDeviationCoord = 10.;
    std::normal_distribution<> coordDistribution {meanValueCoord, standardDeviationCoord};

    double meanValueVel = 0.;
    double standardDeviationVel = 0.02;
    std::normal_distribution<> velDistribution {meanValueVel, standardDeviationVel};

    for (auto & measurement : measurements) {
        measurement.u = {measurement.u[0] + coordDistribution(gen),
                         measurement.u[1] + coordDistribution(gen),
                         measurement.u[2] + coordDistribution(gen),
                         measurement.u[3] + velDistribution(gen),
                         measurement.u[4] + velDistribution(gen),
                         measurement.u[5] + velDistribution(gen)};
    }

    // Here we are constructing a matrix W with error of measurements
    const scalar sigmaCoord = 10.;
    const scalar sigmaVel = 0.02;
    const scalar sigmaCoordSqRev = 1. / (sigmaCoord * sigmaCoord);
    const scalar sigmaVelSqRev = 1. / (sigmaVel * sigmaVel);
    const vector<scalar, RightPart34::size> sigmaVector = {sigmaCoordSqRev, sigmaCoordSqRev, sigmaCoordSqRev,
                                                           sigmaVelSqRev, sigmaVelSqRev, sigmaVelSqRev};
    const diagMatrix<scalar, RightPart34::size, RightPart34::size> W(sigmaVector);
    std::vector<diagMatrix<scalar, RightPart34::size, RightPart34::size>> vecW(measurements.size());
    for (auto &el : vecW) {
        el = W;
    }

    const vector<scalar, RightPart34::size> delta = {0.3, 0.3, 0.3, 0., 0., 0.};

    std::cout << "First position:" << std::endl;
    std::cout << initState.u + delta<< std::endl;

    const vector<scalar, RightPart34::size> derivativeSteps = {0.01, 0.01, 0.01, 0.0008, 0.0008, 0.0008};
    const scalar tolerance = 1e-10;

    const Optimization::DiffCorrection::DiffCorrectionParameters<RightPart34> diffCorrParams =
            {{initState.u + delta, initState.t},
             vecW,
             derivativeSteps,
             intStep,
             tolerance};

    const std::pair<vector<scalar, RightPart34::size>, matrix<scalar, RightPart34::size, RightPart34::size>> pair =
                                                                                                 Optimization::DiffCorrection::Optimize<Integrators::ExplicitRK::ExplicitRungeKutta<Integrators::ExplicitRK::ButcherTableERK4>,
                                                                                                         RightPart34>(measurements,
                                                                                                                      params,
                                                                                                                      diffCorrParams);
    std::cout << "Calculated position:" << std::endl;
    std::cout << pair.first << std::endl;
    std::cout << "Covariance matrix:" << std::endl;
    std::cout << pair.second << std::endl;
}
