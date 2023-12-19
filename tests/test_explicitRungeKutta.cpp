//
// Created by Арсений Плахотнюк on 09.03.2023.
//
#include "gtest/gtest.h"
#include "ballistics/integrators/explicitRungeKutta/ButcherTableERK4.hpp"
#include "ballistics/integrators/explicitRungeKutta/ExplicitRungeKutta.hpp"
#include "ballistics/integrators/RightPart/RightPartOfDiffEq.hpp"
#include "ballistics/orbit/Orbits.hpp"
#include "ballistics/forces/GravityForce.hpp"
#include "ballistics/frame_converter/FrameConverter.hpp"
#include "ballistics/time/Time.hpp"
#include "ballistics/time/TimeConverter.hpp"
#include "utility/EarthRotationResult.hpp"
#include "ballistics/atmosphere/Atmosphere.hpp"
#include "ballistics/ephemeris/Ephemeris.hpp"


TEST(CIRCLE, ORBIT){
    using namespace Ballistics;

//    Vector<scalar, size> u;
//    Time t;

    const scalar step = 10;
    const scalar tFinal = 40000;

    const scalar mu = 3.986004415e14;
    const scalar semimajor = 7000e3;
    const scalar ecc = 0;
    const scalar incl = 0;
    const scalar raan = 0;
    const scalar argPerig = 0;
    const scalar trueAnomaly = 0;

    Ballistics::Orbits::KeplerTrue KT = {incl, raan, argPerig, semimajor, ecc, trueAnomaly};

    Ballistics::Orbits::PosVel PV = Orbits::convert<Ballistics::Orbits::PosVel>(KT, mu);

    Integrators::RightPartPointPotential::State initState = {Vector<scalar, 6>{PV.r(0), PV.r(1), PV.r(2), PV.v(0), PV.v(1), PV.v(2)}, 0.};
    Integrators::RightPartPointPotential::Params params = {mu};
    Integrators::ExplicitRK::IntegrationParameters<Integrators::RightPartPointPotential, scalar> intParams = {step,
                                                                                                      tFinal};
    std::vector<Integrators::RightPartPointPotential::State> solution =
            Integrators::ExplicitRK::ExplicitRungeKutta<Integrators::ExplicitRK::ButcherTableERK4>::
            calc<Integrators::RightPartPointPotential>(initState,
                                                       params,
                                                       intParams);

    const std::string FILE_PATH = __FILE__;
    const std::string DIR_PATH = FILE_PATH.substr(0, FILE_PATH.size() - 27);

    std::fstream file;
    file.open(DIR_PATH + "data_files/integratorRk4.csv", std::ios::out);

    file << "t" << "," << "x" << "," << "y" << "," << "z" << "\n";

    for (const auto &el : solution)
    {
        file << el.t << "," << el.u[0] << "," << el.u[1] << "," << el.u[2] << "\n";
    }

    file.close();

}

TEST(REALPOTENTIAL, ORBIT){
    using namespace Ballistics;

//    Vector<scalar, size> u;
//    Time t;
    using TimeTt = Ballistics::Time::Time<Ballistics::Time::Scale::TT>;

    const scalar step = 100;
    const TimeTt tStart(2459842, 0.);
    const TimeTt tFinal (2459843, 0.);

    const scalar mu = 3.986004415e14;
    const scalar semimajor = 8500e3;
    const scalar ecc = 0;
    const scalar incl = 0.1;
    const scalar raan = 0;
    const scalar argPerig = 0;
    const scalar trueAnomaly = 0;
    const scalar mass = 1;

    Ballistics::Orbits::KeplerMean KT = {incl, raan, argPerig, semimajor, ecc, trueAnomaly};

    Ballistics::Orbits::PosVel PV = Orbits::convert<Ballistics::Orbits::PosVel>(KT, mu);

    Integrators::RightPartRealPotential::State initState = {Vector<scalar, 6>{PV.r(0), PV.r(1), PV.r(2), PV.v(0), PV.v(1), PV.v(2)}, tStart};

    const Ballistics::Frame::FrameConverter frameConverter("../../resources/rotation/earth_rotation.csv");
    const Ballistics::Time::TimeConverter timeConverter = Ballistics::Time::TimeConverter("../../resources/rotation/earth_rotation.csv");

    const Forces::GravityForce force("egm2008", "../../resources/gravity", 20, 20);

    Integrators::RightPartRealPotential::Params params{mass, force, frameConverter, timeConverter};
    Integrators::ExplicitRK::IntegrationParameters<Integrators::RightPartRealPotential, TimeTt> intParams = {step,
                                                                                                      tFinal};

    std::vector<Integrators::RightPartRealPotential::State> solution =
            Integrators::ExplicitRK::ExplicitRungeKutta<Integrators::ExplicitRK::ButcherTableERK4>::
            calc<Integrators::RightPartRealPotential>(initState,
                                                       params,
                                                       intParams);


    const std::string FILE_PATH = __FILE__;
    const std::string DIR_PATH = FILE_PATH.substr(0, FILE_PATH.size() - 27);

    std::fstream file;
    file.open(DIR_PATH + "data_files/integratorRk4RealPotential.csv", std::ios::out);

    std::cout << "write results to file" << std::endl;
    file << "t" << "," << "x" << "," << "y" << "," << "z" << "\n";

    for (const auto &el : solution)
    {
        file << el.t.toScalar() << "," << el.u[0] << "," << el.u[1] << "," << el.u[2] << "\n";
    }

    file.close();

}


TEST(ALLFORCES, ORBIT){
    using namespace Ballistics;

//    Vector<scalar, size> u;
//    Time t;
    using TimeTt = Ballistics::Time::Time<Ballistics::Time::Scale::TT>;

    const scalar step = 100;
    const TimeTt tStart(2459842, 0.);
    const TimeTt tFinal(2459844, 0.);

    const scalar mu = 3.986004415e14;
    const scalar semimajor = 27500e3;
    const scalar ecc = 0.6;
    const scalar incl = M_PI * 63 / 180;
    const scalar raan = 0.;
    const scalar argPerig = 0.;
    const scalar trueAnomaly = 0.;
    const scalar mass = 1e3;
    const scalar S = 1;
    const scalar solarRadius = 695.700e+6;
    const scalar earthRadius = 6.371e+6;
    const scalar TSI = 1366;
    const scalar AU = 149597870700;
    const scalar Cd = 2.2;

    Ballistics::Orbits::KeplerMean KT = {incl, raan, argPerig, semimajor, ecc, trueAnomaly};

    Ballistics::Orbits::PosVel PV = Orbits::convert<Ballistics::Orbits::PosVel>(KT, mu);

    Integrators::RightPartAllForces::State initState = {Vector<scalar, 6>{PV.r(0), PV.r(1), PV.r(2), PV.v(0), PV.v(1), PV.v(2)}, tStart};

    const Ballistics::Frame::FrameConverter frameConverter("../../resources/rotation/earth_rotation.csv");
    const Ballistics::Time::TimeConverter timeConverter = Ballistics::Time::TimeConverter("../../resources/rotation/earth_rotation.csv");

    const Forces::GravityForce force("egm2008", "../../resources/gravity", 20, 20);

    Integrators::RightPartAllForces::Params params{mass,
                                                   force,
                                                   frameConverter,
                                                   timeConverter,
                                                   solarRadius,
                                                   earthRadius,
                                                   S,
                                                   TSI,
                                                   AU,
                                                   {Ephemeris::Body::Sun, Ephemeris::Body::MarsBarycenter,
                                                    Ephemeris::Body::Moon, Ephemeris::Body::JupiterBarycenter,
                                                    Ephemeris::Body::SaturnBarycenter},
                                                   {1.3271244e+20, 4.2828e+13,
                                                    4.9028e+12, 1.26686534e+17,
                                                    3.7931187e+16},
                                                   Ephemeris::Ephemeris(),
                                                   Atmosphere::Atmosphere("../../resources/atmosphere/SW-Last5Years.csv")


    };
    Integrators::ExplicitRK::IntegrationParameters<Integrators::RightPartAllForces, TimeTt> intParams = {step,
                                                                                                             tFinal};

    std::vector<Integrators::RightPartAllForces::State> solution =
            Integrators::ExplicitRK::ExplicitRungeKutta<Integrators::ExplicitRK::ButcherTableERK4>::
            calc<Integrators::RightPartAllForces>(initState,
                                                      params,
                                                      intParams);


    const std::string FILE_PATH = __FILE__;
    const std::string DIR_PATH = FILE_PATH.substr(0, FILE_PATH.size() - 27);

    std::fstream file;
    file.open(DIR_PATH + "data_files/integratorRk4AllForces.csv", std::ios::out);

    std::cout << "write results to file" << std::endl;
    file << "t" << "," << "x" << "," << "y" << "," << "z" << "\n";

    for (const auto &el : solution)
    {
        file << el.t.toScalar() << "," << el.u[0] << "," << el.u[1] << "," << el.u[2] << "\n";
    }

    file.close();

}
