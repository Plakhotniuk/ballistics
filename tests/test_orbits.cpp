//
// Created by Арсений Плахотнюк on 10.02.2023.
//
#include "gtest/gtest.h"
#include "ballistics/orbit/Orbits.hpp"
#include "ballistics/types/BasicTypes.hpp"

using namespace Ballistics;

TEST(ELEMENTS1, CONVERTATIONS2PV){
    scalar tolerance = 1.e-14;
    std::string const csvFileLocation = "../../resources/orbits/KeplerianAndCartesian.csv";
    io::CSVReader<14> data(csvFileLocation);
    scalar mu, x, y, z, vx, vy, vz, i, RAAN, a, e, argPer, trueAnom, meanAnom;

    data.read_header(io::ignore_extra_column, "mu", "x", "y", "z", "vx", "vy", "vz",
                     "i", "RAAN", "a", "e", "argPer", "trueAnom", "meanAnom");
    int k = 1;
    while (data.read_row(mu, x, y, z, vx, vy, vz, i, RAAN, a, e, argPer, trueAnom, meanAnom)) {

        Ballistics::Orbits::PosVel PV = {{x, y, z}, {vx, vy, vz}};
        Ballistics::Orbits::KeplerMean KM = {i, RAAN, argPer, a, e, meanAnom};
        Ballistics::Orbits::KeplerTrue KT = {i, RAAN, argPer, a, e, trueAnom};

        Ballistics::Orbits::PosVel ConvertedFromKM = Ballistics::Orbits::convert<Ballistics::Orbits::PosVel>(KM, mu);
        Ballistics::Orbits::PosVel ConvertedFromKT = Ballistics::Orbits::convert<Ballistics::Orbits::PosVel>(KT, mu);

        ASSERT_NEAR(std::abs((PV.r.norm() - ConvertedFromKM.r.norm()) / ConvertedFromKM.r.norm()), 0, tolerance);
        ASSERT_NEAR(std::abs((PV.v.norm() - ConvertedFromKM.v.norm()) / ConvertedFromKM.v.norm()), 0, tolerance);

        ASSERT_NEAR(std::abs((PV.r.norm() - ConvertedFromKT.r.norm()) / ConvertedFromKT.r.norm()), 0, tolerance);
        ASSERT_NEAR(std::abs((PV.v.norm() - ConvertedFromKT.v.norm()) / ConvertedFromKT.v.norm()), 0, tolerance);

    }
}

bool compareAngles2Pi(const scalar ang1, const scalar ang2, const scalar tolerance=1.e-14){
    scalar temp = std::abs(std::abs(ang1 - ang2) - static_cast<scalar>(2)*M_PI);
    return (std::abs(ang1 - ang2) < tolerance || std::abs(std::abs(ang1 - ang2) - static_cast<scalar>(2)*M_PI) < tolerance);
}

TEST(ELEMENTS2, CONVERTATIONS2KM_ANGLES){
    scalar tolerance = 1.e-14;
    std::string const csvFileLocation = "../../resources/orbits/KeplerianAndCartesian.csv";
    io::CSVReader<14> data(csvFileLocation);
    scalar mu, x, y, z, vx, vy, vz, i, RAAN, a, e, argPer, trueAnom, meanAnom;

    data.read_header(io::ignore_extra_column, "mu", "x", "y", "z", "vx", "vy", "vz",
                     "i", "RAAN", "a", "e", "argPer", "trueAnom", "meanAnom");

    while (data.read_row(mu, x, y, z, vx, vy, vz, i, RAAN, a, e, argPer, trueAnom, meanAnom)) {

            Ballistics::Orbits::KeplerMean KM = {i, RAAN, argPer, a, e, meanAnom};
            Ballistics::Orbits::PosVel PV = {{x, y, z}, {vx, vy, vz}};
            Ballistics::Orbits::KeplerTrue KT = {i, RAAN, argPer, a, e, trueAnom};

            Ballistics::Orbits::KeplerMean ConvertedFromPV = Ballistics::Orbits::convert<Ballistics::Orbits::KeplerMean>(PV, mu);
            Ballistics::Orbits::KeplerMean ConvertedFromKT = Ballistics::Orbits::convert<Ballistics::Orbits::KeplerMean>(KT, mu);

            // долгота восх, средняя и истинная аном, аргумент перицентра
            ASSERT_TRUE(compareAngles2Pi(KM.inclination, ConvertedFromPV.inclination, tolerance));
            ASSERT_TRUE(compareAngles2Pi(KM.ascendingNodeLongitude, ConvertedFromPV.ascendingNodeLongitude, tolerance));
            ASSERT_TRUE(compareAngles2Pi(KM.inclination, ConvertedFromKT.inclination, tolerance));
            ASSERT_TRUE(compareAngles2Pi(KM.ascendingNodeLongitude, ConvertedFromKT.ascendingNodeLongitude, tolerance));
    }
}

TEST(ELEMENTS3, CONVERTATIONS2KM_SEMIMAJORAXIS) {
    scalar tolerance = 1.e-14;
    std::string const csvFileLocation = "../../resources/orbits/KeplerianAndCartesian.csv";
    io::CSVReader<14> data(csvFileLocation);
    scalar mu, x, y, z, vx, vy, vz, i, RAAN, a, e, argPer, trueAnom, meanAnom;

    data.read_header(io::ignore_extra_column, "mu", "x", "y", "z", "vx", "vy", "vz",
                     "i", "RAAN", "a", "e", "argPer", "trueAnom", "meanAnom");
    while (data.read_row(mu, x, y, z, vx, vy, vz, i, RAAN, a, e, argPer, trueAnom, meanAnom)) {

        Ballistics::Orbits::KeplerMean KM = {i, RAAN, argPer, a, e, meanAnom};
        Ballistics::Orbits::PosVel PV = {{x,  y,  z},
                                         {vx, vy, vz}};
        Ballistics::Orbits::KeplerTrue KT = {i, RAAN, argPer, a, e, trueAnom};

        Ballistics::Orbits::KeplerMean ConvertedFromPV = Ballistics::Orbits::convert<Ballistics::Orbits::KeplerMean>(PV,
                                                                                                                     mu);
        Ballistics::Orbits::KeplerMean ConvertedFromKT = Ballistics::Orbits::convert<Ballistics::Orbits::KeplerMean>(KT,
                                                                                                                     mu);

        ASSERT_NEAR(std::abs((KM.semiMajorAxis - ConvertedFromPV.semiMajorAxis)) / std::abs(ConvertedFromPV.semiMajorAxis), 0, tolerance);

        ASSERT_NEAR(std::abs((KM.semiMajorAxis - ConvertedFromKT.semiMajorAxis)) / std::abs(ConvertedFromKT.semiMajorAxis), 0, tolerance);
    }
}

TEST(ELEMENTS4, CONVERTATIONS2KM_ECCENTRICITY) {
    scalar tolerance = 1.e-15;
    std::string const csvFileLocation = "../../resources/orbits/KeplerianAndCartesian.csv";
    io::CSVReader<14> data(csvFileLocation);
    scalar mu, x, y, z, vx, vy, vz, i, RAAN, a, e, argPer, trueAnom, meanAnom;

    data.read_header(io::ignore_extra_column, "mu", "x", "y", "z", "vx", "vy", "vz",
                     "i", "RAAN", "a", "e", "argPer", "trueAnom", "meanAnom");
    while (data.read_row(mu, x, y, z, vx, vy, vz, i, RAAN, a, e, argPer, trueAnom, meanAnom)) {

        Ballistics::Orbits::KeplerMean KM = {i, RAAN, argPer, a, e, meanAnom};
        Ballistics::Orbits::PosVel PV = {{x,  y,  z},
                                         {vx, vy, vz}};
        Ballistics::Orbits::KeplerTrue KT = {i, RAAN, argPer, a, e, trueAnom};

        Ballistics::Orbits::KeplerMean KMConvertedFromPV = Ballistics::Orbits::convert<Ballistics::Orbits::KeplerMean>(PV,
                                                                                                                     mu);
        Ballistics::Orbits::KeplerMean KMConvertedFromKT = Ballistics::Orbits::convert<Ballistics::Orbits::KeplerMean>(KT,
                                                                                                                     mu);
        ASSERT_NEAR(KM.eccentricity, KMConvertedFromPV.eccentricity, tolerance);
        ASSERT_NEAR(KM.eccentricity, KMConvertedFromKT.eccentricity, tolerance);
    }
}


TEST(ELEMENTS5, CONVERTATIONS2KT_ANGLES){
    scalar tolerance = 1.e-15;
    std::string const csvFileLocation = "../../resources/orbits/KeplerianAndCartesian.csv";
    io::CSVReader<14> data(csvFileLocation);
    scalar mu, x, y, z, vx, vy, vz, i, RAAN, a, e, argPer, trueAnom, meanAnom;

    data.read_header(io::ignore_extra_column, "mu", "x", "y", "z", "vx", "vy", "vz",
                     "i", "RAAN", "a", "e", "argPer", "trueAnom", "meanAnom");
    while (data.read_row(mu, x, y, z, vx, vy, vz, i, RAAN, a, e, argPer, trueAnom, meanAnom)) {

        Ballistics::Orbits::KeplerTrue KT = {i, RAAN, argPer, a, e, trueAnom};
        Ballistics::Orbits::KeplerMean KM = {i, RAAN, argPer, a, e, meanAnom};
        Ballistics::Orbits::PosVel PV = {{x, y, z}, {vx, vy, vz}};

        Ballistics::Orbits::KeplerTrue ConvertedFromPV = Ballistics::Orbits::convert<Ballistics::Orbits::KeplerTrue>(PV, mu);
        Ballistics::Orbits::KeplerTrue ConvertedFromKM = Ballistics::Orbits::convert<Ballistics::Orbits::KeplerTrue>(KM, mu);

        ASSERT_NEAR(KT.inclination, ConvertedFromPV.inclination, tolerance);
        ASSERT_TRUE(compareAngles2Pi(KT.ascendingNodeLongitude, ConvertedFromPV.ascendingNodeLongitude, tolerance));
//        ASSERT_TRUE(compareAngles2Pi(KT.periapsisArgument + KT.trueAnomaly + KT.ascendingNodeLongitude, ConvertedFromPV.trueAnomaly + ConvertedFromPV.periapsisArgument + ConvertedFromPV.ascendingNodeLongitude, tolerance));
//        ASSERT_NEAR(KT.trueAnomaly, ConvertedFromPV.trueAnomaly, tolerance);

        ASSERT_NEAR(KT.inclination, ConvertedFromKM.inclination, tolerance);
        ASSERT_TRUE(compareAngles2Pi(KT.ascendingNodeLongitude, ConvertedFromKM.ascendingNodeLongitude, tolerance));
//        ASSERT_NEAR(KT.periapsisArgument, ConvertedFromKM.periapsisArgument, tolerance);
//        ASSERT_NEAR(KT.trueAnomaly, ConvertedFromKM.trueAnomaly, tolerance);

    }
}

TEST(ELEMENTS6, CONVERTATIONS2KT_SEMIMAJORAXIS) {
    scalar tolerance = 1.e-14;
    std::string const csvFileLocation = "../../resources/orbits/KeplerianAndCartesian.csv";
    io::CSVReader<14> data(csvFileLocation);
    scalar mu, x, y, z, vx, vy, vz, i, RAAN, a, e, argPer, trueAnom, meanAnom;

    data.read_header(io::ignore_extra_column, "mu", "x", "y", "z", "vx", "vy", "vz",
                     "i", "RAAN", "a", "e", "argPer", "trueAnom", "meanAnom");
//    int k = 1;
    while (data.read_row(mu, x, y, z, vx, vy, vz, i, RAAN, a, e, argPer, trueAnom, meanAnom)) {

        Ballistics::Orbits::KeplerTrue KT = {i, RAAN, argPer, a, e, trueAnom};
        Ballistics::Orbits::KeplerMean KM = {i, RAAN, argPer, a, e, meanAnom};
        Ballistics::Orbits::PosVel PV = {{x,  y,  z},
                                         {vx, vy, vz}};

        Ballistics::Orbits::KeplerTrue ConvertedFromPV = Ballistics::Orbits::convert<Ballistics::Orbits::KeplerTrue>(PV, mu);
        Ballistics::Orbits::KeplerTrue ConvertedFromKM = Ballistics::Orbits::convert<Ballistics::Orbits::KeplerTrue>(KM, mu);

        ASSERT_NEAR(KT.semiMajorAxis, ConvertedFromPV.semiMajorAxis, tolerance*KT.semiMajorAxis);

        ASSERT_NEAR(KT.semiMajorAxis, ConvertedFromKM.semiMajorAxis, tolerance*KT.semiMajorAxis);

    }
}

TEST(ELEMENTS7, CONVERTATIONS2KT_ECCENTRICITY) {
    scalar tolerance = 1.e-15;
    std::string const csvFileLocation = "../../resources/orbits/KeplerianAndCartesian.csv";
    io::CSVReader<14> data(csvFileLocation);
    scalar mu, x, y, z, vx, vy, vz, i, RAAN, a, e, argPer, trueAnom, meanAnom;

    data.read_header(io::ignore_extra_column, "mu", "x", "y", "z", "vx", "vy", "vz",
                     "i", "RAAN", "a", "e", "argPer", "trueAnom", "meanAnom");

    while (data.read_row(mu, x, y, z, vx, vy, vz, i, RAAN, a, e, argPer, trueAnom, meanAnom)) {

        Ballistics::Orbits::KeplerTrue KT = {i, RAAN, argPer, a, e, trueAnom};
        Ballistics::Orbits::KeplerMean KM = {i, RAAN, argPer, a, e, meanAnom};
        Ballistics::Orbits::PosVel PV = {{x,  y,  z},
                                         {vx, vy, vz}};

        Ballistics::Orbits::KeplerTrue ConvertedFromPV = Ballistics::Orbits::convert<Ballistics::Orbits::KeplerTrue>(PV, mu);
        Ballistics::Orbits::KeplerTrue ConvertedFromKM = Ballistics::Orbits::convert<Ballistics::Orbits::KeplerTrue>(KM, mu);

        ASSERT_NEAR(KT.eccentricity, ConvertedFromPV.eccentricity, tolerance);

        ASSERT_NEAR(KT.eccentricity, ConvertedFromKM.eccentricity, tolerance);
    }
}


