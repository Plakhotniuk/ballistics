//
// Created by Арсений Плахотнюк on 20.12.2022.
//
#include "ballistics/time/Time.hpp"
#include "ballistics/time/TimeConverter.hpp"
#include "ballistics/types/BasicTypes.hpp"
#include "gtest/gtest.h"
#include "third_party/include/rapidcsv.h"
#include "ballistics/ephemeris/Ephemeris.hpp"

bool isEqualByRelativeError(const Eigen::Vector3d &l, const Eigen::Vector3d &r, const scalar tolerance)
{
    return (l - r).norm()  <= tolerance * r.norm();
}


TEST(EPHEMERIS, TEST1){
    const std::string FILE_PATH = __FILE__;
    const std::string dir_path = FILE_PATH.substr(0, FILE_PATH.size() - 18);
    std::string filename = dir_path + "utility/ephems.csv";

    rapidcsv::Document doc(filename);
    std::vector<integer> bodyNumCol = doc.GetColumn<integer>("body");
    std::vector<scalar> jdDayCol = doc.GetColumn<scalar>("jdDay");
    std::vector<scalar> jdPartCol = doc.GetColumn<scalar>("jdDayPart");
    std::vector<scalar> xCol = doc.GetColumn<scalar>("x");
    std::vector<scalar> yCol = doc.GetColumn<scalar>("y");
    std::vector<scalar> zCol = doc.GetColumn<scalar>("z");
    std::vector<scalar> vxCol = doc.GetColumn<scalar>("vx");
    std::vector<scalar> vyCol = doc.GetColumn<scalar>("vy");
    std::vector<scalar> vzCol = doc.GetColumn<scalar>("vz");


    auto converter = Ballistics::Time::TimeConverter("../../resources/rotation/earth_rotation.csv");
    Ballistics::Ephemeris::Ephemeris Ephems;

    integer n = xCol.size();
    scalar tolerance = 1e-14;

    for (int i = 0; i < n; ++i)
    {
        const Ballistics::Time::Time<Ballistics::Time::Scale::TT> timeTT(jdDayCol[i], jdPartCol[i]);
        const Ballistics::Time::Time<Ballistics::Time::Scale::TDB> timeTDB = converter.convert<Ballistics::Time::Scale::TDB, Ballistics::Time::Scale::TT>(timeTT).value();
        const Vector3d position = {xCol[i], yCol[i], zCol[i]};
        const Vector3d velocity = {vxCol[i], vyCol[i], vzCol[i]};

        ASSERT_TRUE(isEqualByRelativeError(Ephems.calcCoordinates(timeTDB, bodyNumCol[i], 3)*1000, position, tolerance));

        ASSERT_TRUE(isEqualByRelativeError(Ephems.calcVelocities(timeTDB, bodyNumCol[i], 3)*1000, velocity, tolerance));
    }

}
