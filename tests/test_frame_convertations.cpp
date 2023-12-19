//
// Created by Арсений Плахотнюк on 08.12.2022.
//

#include "gtest/gtest.h"
#include "ballistics/frame_converter/FrameConverter.hpp"
#include "ballistics/time/Time.hpp"
#include "ballistics/time/TimeConverter.hpp"
#include "utility/EarthRotationResult.hpp"

TEST(FRAME, CONVERT){
    using namespace Ballistics::Time;
    using namespace Ballistics;
    struct Line
    {
        scalar mjd;
        Eigen::Vector3d v1;
        Eigen::Vector3d v2;
        Eigen::Vector3d v3;
    };

    Ballistics::Frame::FrameConverter FrameConverter("../../resources/rotation/earth_rotation.csv");
    Ballistics::Time::TimeConverter TimeConverter = Ballistics::Time::TimeConverter("../../resources/rotation/earth_rotation.csv");
    Eigen::Vector3d initVec1 = {6.7e+6, 0., 0.};
    Eigen::Vector3d initVec2 = {0., 6.7e+6, 0.};
    Eigen::Vector3d initVec3 = {0., 0., 6.7e+6};
    Ballistics::scalar tolerance = 1.e-10;

    for (const auto el : earthRotationResult)
    {
        Line line = {el[0], {el[1], el[2], el[3]}, {el[4], el[5], el[6]}, {el[7], el[8], el[9]}};
        auto const timeTT = Ballistics::Time::Time<Ballistics::Time::Scale::TT>(line.mjd, 0.);
        auto const timeUTC = TimeConverter.convert<Ballistics::Time::Scale::UTC>(timeTT).value();
        auto const timeUT1 = TimeConverter.convert<Ballistics::Time::Scale::UT1>(timeTT).value();

        Eigen::Quaterniond Q = FrameConverter.calcQuaternion<Ballistics::Frame::Frame::GCRS, Ballistics::Frame::Frame::ITRS>(timeTT,
                                                                                                     timeUT1,
                                                                                                     timeUTC);

        Eigen::Vector3d vec1 = Q * initVec1;
        Eigen::Vector3d vec2 = Q * initVec2;
        Eigen::Vector3d vec3 = Q * initVec3;

        ASSERT_NEAR(vec1.norm(), line.v1.norm(), tolerance);
        ASSERT_NEAR(vec2.norm(), line.v2.norm(), tolerance);
        ASSERT_NEAR(vec3.norm(), line.v3.norm(), tolerance);

    }
}
