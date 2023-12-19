//
// Created by Арсений Плахотнюк on 20.10.2022.
//
#include "third_party/geographiclib/include/GeographicLib/GravityModel.hpp"
#include "third_party/geographiclib/include/GeographicLib/Math.hpp"
#include <fstream>
#include <iostream>
#include "gtest/gtest.h"

TEST(GEOLIB, TEST)
{

    auto gravity_model = GeographicLib::GravityModel("egm2008", "../../resources/gravity");
//    GeographicLib::Math::real lat = -90;
//    GeographicLib::Math::real lon = -180;
//    GeographicLib::Math::real h = 6800e3;
//
//    std::ofstream myfile("gravity_data.csv");
//    for(int i = 0; i < 181; ++i){
//        for(int j = 0; j < 361; ++j){
//            myfile.write()
//        }
//    }
//    gravity_model.Gravity();
    double h = 6800e3;
    GeographicLib::Math::real a = 0;
    std::ofstream out;
    out.open("gravity_data.txt");

    for (int lon = -180; lon < 180; lon++){
        for (int lat = -90; lat < 90; lat++){
            GeographicLib::Math::real gx, gy, gz;
            a = gravity_model.Gravity(lat, lon, h, gx, gy, gz);
            auto u = gravity_model.V(h, h, h, gx, gy, gz);
            out << lon << ' ' << lat << ' ' << u << '\n';
        }
    }
    ASSERT_TRUE(1);
}