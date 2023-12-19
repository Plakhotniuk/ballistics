//
// Created by Арсений Плахотнюк on 20.12.2022.
//
#include "../src/ballistics/atmosphere/Atmosphere.hpp"
#include <iostream>
#include <vector>
#include "gtest/gtest.h"

TEST(ATMOSPHERE1, TEST1){
    const std::string FILE_PATH = __FILE__;
    const std::string dir_path = FILE_PATH.substr(0, FILE_PATH.size() - 19);
    std::string filename = dir_path + "utility/gost_density.csv";

    Ballistics::Atmosphere::Atmosphere atmosphere = Ballistics::Atmosphere::Atmosphere("../../resources/atmosphere/SW-Last5Years.csv");
    rapidcsv::Document doc(filename);
    std::vector<scalar> xCol = doc.GetColumn<scalar>("x");
    std::vector<scalar> yCol = doc.GetColumn<scalar>("y");
    std::vector<scalar> zCol = doc.GetColumn<scalar>("z");
    std::vector<scalar> rhoCol = doc.GetColumn<scalar>("rho");

    integer n = xCol.size();
    scalar tolerance = 7e-2;

    for (int i = 0; i < n ; ++i)
    {
        ASSERT_TRUE((abs(atmosphere.getDensity({xCol[i], yCol[i], zCol[i]}, {2459277.5, 0.0}) - rhoCol[i]) / rhoCol[i]) <= tolerance);
    }
}