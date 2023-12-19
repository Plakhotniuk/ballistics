//
// Created by Арсений Плахотнюк on 01.10.2022.
//
#include "gtest/gtest.h"
#include "ballistics/time/Time.hpp"
#include "ballistics/time/TimeConverter.hpp"
#include "ballistics/types/BasicTypes.hpp"
#include "tests/utility/TimeResult.hpp"
#include "g3log/logworker.hpp"
#include "tests/log_files/logger.hpp"

using namespace Ballistics;
// UT1->
TEST(CONVERTATION_UT1_UTC, TIME){
    INITIALIZE_LOG(27);
    scalar tol = 1e-15;
    const std::string dir_path = FILE_PATH.substr(0, FILE_PATH.size() - 27);
    std::string filename = dir_path + "/test_data/earth_rotation.csv";
    Ballistics::Time::TimeConverter converter(filename);
    for(int i = 0; i < 500; ++i)
    {
        Ballistics::Time::Time<Ballistics::Time::Scale::UTC> time(timeResult[i][3], timeResult[i][4]);
        Expected<Ballistics::Time::Time<Ballistics::Time::Scale::UT1>, Error> time_converted = converter.convert<Ballistics::Time::Scale::UT1, Ballistics::Time::Scale::UTC>(time);
        Ballistics::Time::Time<Ballistics::Time::Scale::UT1> time2compare(timeResult[i][1], timeResult[i][2]);
        ASSERT_TRUE(time_converted.has_value());
        ASSERT_NEAR(time2compare.get_jdDay(), time_converted.value().get_jdDay(), tol);
        ASSERT_NEAR(time2compare.get_jdPart(), time_converted.value().get_jdPart(), tol);
    }

}

TEST(CONVERTATION_UT1_TT, TIME){
    INITIALIZE_LOG(27)
    scalar tol = 1e-15;
    const std::string dir_path = FILE_PATH.substr(0, FILE_PATH.size() - 27);
    std::string filename = dir_path + "/test_data/earth_rotation.csv";
    Ballistics::Time::TimeConverter converter(filename);
    for(int i = 0; i < 500; ++i)
    {
        Ballistics::Time::Time<Ballistics::Time::Scale::UT1> time(timeResult[i][1], timeResult[i][2]);
        Expected<Ballistics::Time::Time<Ballistics::Time::Scale::TT>, Error> time_converted =
                converter.convert<Ballistics::Time::Scale::TT, Ballistics::Time::Scale::UT1>(time);
        Ballistics::Time::Time<Ballistics::Time::Scale::TT> time2compare(timeResult[i][7], timeResult[i][8]);
        ASSERT_TRUE(time_converted.has_value());
        ASSERT_TRUE(time2compare == time_converted.value());
    }

}

TEST(CONVERTATION_UT1_TAI, TIME){
    INITIALIZE_LOG(27);
    scalar tol = 1e-15;
    const std::string dir_path = FILE_PATH.substr(0, FILE_PATH.size() - 27);
    std::string filename = dir_path + "/test_data/earth_rotation.csv";
    Ballistics::Time::TimeConverter converter(filename);
    for(int i = 0; i < 500; ++i)
    {
        Ballistics::Time::Time<Ballistics::Time::Scale::UT1> time(timeResult[i][1], timeResult[i][2]);
        Expected<Ballistics::Time::Time<Ballistics::Time::Scale::TAI>, Error> time_converted =
                converter.convert<Ballistics::Time::Scale::TAI, Ballistics::Time::Scale::UT1>(time);
        Ballistics::Time::Time<Ballistics::Time::Scale::TAI> time2compare(timeResult[i][5], timeResult[i][6]);
        ASSERT_TRUE(time_converted.has_value());
        ASSERT_TRUE(time2compare == time_converted.value());
    }

}

TEST(CONVERTATION_UT1_TCB, TIME){
    INITIALIZE_LOG(27);
    scalar tol = 1e-15;
    const std::string dir_path = FILE_PATH.substr(0, FILE_PATH.size() - 27);
    std::string filename = dir_path + "/test_data/earth_rotation.csv";
    Ballistics::Time::TimeConverter converter(filename);
    for(int i = 0; i < 500; ++i)
    {
        Ballistics::Time::Time<Ballistics::Time::Scale::UT1> time(timeResult[i][1], timeResult[i][2]);
        Expected<Ballistics::Time::Time<Ballistics::Time::Scale::TCB>, Error> time_converted =
                converter.convert<Ballistics::Time::Scale::TCB, Ballistics::Time::Scale::UT1>(time);
        Ballistics::Time::Time<Ballistics::Time::Scale::TCB> time2compare(timeResult[i][11], timeResult[i][12]);
        ASSERT_TRUE(time_converted.has_value());
        ASSERT_TRUE(time2compare == time_converted.value());
    }

}

TEST(CONVERTATION_UT1_TCG, TIME){
    INITIALIZE_LOG(27);
    scalar tol = 1e-15;
    const std::string dir_path = FILE_PATH.substr(0, FILE_PATH.size() - 27);
    std::string filename = dir_path + "/test_data/earth_rotation.csv";
    Ballistics::Time::TimeConverter converter(filename);
    for(int i = 0; i < 500; ++i)
    {
        Ballistics::Time::Time<Ballistics::Time::Scale::UT1> time(timeResult[i][1], timeResult[i][2]);
        Expected<Ballistics::Time::Time<Ballistics::Time::Scale::TCG>, Error> time_converted =
                converter.convert<Ballistics::Time::Scale::TCG, Ballistics::Time::Scale::UT1>(time);
        Ballistics::Time::Time<Ballistics::Time::Scale::TCG> time2compare(timeResult[i][9], timeResult[i][10]);
        ASSERT_TRUE(time_converted.has_value());
        ASSERT_TRUE(time2compare == time_converted.value());
    }

}

// UTC->
TEST(CONVERTATION_UTC_TT, TIME){
    INITIALIZE_LOG(27);
    scalar tol = 1e-15;
    const std::string dir_path = FILE_PATH.substr(0, FILE_PATH.size() - 27);
    std::string filename = dir_path + "/test_data/earth_rotation.csv";
    Ballistics::Time::TimeConverter converter(filename);
    for(int i = 0; i < 500; ++i)
    {
        Ballistics::Time::Time<Ballistics::Time::Scale::UTC> time(timeResult[i][3], timeResult[i][4]);
        Expected<Ballistics::Time::Time<Ballistics::Time::Scale::TT>, Error> time_converted = converter.convert<Ballistics::Time::Scale::TT, Ballistics::Time::Scale::UTC>(time);
        Ballistics::Time::Time<Ballistics::Time::Scale::TT> time2compare(timeResult[i][7], timeResult[i][8]);
        ASSERT_TRUE(time_converted.has_value());
        ASSERT_TRUE(time2compare == time_converted.value());
    }

}

TEST(CONVERTATION_UTC_TAI, TIME){
    INITIALIZE_LOG(27);
    scalar tol = 1e-15;
    const std::string dir_path = FILE_PATH.substr(0, FILE_PATH.size() - 27);
    std::string filename = dir_path + "/test_data/earth_rotation.csv";
    Ballistics::Time::TimeConverter converter(filename);
    std::cout<<std::numeric_limits<double>::epsilon();
    for(int i = 0; i < 500; ++i)
    {
        Ballistics::Time::Time<Ballistics::Time::Scale::UTC> time(timeResult[i][3], timeResult[i][4]);
        Expected<Ballistics::Time::Time<Ballistics::Time::Scale::TAI>, Error> time_converted =
                converter.convert<Ballistics::Time::Scale::TAI, Ballistics::Time::Scale::UTC>(time);
        Ballistics::Time::Time<Ballistics::Time::Scale::TAI> time2compare(timeResult[i][5], timeResult[i][6]);
        ASSERT_TRUE(time_converted.has_value());
        ASSERT_TRUE(time2compare == time_converted.value());
    }

}

TEST(CONVERTATION_UTC_UT1, TIME){
    INITIALIZE_LOG(27);
    scalar tol = 1e-15;
    const std::string dir_path = FILE_PATH.substr(0, FILE_PATH.size() - 27);
    std::string filename = dir_path + "/test_data/earth_rotation.csv";
    Ballistics::Time::TimeConverter converter(filename);
    for(int i = 0; i < 500; ++i)
    {
        Ballistics::Time::Time<Ballistics::Time::Scale::UTC> time(timeResult[i][3], timeResult[i][4]);
        Expected<Ballistics::Time::Time<Ballistics::Time::Scale::UT1>, Error> time_converted =
                converter.convert<Ballistics::Time::Scale::UT1, Ballistics::Time::Scale::UTC>(time);
        Ballistics::Time::Time<Ballistics::Time::Scale::UT1> time2compare(timeResult[i][1], timeResult[i][2]);
        ASSERT_TRUE(time_converted.has_value());
        ASSERT_TRUE(time2compare == time_converted.value());
    }

}

TEST(CONVERTATION_UTC_TCB, TIME){
    INITIALIZE_LOG(27);
    scalar tol = 1e-15;
    const std::string dir_path = FILE_PATH.substr(0, FILE_PATH.size() - 27);
    std::string filename = dir_path + "/test_data/earth_rotation.csv";
    Ballistics::Time::TimeConverter converter(filename);
    for(int i = 0; i < 500; ++i)
    {
        Ballistics::Time::Time<Ballistics::Time::Scale::UTC> time(timeResult[i][3], timeResult[i][4]);
        Expected<Ballistics::Time::Time<Ballistics::Time::Scale::TCB>, Error> time_converted =
                converter.convert<Ballistics::Time::Scale::TCB, Ballistics::Time::Scale::UTC>(time);
        Ballistics::Time::Time<Ballistics::Time::Scale::TCB> time2compare(timeResult[i][11], timeResult[i][12]);
        ASSERT_TRUE(time_converted.has_value());
        ASSERT_TRUE(time2compare == time_converted.value());
    }

}

TEST(CONVERTATION_UTC_TDB, TIME){
    INITIALIZE_LOG(27);
    scalar tol = 1e-15;
    const std::string dir_path = FILE_PATH.substr(0, FILE_PATH.size() - 27);
    std::string filename = dir_path + "/test_data/earth_rotation.csv";
    Ballistics::Time::TimeConverter converter(filename);
    for(int i = 0; i < 500; ++i)
    {
        Ballistics::Time::Time<Ballistics::Time::Scale::UTC> time(timeResult[i][3], timeResult[i][4]);
        Expected<Ballistics::Time::Time<Ballistics::Time::Scale::TDB>, Error> time_converted =
                converter.convert<Ballistics::Time::Scale::TDB, Ballistics::Time::Scale::UTC>(time);
        Ballistics::Time::Time<Ballistics::Time::Scale::TDB> time2compare(timeResult[i][13], timeResult[i][14]);
        ASSERT_TRUE(time_converted.has_value());
        ASSERT_TRUE(time2compare == time_converted.value());
    }

}

// TCG->
TEST(CONVERTATION_TCG_TT, TIME){
    INITIALIZE_LOG(27);
    scalar tol = 1e-15;
    const std::string dir_path = FILE_PATH.substr(0, FILE_PATH.size() - 27);
    std::string filename = dir_path + "/test_data/earth_rotation.csv";
    Ballistics::Time::TimeConverter converter(filename);
    for(int i = 0; i < 500; ++i) {
        Ballistics::Time::Time<Ballistics::Time::Scale::TCG> time(timeResult[i][9], timeResult[i][10]);
        Expected<Ballistics::Time::Time<Ballistics::Time::Scale::TT>, Error> time_converted = converter.convert<Ballistics::Time::Scale::TT, Ballistics::Time::Scale::TCG>(
                time);
        Ballistics::Time::Time<Ballistics::Time::Scale::TT> time2compare(timeResult[i][7], timeResult[i][8]);
        ASSERT_TRUE(time_converted.has_value());
        ASSERT_TRUE(time2compare == time_converted.value());
    }
}

TEST(CONVERTATION_TCG_TDB, TIME){
    INITIALIZE_LOG(27);
    scalar tol = 1e-15;
    const std::string dir_path = FILE_PATH.substr(0, FILE_PATH.size() - 27);
    std::string filename = dir_path + "/test_data/earth_rotation.csv";
    Ballistics::Time::TimeConverter converter(filename);
    for(int i = 0; i < 500; ++i) {
        Ballistics::Time::Time<Ballistics::Time::Scale::TCG> time(timeResult[i][9], timeResult[i][10]);
        Expected<Ballistics::Time::Time<Ballistics::Time::Scale::TDB>, Error> time_converted =
                converter.convert<Ballistics::Time::Scale::TDB, Ballistics::Time::Scale::TCG>(time);
        Ballistics::Time::Time<Ballistics::Time::Scale::TDB> time2compare(timeResult[i][13], timeResult[i][14]);
        ASSERT_TRUE(time_converted.has_value());
        ASSERT_TRUE(time2compare == time_converted.value());
    }
}

TEST(CONVERTATION_TCG_UT1, TIME){
    INITIALIZE_LOG(27);
    scalar tol = 1e-15;
    const std::string dir_path = FILE_PATH.substr(0, FILE_PATH.size() - 27);
    std::string filename = dir_path + "/test_data/earth_rotation.csv";
    Ballistics::Time::TimeConverter converter(filename);
    for(int i = 0; i < 500; ++i) {
        Ballistics::Time::Time<Ballistics::Time::Scale::TCG> time(timeResult[i][9], timeResult[i][10]);
        Expected<Ballistics::Time::Time<Ballistics::Time::Scale::UT1>, Error> time_converted =
                converter.convert<Ballistics::Time::Scale::UT1, Ballistics::Time::Scale::TCG>(time);
        Ballistics::Time::Time<Ballistics::Time::Scale::UT1> time2compare(timeResult[i][1], timeResult[i][2]);
        ASSERT_TRUE(time_converted.has_value());
        ASSERT_TRUE(time2compare == time_converted.value());
    }
}

// TAI->
TEST(CONVERTATION_TAI_TCB, TIME){
    INITIALIZE_LOG(27);
    scalar tol = 1e-15;
    const std::string dir_path = FILE_PATH.substr(0, FILE_PATH.size() - 27);
    std::string filename = dir_path + "/test_data/earth_rotation.csv";
    Ballistics::Time::TimeConverter converter(filename);
    for(int i = 0; i < 500; ++i) {
        Ballistics::Time::Time<Ballistics::Time::Scale::TAI> time(timeResult[i][5], timeResult[i][6]);

        Expected<Ballistics::Time::Time<Ballistics::Time::Scale::TCB>, Error> time_converted =
                converter.convert<Ballistics::Time::Scale::TCB, Ballistics::Time::Scale::TAI>(time);
        Ballistics::Time::Time<Ballistics::Time::Scale::TCB> time2compare(timeResult[i][11], timeResult[i][12]);
        ASSERT_TRUE(time_converted.has_value());
        ASSERT_TRUE(time2compare == time_converted.value());
    }
}

TEST(CONVERTATION_TAI_TT, TIME){
    INITIALIZE_LOG(27);
    scalar tol = 1e-10;
    const std::string dir_path = FILE_PATH.substr(0, FILE_PATH.size() - 27);
    std::string filename = dir_path + "/test_data/earth_rotation.csv";
    Ballistics::Time::TimeConverter converter(filename);
    for(int i = 0; i < 500; ++i) {
        Ballistics::Time::Time<Ballistics::Time::Scale::TAI> time(timeResult[i][5], timeResult[i][6]);
        Expected<Ballistics::Time::Time<Ballistics::Time::Scale::TT>, Error> time_converted = converter.convert<Ballistics::Time::Scale::TT, Ballistics::Time::Scale::TAI>(
                time);
        Ballistics::Time::Time<Ballistics::Time::Scale::TT> time2compare(timeResult[i][7], timeResult[i][8]);
        ASSERT_TRUE(time_converted.has_value());
        ASSERT_TRUE(time2compare == time_converted.value());
    }
}

TEST(CONVERTATION_TAI_UTC, TIME){
    INITIALIZE_LOG(27);
    scalar tol = 1e-10;
    const std::string dir_path = FILE_PATH.substr(0, FILE_PATH.size() - 27);
    std::string filename = dir_path + "/test_data/earth_rotation.csv";
    Ballistics::Time::TimeConverter converter(filename);
    for(int i = 0; i < 500; ++i) {
        Ballistics::Time::Time<Ballistics::Time::Scale::TAI> time(timeResult[i][5], timeResult[i][6]);
        Expected<Ballistics::Time::Time<Ballistics::Time::Scale::UTC>, Error> time_converted = converter.convert<Ballistics::Time::Scale::UTC, Ballistics::Time::Scale::TAI>(
                time);
        Ballistics::Time::Time<Ballistics::Time::Scale::UTC> time2compare(timeResult[i][3], timeResult[i][4]);
        ASSERT_TRUE(time_converted.has_value());
        ASSERT_TRUE(time2compare == time_converted.value());
    }
}

TEST(CONVERTATION_TAI_TCG, TIME){
    INITIALIZE_LOG(27);
    scalar tol = 1e-10;
    const std::string dir_path = FILE_PATH.substr(0, FILE_PATH.size() - 27);
    std::string filename = dir_path + "/test_data/earth_rotation.csv";
    Ballistics::Time::TimeConverter converter(filename);
    for(int i = 0; i < 500; ++i) {
        Ballistics::Time::Time<Ballistics::Time::Scale::TAI> time(timeResult[i][5], timeResult[i][6]);
        Expected<Ballistics::Time::Time<Ballistics::Time::Scale::TCG>, Error> time_converted = converter.convert<Ballistics::Time::Scale::TCG, Ballistics::Time::Scale::TAI>(
                time);
        Ballistics::Time::Time<Ballistics::Time::Scale::TCG> time2compare(timeResult[i][9], timeResult[i][10]);
        ASSERT_TRUE(time_converted.has_value());
        ASSERT_TRUE(time2compare == time_converted.value());
    }
}

// TCB->
TEST(CONVERTATION_TCB_TAI, TIME){
    INITIALIZE_LOG(27);
    scalar tol = 1e-15;
    const std::string dir_path = FILE_PATH.substr(0, FILE_PATH.size() - 27);
    std::string filename = dir_path + "/test_data/earth_rotation.csv";
    Ballistics::Time::TimeConverter converter(filename);
    for(int i = 0; i < 500; ++i) {
        Ballistics::Time::Time<Ballistics::Time::Scale::TCB> time(timeResult[i][11], timeResult[i][12]);

        Expected<Ballistics::Time::Time<Ballistics::Time::Scale::TAI>, Error> time_converted =
                converter.convert<Ballistics::Time::Scale::TAI, Ballistics::Time::Scale::TCB>(time);
        Ballistics::Time::Time<Ballistics::Time::Scale::TAI> time2compare(timeResult[i][5], timeResult[i][6]);
        ASSERT_TRUE(time_converted.has_value());
        ASSERT_TRUE(time2compare == time_converted.value());
    }
}

TEST(CONVERTATION_TCB_TDB, TIME){
    INITIALIZE_LOG(27);
    scalar tol = 1e-10;
    const std::string dir_path = FILE_PATH.substr(0, FILE_PATH.size() - 27);
    std::string filename = dir_path + "/test_data/earth_rotation.csv";
    Ballistics::Time::TimeConverter converter(filename);
    for(int i = 0; i < 500; ++i) {
        Ballistics::Time::Time<Ballistics::Time::Scale::TCB> time(timeResult[i][11], timeResult[i][12]);
        Expected<Ballistics::Time::Time<Ballistics::Time::Scale::TDB>, Error> time_converted = converter.convert<Ballistics::Time::Scale::TDB, Ballistics::Time::Scale::TCB>(
                time);
        Ballistics::Time::Time<Ballistics::Time::Scale::TDB> time2compare(timeResult[i][13], timeResult[i][14]);
        ASSERT_TRUE(time_converted.has_value());
        ASSERT_TRUE(time2compare == time_converted.value());
    }
}

TEST(CONVERTATION_TCB_TT, TIME){
    INITIALIZE_LOG(27);
    scalar tol = 1e-9;
    const std::string dir_path = FILE_PATH.substr(0, FILE_PATH.size() - 27);
    std::string filename = dir_path + "/test_data/earth_rotation.csv";
    Ballistics::Time::TimeConverter converter(filename);
    for(int i = 0; i < 500; ++i) {
        Ballistics::Time::Time<Ballistics::Time::Scale::TCB> time(timeResult[i][11], timeResult[i][12]);
        Expected<Ballistics::Time::Time<Ballistics::Time::Scale::TT>, Error> time_converted = converter.convert<Ballistics::Time::Scale::TT, Ballistics::Time::Scale::TCB>(
                time);
        Ballistics::Time::Time<Ballistics::Time::Scale::TT> time2compare(timeResult[i][7], timeResult[i][8]);
        ASSERT_TRUE(time_converted.has_value());
        ASSERT_TRUE(time2compare == time_converted.value());
    }
}

// TDB->
TEST(CONVERTATION_TDB_TCB, TIME){
    INITIALIZE_LOG(27);
    scalar tol = 1e-15;
    const std::string dir_path = FILE_PATH.substr(0, FILE_PATH.size() - 27);
    std::string filename = dir_path + "/test_data/earth_rotation.csv";
    Ballistics::Time::TimeConverter converter(filename);
    for(int i = 0; i < 500; ++i) {
        Ballistics::Time::Time<Ballistics::Time::Scale::TDB> time(timeResult[i][13], timeResult[i][14]);

        Expected<Ballistics::Time::Time<Ballistics::Time::Scale::TCB>, Error> time_converted =
                converter.convert<Ballistics::Time::Scale::TCB, Ballistics::Time::Scale::TDB>(time);
        Ballistics::Time::Time<Ballistics::Time::Scale::TCB> time2compare(timeResult[i][11], timeResult[i][12]);
        ASSERT_TRUE(time_converted.has_value());
        ASSERT_TRUE(time2compare == time_converted.value());
    }
}

TEST(CONVERTATION_TDB_TT, TIME){
    INITIALIZE_LOG(27)
    scalar tol = 1e-15;
    const std::string dir_path = FILE_PATH.substr(0, FILE_PATH.size() - 27);
    std::string filename = dir_path + "/test_data/earth_rotation.csv";
    Ballistics::Time::TimeConverter converter(filename);
    for(int i = 0; i < 500; ++i) {
        Ballistics::Time::Time<Ballistics::Time::Scale::TDB> time(timeResult[i][13], timeResult[i][14]);
        Expected<Ballistics::Time::Time<Ballistics::Time::Scale::TT>, Error> time_converted = converter.convert<Ballistics::Time::Scale::TT, Ballistics::Time::Scale::TDB>(
                time);
        Ballistics::Time::Time<Ballistics::Time::Scale::TT> time2compare(timeResult[i][7], timeResult[i][8]);
        ASSERT_TRUE(time_converted.has_value());
        ASSERT_TRUE(time2compare == time_converted.value());
    }
}

// TT->
TEST(CONVERTATION_TT_TDB, TIME){
    INITIALIZE_LOG(27);
    scalar tol = 1e-15;
    const std::string dir_path = FILE_PATH.substr(0, FILE_PATH.size() - 27);
    std::string filename = dir_path + "/test_data/earth_rotation.csv";
    Ballistics::Time::TimeConverter converter(filename);
    for(int i = 0; i < 500; ++i) {
        Ballistics::Time::Time<Ballistics::Time::Scale::TT> time(timeResult[i][7], timeResult[i][8]);

        Expected<Ballistics::Time::Time<Ballistics::Time::Scale::TDB>, Error> time_converted =
                converter.convert<Ballistics::Time::Scale::TDB, Ballistics::Time::Scale::TT>(time);
        Ballistics::Time::Time<Ballistics::Time::Scale::TDB> time2compare(timeResult[i][13], timeResult[i][14]);
        ASSERT_TRUE(time_converted.has_value());
        ASSERT_TRUE(time2compare == time_converted.value());
    }
}

TEST(CONVERTATION_TT_TCG, TIME){
    INITIALIZE_LOG(27);
    scalar tol = 1e-10;
    const std::string dir_path = FILE_PATH.substr(0, FILE_PATH.size() - 27);
    std::string filename = dir_path + "/test_data/earth_rotation.csv";
    Ballistics::Time::TimeConverter converter(filename);
    for(int i = 0; i < 500; ++i) {
        Ballistics::Time::Time<Ballistics::Time::Scale::TT> time(timeResult[i][7], timeResult[i][8]);

        Expected<Ballistics::Time::Time<Ballistics::Time::Scale::TCG>, Error> time_converted =
                converter.convert<Ballistics::Time::Scale::TCG, Ballistics::Time::Scale::TT>(time);
        Ballistics::Time::Time<Ballistics::Time::Scale::TCG> time2compare(timeResult[i][9], timeResult[i][10]);
        ASSERT_TRUE(time_converted.has_value());
        ASSERT_TRUE(time2compare == time_converted.value());
    }
}

TEST(CONVERTATION_TT_UT1, TIME){
    INITIALIZE_LOG(27);
    scalar tol = 1e-15;
    const std::string dir_path = FILE_PATH.substr(0, FILE_PATH.size() - 27);
    std::string filename = dir_path + "/test_data/earth_rotation.csv";
    Ballistics::Time::TimeConverter converter(filename);
    for(int i = 0; i < 500; ++i) {
        Ballistics::Time::Time<Ballistics::Time::Scale::TT> time(timeResult[i][7], timeResult[i][8]);
        Expected<Ballistics::Time::Time<Ballistics::Time::Scale::UT1>, Error> time_converted = converter.convert<Ballistics::Time::Scale::UT1, Ballistics::Time::Scale::TT>(
                time);
        Ballistics::Time::Time<Ballistics::Time::Scale::UT1> time2compare(timeResult[i][1], timeResult[i][2]);
        ASSERT_TRUE(time_converted.has_value());
        ASSERT_TRUE(time2compare == time_converted.value());
    }
}

TEST(CONVERTATION_TT_TAI, TIME){
    INITIALIZE_LOG(27);
    scalar tol = 1e-15;
    const std::string dir_path = FILE_PATH.substr(0, FILE_PATH.size() - 27);
    std::string filename = dir_path + "/test_data/earth_rotation.csv";
    Ballistics::Time::TimeConverter converter(filename);
    for(int i = 0; i < 500; ++i) {
        Ballistics::Time::Time<Ballistics::Time::Scale::TT> time(timeResult[i][7], timeResult[i][8]);
        Expected<Ballistics::Time::Time<Ballistics::Time::Scale::TAI>, Error> time_converted = converter.convert<Ballistics::Time::Scale::TAI, Ballistics::Time::Scale::TT>(
                time);
        Ballistics::Time::Time<Ballistics::Time::Scale::TAI> time2compare(timeResult[i][5], timeResult[i][6]);
        ASSERT_TRUE(time_converted.has_value());
        ASSERT_TRUE(time2compare == time_converted.value());
    }
}

TEST(CONVERTATION_TT_UTC, TIME){
    INITIALIZE_LOG(27);
    scalar tol = 1e-15;
    const std::string dir_path = FILE_PATH.substr(0, FILE_PATH.size() - 27);
    std::string filename = dir_path + "/test_data/earth_rotation.csv";
    Ballistics::Time::TimeConverter converter(filename);
    for(int i = 0; i < 500; ++i) {
        Ballistics::Time::Time<Ballistics::Time::Scale::TT> time(timeResult[i][7], timeResult[i][8]);
        Expected<Ballistics::Time::Time<Ballistics::Time::Scale::UTC>, Error> time_converted =
                converter.convert<Ballistics::Time::Scale::UTC, Ballistics::Time::Scale::TT>(time);
        Ballistics::Time::Time<Ballistics::Time::Scale::UTC> time2compare(timeResult[i][3], timeResult[i][4]);
        ASSERT_TRUE(time_converted.has_value());
        ASSERT_TRUE(time2compare == time_converted.value());
    }
}


TEST(CONVERT_JD_MJD, TIME){
    scalar tol = 1e-15;
    Ballistics::Time::Time<Ballistics::Time::Scale::TAI> time(timeResult[10][5], timeResult[10][6]);
    ASSERT_NEAR(timeResult[10][5] + timeResult[10][6] - 2400000.5, time.from_jd_to_mjd(), tol);
}

