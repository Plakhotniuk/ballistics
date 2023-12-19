################################
# GTest
################################
cmake_minimum_required(VERSION 3.17)
project(Ballistics)
ADD_SUBDIRECTORY(googletest-master)
enable_testing()

include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_subdirectory(utility)

file(GLOB_RECURSE TEST_SOURCES test_*.cpp)  # В переменной TEST_SOURCES лежат все тесты, имена файлов которые начинаются с test_

foreach(TEST_FILE IN LISTS TEST_SOURCES)
    get_filename_component(TEST_FILE_NAME ${TEST_FILE} NAME_WE)
    set(EXECUTABLE_NAME run_${TEST_FILE_NAME})
    add_executable(${EXECUTABLE_NAME} ${TEST_FILE})
    target_link_libraries(
            ${EXECUTABLE_NAME}
            ballistics
            gtest
            gtest_main
            sofa
            calceph
            GeographicLib
            g3log
            fast-cpp-csv-parser
    )
    message("TEST_FILE: ${TEST_FILE}")
    add_test(NAME ${EXECUTABLE_NAME} COMMAND ${EXECUTABLE_NAME})
endforeach()
