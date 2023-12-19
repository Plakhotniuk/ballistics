cmake_minimum_required(VERSION 3.16)

project(Ballistics)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(src)

add_subdirectory(third_party)

message("TESTING: ${WITH_TESTS}")

if (WITH_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif ()
