cmake_minimum_required(VERSION 3.17)
project(Ballistics)
set(CMAKE_CXX_STANDARD 17)

add_subdirectory(sofa)
add_subdirectory(g3log)
add_subdirectory(geographiclib)
add_subdirectory(calceph)

add_library(Eigen INTERFACE)
target_include_directories(Eigen INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

add_library(fast-cpp-csv-parser INTERFACE)
target_include_directories(fast-cpp-csv-parser INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

add_library(expected INTERFACE)
target_include_directories(expected INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
