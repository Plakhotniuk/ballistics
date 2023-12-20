set(MODULE_NAME ballistics)

message("ballistics_build" ${MODULE_NAME})

file(GLOB_RECURSE SOURCES *.h *.cpp *.hpp)

add_library(${MODULE_NAME} STATIC ${SOURCES})

set_target_properties(${MODULE_NAME} PROPERTIES LINKER_LANGUAGE CXX)

target_include_directories(${MODULE_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_include_directories(${MODULE_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)


target_link_libraries(${MODULE_NAME} PUBLIC GeographicLib sofa Eigen g3log fast-cpp-csv-parser calceph)
