set(MODULE_NAME Utility)
file(GLOB_RECURSE SOURCES *.h *.cpp *.hpp)

add_library(${MODULE_NAME} STATIC ${SOURCES})
target_link_libraries(${MODULE_NAME} PUBLIC)
set_target_properties(${MODULE_NAME} PROPERTIES LINKER_LANGUAGE CXX)

target_include_directories(${MODULE_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
