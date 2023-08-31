cmake_minimum_required(VERSION 3.16)

project(vawt)

add_library(vawt aerofoil.hpp aerofoil.cpp)

find_package(Boost REQUIRED)

target_include_directories(vawt 
    PUBLIC ../external/csv-parser/single_include
    PUBLIC ../external/libinterpolate/src/libInterpolate
)

target_link_libraries(vawt
    csv
    Boost::boost
    Interpolate
)