cmake_minimum_required(VERSION 3.5)

project(vawt)

add_executable(vawt-cpp main.cpp)

add_subdirectory(external/libinterpolate)
add_subdirectory(external/csv-parser)
add_subdirectory(vawt)

target_link_libraries(vawt-cpp csv)

#install(TARGETS vawt RUNTIME DESTINATION bin)
