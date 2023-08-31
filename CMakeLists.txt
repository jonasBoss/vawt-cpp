cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 20)

project(vawt)

add_executable(vawt-cpp main.cpp)

add_subdirectory(external/libinterpolate)
add_subdirectory(external/csv-parser)
add_subdirectory(vawt)

target_include_directories(vawt-cpp PUBLIC vawt)

target_link_libraries(vawt-cpp vawt)

#install(TARGETS vawt RUNTIME DESTINATION bin)
