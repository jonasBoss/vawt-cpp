cmake_minimum_required(VERSION 3.0)

project(vawt-cpp)

add_executable(vawt-cpp main.cpp)

install(TARGETS vawt-cpp RUNTIME DESTINATION bin)
