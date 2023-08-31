cmake_minimum_required(VERSION 3.5)

project(vawt)

add_executable(vawt-cpp main.cpp)

add_subdirectory(vawt)

#install(TARGETS vawt RUNTIME DESTINATION bin)
