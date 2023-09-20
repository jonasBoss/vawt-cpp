cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 20)

project(vawt)

include(CTest)

add_executable(vawt-cpp main.cpp)
add_executable(vawt-test tests.cpp)
add_test(NAME vawt-test COMMAND vawt-test)

find_package(Boost REQUIRED)

add_subdirectory(external/libinterpolate)
add_subdirectory(external/csv-parser)
add_subdirectory(vawt)

target_include_directories(vawt-cpp PUBLIC vawt)
target_link_libraries(vawt-cpp PUBLIC vawt)

target_include_directories(vawt-test PUBLIC vawt)
target_link_libraries(vawt-test PUBLIC vawt csv Boost::boost)

#install(TARGETS vawt RUNTIME DESTINATION bin)
