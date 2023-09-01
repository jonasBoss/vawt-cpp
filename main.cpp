#include <iostream>
#include <ostream>
#include <aerofoil.hpp>

int main(int argc, char **argv) {
    std::cout << "Hello, world!" << std::endl;

    vawt::AerofoilBuilder* builder = new vawt::AerofoilBuilder;

    builder->
        load_data("examples/NACA0018/NACA0018Re0080.data", 80000.0)
        .load_data("examples/NACA0018/NACA0018Re0040.data", 40000.0)
        .load_data("examples/NACA0018/NACA0018Re0160.data", 160000.0)
        .build();

    std::cout << "bye" << std::endl;
    
    return 0;
}
