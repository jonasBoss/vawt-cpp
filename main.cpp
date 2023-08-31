#include <iostream>
#include <ostream>
#include <aerofoil.hpp>

int main(int argc, char **argv) {
    std::cout << "Hello, world!" << std::endl;

    vawt::AerofoilBuilder* builder = new vawt::AerofoilBuilder;

    builder->load_data("examples/NACA0018/NACA0018Re0040.data", 20000.0);

    std::cout << "bye" << std::endl;
    
    return 0;
}
