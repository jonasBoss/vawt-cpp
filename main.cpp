#include <algorithm>
#include <boost/math/constants/constants.hpp>
#include <iostream>
#include <ostream>
#include <aerofoil.hpp>
#include <vector>

const double TO_RAD = boost::math::double_constants::pi / 180;


int main(int argc, char **argv) {
    std::cout << "Hello, world!" << std::endl;

    vawt::AerofoilBuilder* builder = new vawt::AerofoilBuilder;

    auto aerofoil = builder->
        load_data("examples/NACA0018/NACA0018Re0080.data", 80000.0)
        .load_data("examples/NACA0018/NACA0018Re0040.data", 40000.0)
        .load_data("examples/NACA0018/NACA0018Re0160.data", 160000.0)
        .build();

    std::vector<double> alpha(45);
    std::vector<double> re(20);

    double step = 90.0 / 44.0 * TO_RAD; 
    std::generate(alpha.begin(), alpha.end(), [i = 0.0, step]() mutable {
        double current = i;
        i +=step;
        return current;
    });

    step = (180000.0 - 30000.0) / 19; 
    std::generate(re.begin(), re.end(), [i = 30000.0, step]() mutable {
        double current = i;
        i +=step;
        return current;
    });

    std::cout << "cl = [";
    for (auto alpha: alpha){
        for(auto re: re) {
            auto cl = aerofoil->cl_cl(alpha, re).cl();
            std::cout << cl << ", ";
        }
    }
    std::cout << "]" << std::endl;

    std::cout << "alpha = [";
    for (auto alpha: alpha){
        for(auto re: re) {
            std::cout << alpha << ", ";
        }
    }
    std::cout << "]" << std::endl;

    std::cout << "re = [";
    for (auto alpha: alpha){
        for(auto re: re) {
            std::cout << re << ", ";
        }
    }
    std::cout << "]" << std::endl;

    
    return 0;
}
