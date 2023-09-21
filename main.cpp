#include <boost/math/constants/constants.hpp>
#include <vawt.hpp>
#include <chrono>
#include <iostream>

using namespace vawt;
const double TO_RAD = boost::math::double_constants::pi / 180;

static std::shared_ptr<Aerofoil> load_naca0018() {
    vawt::AerofoilBuilder* builder = new vawt::AerofoilBuilder;
    return builder->load_data("examples/NACA0018/NACA0018Re0080.data", 80'000.0)
        .load_data("examples/NACA0018/NACA0018Re0040.data", 40'000.0)
        .load_data("examples/NACA0018/NACA0018Re0160.data", 160'000.0)
        .set_aspect_ratio(12.8)
        .update_aspect_ratio(true)
        .symmetric(true)
        .build();
}

static VAWTSolver setup_solver(std::shared_ptr<Aerofoil> foil) {
    auto testcase = VAWTSolver(foil);
    testcase.re(31'300.0).solidity(0.3525).n_streamtubes(72).tsr(3.25);
    return testcase;
}

int main(int argc, char** argv) {
    auto foil = load_naca0018();
    auto testcase = setup_solver(foil);

    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i<10'000; i++){
        testcase.solve(0.0);
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Duration for 10'000 solutions: " << duration.count() << " microseconds" << std::endl;

    return 0;
}
