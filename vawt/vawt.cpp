#include "vawt.hpp"
#include "Interpolators/_1D/LinearInterpolator.hpp"
#include "streamtube.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <tuple>
#include <vector>

using namespace vawt;
const double PI = boost::math::double_constants::pi;

VAWTSolution VAWTSolver::solve(double beta) {
    return this->solve([beta](double theta) { return beta; });
}

VAWTSolution VAWTSolver::solve(std::function<double(double)> beta) {
    return this->map_streamtubes([beta, this](VAWTCase case_, double theta_up,
                                              double theta_down) {
        double beta_up = beta(theta_up);
        double beta_down = beta(beta_down);
        double a_up =
            StreamTube(theta_up, beta_up, 0.0).solve_a(case_, this->_epsilon);
        double a_down = StreamTube(theta_down, beta_down, a_up)
                            .solve_a(case_, this->_epsilon);
        return std::tuple(beta_up, beta_down, a_up, a_down);
    });
}

VAWTSolution VAWTSolver::map_streamtubes(
    std::function<std::tuple<double, double, double, double>(VAWTCase, double,
                                                             double)>
        solve_fn) {
    auto d_t_half = PI / (double)this->_n_streamtubes;
    std::vector<double> theta(this->_n_streamtubes);
    std::vector<double> beta(this->_n_streamtubes, 0.0);
    std::vector<double> a(this->_n_streamtubes, 0.0);
    std::vector<double> a_0(this->_n_streamtubes, 0.0);

    std::generate(theta.begin(), theta.end(),
                  [i = d_t_half, d_t = 2.0 * d_t_half]() mutable {
                      double current = i;
                      i += d_t;
                      return current;
                  });

    auto case_ = this->get_case();

    for (uint i = 0; i < this->_n_streamtubes / 2; i++) {
        uint i_down = this->_n_streamtubes - 1 - i;

        double theta_up = theta[i];
        double theta_down = theta[i_down];
        auto [beta_up, beta_down, a_up, a_down] =
            solve_fn(case_, theta_up, theta_down);

        beta[i] = beta_up;
        beta[i_down] = beta_down;
        a[i] = a_up;
        a[i_down] = a_down;
        a_0[i_down] = a_up;
    }
    return VAWTSolution(case_, this->_n_streamtubes, theta, beta, a, a_0,
                        this->_epsilon);
}

VAWTCase VAWTSolver::get_case() {
    return VAWTCase{this->_re, this->_tsr, this->_solidity, this->aerofoil};
}

StreamTubeSolution VAWTSolution::solution(double theta) {
    auto a_0 = this->a_0(theta);
    auto a = this->a(theta);
    auto beta = this->beta(theta);
    auto tube = StreamTube(theta, beta, a_0);
    return StreamTubeSolution(this->case_, tube, a);
}

double VAWTSolution::c_torque() {
    double ct = 0.0;
    for (int i = 0; i < this->_theta.size(); i++) {
        auto theta = this->_theta[i];
        auto beta = this->_beta[i];
        auto a = this->_a[i];
        auto a_0 = this->_a_0[i];
        auto tube = StreamTube(theta, beta, a_0);
        auto solution = StreamTubeSolution(this->case_, tube, a);
        ct += solution.c_tan() * pow(solution.w(), 2);
    }
    return ct * this->case_.solidity / (double)this->n_streamtubes;
}
double VAWTSolution::beta(double theta) {
    return _1D::LinearInterpolator<double>(this->_theta, this->_beta)(theta);
}
double VAWTSolution::a(double theta) {
    return _1D::LinearInterpolator<double>(this->_theta, this->_a)(theta);
}
double VAWTSolution::a_0(double theta) {
    return _1D::LinearInterpolator<double>(this->_theta, this->_a_0)(theta);
}
double VAWTSolution::thrust_error(double theta) {
    return this->solution(theta).thrust_error();
}
double VAWTSolution::c_tan(double theta) {
    return this->solution(theta).c_tan();
}
double VAWTSolution::w(double theta) { return this->solution(theta).w(); }
double VAWTSolution::alpha(double theta) {
    return this->solution(theta).alpha();
}
double VAWTSolution::re(double theta) { return this->solution(theta).re(); }
