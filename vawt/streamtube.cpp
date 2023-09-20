#include "streamtube.hpp"
#include "private_stuff.hpp"
#include "vawt.hpp"
#include <boost/math/constants/constants.hpp>
#include <cmath>

using namespace std;

namespace vawt {

const double PI = boost::math::double_constants::pi;

StreamTube::Velocity StreamTube::Velocity::from_tangetial(double x, double y, double theta) {
    auto [a, b] = rot_vec(x, y, theta);
    return Velocity(a, b);
}

std::pair<double, double> StreamTube::Velocity::to_foil(double theta, double beta) {
    return rot_vec(x, y, -theta - beta);
}

double StreamTube::Velocity::magnitude() {
    return sqrt(pow(x,2)+pow(y,2));
}

tuple<double, double, double> StreamTube::w_alpha_re(double a, VAWTCase case_) {
    auto w = this->w_vec(a, case_);
    auto [w_x_foil, w_y_foil] = w.to_foil(this->theta, this->beta);
    auto alpha = atan2(w_y_foil, w_x_foil) + PI / 2.0;
    auto w_norm = w.magnitude();
    auto re = case_.re * w_norm;
    return tuple(w_norm, alpha, re);
}

double StreamTube::c_tan(double a, VAWTCase case_) {
    auto [w, alpha, re] = this->w_alpha_re(a, case_);
    return get<1>(
        case_.aerofoil->cl_cd(alpha, re).to_tangential(alpha, this->beta));
}

double StreamTube::a_strickland(VAWTCase case_) {
    double a = 0.0;
    for (int i = 0; i < 10; i++) {
        auto c_s = this->foil_thrust(a, case_);
        auto a_new = 0.25 * c_s + pow(a, 2);
        if (a_new < 1.0) {
            a = a_new;
        } else {
            a = 1.0;
        }
    }
    return a;
}

double StreamTube::foil_thrust(double a, VAWTCase case_) {
    auto [w, alpha, re] = this->w_alpha_re(a, case_);

    auto cl_cd = case_.aerofoil->cl_cd(alpha, re);
    auto [_, force_coeff] = cl_cd.to_global(alpha, this->beta, this->theta);
    return -force_coeff * pow(w / this->c_0(), 2) * case_.solidity /
           (PI * abs(sin(this->theta)));
}

double StreamTube::wind_thrust(double a) {
    if (a < 0.4) {
        return 4.0 * a * (1.0 - a);
    } else {
        return 26.0 / 15.0 * a + 4.0 / 15.0;
    }
}

double StreamTube::solve_a(VAWTCase case_, double epsilon) {
    double a_left = -2.0;
    double a_right = 2.0;
    double err_left = 0.0;
    double err_right = 0.0;
    err_left = this->thrust_error(a_left, case_);
    err_right = this->thrust_error(a_right, case_);
    if (err_left * err_right > 0.0) {
        return this->a_strickland(case_);
    }
    while ((a_right - a_left) > epsilon) {
        double a = a_left + (a_right - a_left) / 2.0;
        double err = this->thrust_error(a, case_);

        if (err_left * err <= 0.0) {
            a_right = a;
            err_right = err;
        } else {
            a_left = a;
            err_left = err;
        }
    }
    return a_left + (a_right - a_left) / 2.0;
}
} // namespace vawt