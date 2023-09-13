#include "vawt.hpp"
#include "streamtube.hpp"

using namespace vawt;

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