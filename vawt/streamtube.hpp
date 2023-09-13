#pragma once

#include "vawt.hpp"
#include <tuple>
#include <utility>
namespace vawt {

class StreamTube {
  private:
    double a_0;
    double theta;
    double beta;

    /**
     * @brief the difference between the wind thrust and the foil force for a
     * given induction factor a.
     * for good solutions this should be small
     * @param a - induction factor
     * @param case_ - case settings
     * @return double
     */
    double thrust_error(double a, VAWTCase case_) {
        return this->foil_thrust(a, case_) - StreamTube::wind_thrust(a);
    }

    /**
     * @brief the relative velocity magnitude `w`, the angle of attack `alpha`
     * in radians and the local reynolds number `re` at the foil for a given
     * induction factor a
     *
     * @param a
     * @param case_
     * @return std::tuple<double, double, double>
     */
    std::tuple<double, double, double> w_alpha_re(double a, VAWTCase case_);

    /**
     * @brief tangential foil coefficient
     *
     * @param a
     * @param case_
     * @return double
     */
    double c_tan(double a, VAWTCase case_);
    double a_strickland(VAWTCase case_);
    double foil_thrust(double a, VAWTCase case_);

    /**
     * @brief Thrust coefficient by momentum theory or Glauert empirical formula

     * A crude straight line approximation for Glauert formula is used
     * between 0.4 < a < 1.0,  0.96 < CtubeThru < 2.0
     *
     * @param a
     * @return double
     */
    static double wind_thrust(double a);

    /**
     * @brief reference windspeed
     *
     * @return double
     */
    double c_0() { return 1.0 - 2.0 * this->a_0; }

    class Velocity {
      private:
        double x, y;
        Velocity(double x, double y) : x(x), y(y) {}

      public:
        static Velocity from_global(double x, double y) {
            return Velocity(x, y);
        }
        static Velocity from_tangetial(double x, double y, double theta);
        Velocity operator-(Velocity rhs) {
            return Velocity(this->x - rhs.x, this->y - rhs.y);
        }
        std::pair<double, double> to_foil(double theta, double beta);
        double magnitude();
    };

    /**
     * @brief windspeed at the foil
     *
     * @param a
     * @return Velocity
     */
    Velocity c_1_vec(double a) {
        return Velocity::from_global(0.0, -this->c_0() * (1.0 - a));
    }

    /**
     * @brief relative velocity at foil in global xy coordinates
     *
     * @param a
     * @param case_
     * @return Velocity
     */
    Velocity w_vec(double a, VAWTCase case_) {
        return this->c_1_vec(a) -
               Velocity::from_tangetial(0.0, case_.tsr, this->theta);
    }

  public:
    /**
     * @brief Construct a new StreamTube object
     *
     * @param theta - stramtube position in turbine (radians)
     * @param beta - foil pitch angle (radians)
     * @param a_0 - upstream induction factor when `theta < PI` this is probably
     * 0
     */
    StreamTube(double theta, double beta, double a_0) {
        this->a_0 = a_0;
        this->beta = beta;
        this->theta = theta;
    }

    /**
     * @brief solve the streamtube for induction factor a
     *
     * @param case_
     * @param epsilon
     * @return double
     */
    double solve_a(VAWTCase case_, double epsilon);
};

} // namespace vawt