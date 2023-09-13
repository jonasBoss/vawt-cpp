#pragma once

#include "aerofoil.hpp"
#include <functional>
#include <memory>
#include <sys/types.h>
#include <vector>

namespace vawt {

class VAWTSolution;
class StreamTubeSolution;

class VAWTSolver {
  private:
    std::shared_ptr<Aerofoil> aerofoil;
    uint _n_streamtubes = 50;
    double _tsr = 2.0;
    double _re = 60'000.0;
    double _solidity = 0.1;
    double _epsilon = 0.01;

  public:
    /**
     * @brief create a new Solver with the following default values:
     *
     * - `n_streamtubes = 50` Number of streamtubes over the whole turbine
     * - `tsr = 2.0` Tipspeed ratio of the turbine
     * - `re = 60_000.0` Reynolds number of the turbine
     * - `solidity = 0.1` Solidity of the Turbine
     * - `epsilon = 0.01` the solution accuracy for a
     * - `particles = 8` the number of particles for beta optimization
     * - `iterations = 30` the number of iterations for beta optimization
     *
     * @param aerofoil
     */
    VAWTSolver(std::shared_ptr<Aerofoil> aerofoil) {
        this->aerofoil = aerofoil;
    }

    /**
     * @brief  update the number of streamtubes for the solution if n is not a
     * multiple of 2 `n+1` is used.
     *
     * @param n
     * @return VAWTSolver&
     */
    VAWTSolver& n_streamtubes(uint n) {
        if (n % 2 != 0)
            n++;
        this->_n_streamtubes = n;
        return *this;
    }

    /**
     * @brief update the tipspeed ratio for the solution
     *
     * @param tsr
     * @return VAWTSolver&
     */
    VAWTSolver& tsr(double tsr) {
        this->_tsr = tsr;
        return *this;
    }

    /**
     * @brief update the raynolds number for the solution
     *
     * @param re
     * @return VAWTSolver&
     */
    VAWTSolver& re(double re) {
        this->_re = re;
        return *this;
    }

    /**
     * @brief update the turbine solidity for the solution
     *
     * @param solidity
     * @return VAWTSolver&
     */
    VAWTSolver& solidity(double solidity) {
        this->_solidity = solidity;
        return *this;
    }

    VAWTSolver& epsilon(double epsilon) {
        this->_epsilon = epsilon;
        return *this;
    }

    VAWTSolution solve(double beta);
    VAWTSolution solve(std::function<double(double)> beta);
};

/**
 * @brief Turbine settings for the VAWT case
 */
struct VAWTCase {

    /**
     * @brief Raynoldsnumber of the turbine
     */
    double re;

    /**
     * @brief Tipspeed ratio of the turbine
     */
    double tsr;

    /**
     * @brief Turbine solidity
     */
    double solidity;

    /**
     * @brief Aerofoil
     */
    std::shared_ptr<Aerofoil> aerofoil;
};

class VAWTSolution {
  private:
    VAWTCase case_;
    uint n_streamtubes;
    std::vector<double> _theta;
    std::vector<double> _beta;
    std::vector<double> _a;
    std::vector<double> _a_0;
    double _epsilon;
    StreamTubeSolution solution(double theta);

  public:
    /**
     * @brief Torque ceofficient of the turbine
     *
     * @return double
     */
    double c_torque();

    /**
     * @brief Power coefficient of the turbine
     *
     * @return double
     */
    double c_power() { return this->c_torque() * this->case_.tsr; }

    /**
     * @brief the pitch angle `beta` at the location `theta`
     *
     * @param theta
     * @return double
     */
    double beta(double theta);

    /**
     * @brief the induction factor `a` at the location `theta`
     *
     * @param theta
     * @return double
     */
    double a(double theta);

    /**
     * @brief the upstream induction factor `a_0` at the location `theta`
     *
     * @param theta
     * @return double
     */
    double a_0(double theta);

    /**
     * @brief the difference between the wind thrust and the foil force
     * (solution error) at the location `theta`
     *
     * @param theta
     * @return double
     */
    double thrust_error(double theta);

    /**
     * @brief tangential foil coefficient at the location `theta`
     *
     * coefficient of lift and drag evaluated in tangential direction
     *
     * @param theta
     * @return double
     */
    double c_tan(double theta);
    double epsilon() { return this->_epsilon; }

    /**
     * @brief the relative windspeed at the foil at location `theta`
     *
     * @param theta
     * @return double
     */
    double w(double theta);

    /**
     * @brief the angle of attac at the foil at location `theta`
     *
     * @param theta
     * @return double
     */
    double alpha(double theta);

    /**
     * @brief the local reynolds number at the foil at location `theta`
     *
     * @param theta
     * @return double
     */
    double re(double theta);
};

} // namespace vawt
