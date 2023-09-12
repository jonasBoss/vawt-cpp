#pragma once

#include "aerofoil.hpp"
#include <memory>
#include <sys/types.h>

namespace vawt {
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

} // namespace vawt
