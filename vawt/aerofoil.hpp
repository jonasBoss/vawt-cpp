#ifndef AEROFOIL_HPP
#define AEROFOIL_HPP

#include "Interpolators/_2D/BilinearInterpolator.hpp"
#include <string_view>
#include <tuple>
#include <utility>
#include <list>
#include <vector>
#include <limits>
#include <Interpolate.hpp>

namespace vawt {

/**
 * @brief Aerofoil coefficients of lift and drag
 * 
 */
class ClCd{
private:
    double _cl;
    double _cd;

public:
    /**
     * @brief Coefficient of Lift
     * 
     * @return double 
     */
    double cl() {return this->_cl;};

    /**
     * @brief Coefficient of Drag
     * 
     * @return double 
     */
    double cd() {return this->_cd;};

    /**
     * @brief Convert coefficients to normal and tangential turbine coordinates
     * 
     * @param alpha - the foil angle of attac  
     * @param beta - the pitch agle from turbine tangent to wing chord
     * @return std::pair<double, double> - (normal coefficinent, tangential coefficient)
     */
    std::pair<double, double> to_tangential(double alpha, double beta);

    /**
     * @brief Convert coefficients to global xy direction
     * 
     * @param alpha - the foil angle of attac  
     * @param beta - the pitch agle from turbine tangent to wing chord
     * @param theta - the position angle at the turbine
     * @return std::pair<double, double> - (x-coefficient, y-coefficient)
     */
    std::pair<double, double> to_global(double alpha, double beta, double theta);
};

class Aerofoil{
private:
    bool symmetric;
    _2D::BilinearInterpolator<double> cl;
    _2D::BilinearInterpolator<double> cd;

public:
    ClCd cl_cl(double alpha, double re);
};

class AerofoilBuilder {
private:
    std::vector<double> alpha, re, cl, cd;
    bool symmetric = false;
    double aspect_ratio = std::numeric_limits<double>::infinity();

public:
    AerofoilBuilder& add_data(std::string_view file, double re);

};

}

#endif // AEROFOIL_HPP