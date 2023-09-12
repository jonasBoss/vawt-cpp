#pragma once

#include <cmath>
#include <math.h>
#include <utility>

/**
 * @brief apply rotation around alpha to x, y
 * 
 * @param x 
 * @param y 
 * @param alpha 
 * @return std::pair<double, double> 
 */
inline std::pair<double, double> rot_vec(double x, double y, double alpha) {
    return std::pair<double, double>(cos(alpha) * x + sin(-alpha) * y,
                                     sin(alpha) * x + cos(alpha) * y);
}