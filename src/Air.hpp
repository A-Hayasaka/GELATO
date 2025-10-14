//
// The MIT License
//
// Copyright (c) 2024 Interstellar Technologies Inc.
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files
// (the "Software"), to deal in the Software without restriction,
// including without limitation the rights to use, copy, modify, merge,
// publish, distribute, sublicense, and/or sell copies of the Software,
// and to permit persons to whom the Software is furnished to do so,
// subject to the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
// CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
// TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
// SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//

/**
 * @file Air.hpp
 * @brief U.S. Standard Atmosphere 1976 model implementation
 * 
 * This header provides atmospheric property calculations according to the
 * U.S. Standard Atmosphere 1976 model. Includes functions for calculating
 * temperature, pressure, density, and speed of sound at various altitudes
 * up to 120 km.
 */

#include <Eigen/Core>
#include <cmath>
#include <vector>

#ifndef SRC_AIR_HPP_
#define SRC_AIR_HPP_

/**
 * @brief Structure containing atmospheric layer parameters
 * 
 * Parameters for a specific atmospheric layer in the U.S. Standard Atmosphere 1976 model.
 */
struct AirParams {
  double Hb;   ///< Base geopotential altitude [m]
  double Lmb;  ///< Temperature gradient [K/m]
  double Tmb;  ///< Base temperature [K]
  double Pb;   ///< Base pressure [Pa]
  double R;    ///< Specific gas constant [J/(kg·K)]
};

/**
 * @brief U.S. Standard Atmosphere 1976 implementation
 * 
 * Provides methods to calculate atmospheric properties (temperature, pressure,
 * density, speed of sound) as a function of altitude according to the
 * U.S. Standard Atmosphere 1976 model.
 */
class Air {
  static const double Rstar, g0, r0;  ///< Physical constants
  static const std::vector<double> hb, lmb, tmb, pb, mb;  ///< Layer parameters

 public:
  /**
   * @brief Convert geometric altitude to geopotential altitude
   * @param geometric_altitude Geometric altitude above sea level [m]
   * @return Geopotential altitude [m]
   */
  static double geopotential_altitude(double geometric_altitude);
  
  /**
   * @brief Get atmospheric layer parameters for given altitude
   * @param altitude Geopotential altitude [m]
   * @return AirParams structure containing layer parameters
   */
  static AirParams us76_params(double altitude);
  
  /**
   * @brief Calculate atmospheric temperature at given altitude
   * @param altitude Geopotential altitude [m]
   * @return Temperature [K]
   */
  static double temperature(double altitude);
  
  /**
   * @brief Calculate atmospheric pressure at given altitude
   * @param altitude Geopotential altitude [m]
   * @return Pressure [Pa]
   */
  static double pressure(double altitude);
  
  /**
   * @brief Calculate atmospheric density at given altitude
   * @param altitude Geopotential altitude [m]
   * @return Density [kg/m³]
   */
  static double density(double altitude);
  
  /**
   * @brief Calculate speed of sound at given altitude
   * @param altitude Geopotential altitude [m]
   * @return Speed of sound [m/s]
   */
  static double speed_of_sound(double altitude);
};

#endif  // SRC_AIR_HPP_
