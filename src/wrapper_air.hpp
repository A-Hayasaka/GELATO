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
 * @file wrapper_air.hpp
 * @brief Python binding wrapper functions for atmospheric calculations
 * 
 * This header provides wrapper functions for the Air class to facilitate
 * Python bindings via PyBind11. These functions serve as a simplified
 * interface for atmospheric property calculations.
 */

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <Eigen/Core>
#include <cmath>

#include "Air.hpp"

namespace py = pybind11;

#ifndef SRC_WRAPPER_AIR_HPP_
#define SRC_WRAPPER_AIR_HPP_

using vec3d = Eigen::Matrix<double, 3, 1>;   ///< 3D vector type
using vec4d = Eigen::Matrix<double, 4, 1>;   ///< 4D vector type (for quaternions)
using vecXd = Eigen::Matrix<double, -1, 1>;  ///< Dynamic size vector type
using matXd = Eigen::Matrix<double, -1, -1, Eigen::RowMajor>;  ///< Dynamic size matrix type

/**
 * @brief Wrapper for geopotential altitude calculation
 * @param z Geometric altitude [m]
 * @return Geopotential altitude [m]
 */
double geopotential_altitude(double z) { return Air::geopotential_altitude(z); }

/**
 * @brief Wrapper for atmospheric temperature calculation
 * @param altitude_m Geopotential altitude [m]
 * @return Temperature [K]
 */
double airtemperature_at(double altitude_m) {
  return Air::temperature(altitude_m);
}

/**
 * @brief Wrapper for atmospheric pressure calculation
 * @param altitude_m Geopotential altitude [m]
 * @return Pressure [Pa]
 */
double airpressure_at(double altitude_m) { return Air::pressure(altitude_m); }

/**
 * @brief Wrapper for atmospheric density calculation
 * @param altitude_m Geopotential altitude [m]
 * @return Density [kg/m³]
 */
double airdensity_at(double altitude_m) { return Air::density(altitude_m); }

/**
 * @brief Wrapper for speed of sound calculation
 * @param altitude_m Geopotential altitude [m]
 * @return Speed of sound [m/s]
 */
double speed_of_sound(double altitude_m) {
  return Air::speed_of_sound(altitude_m);
}

#endif  // SRC_WRAPPER_AIR_HPP_
