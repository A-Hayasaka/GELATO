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
 * @file Earth.hpp
 * @brief WGS84 Earth model and geodetic coordinate utilities
 * 
 * This header provides WGS84 Earth model constants and functions for
 * converting between ECEF and geodetic coordinates. Also includes
 * Vincenty formula for calculating geodesic distances on the Earth's
 * ellipsoidal surface.
 */

#include <Eigen/Core>
#include <cmath>
#include <utility>

#ifndef SRC_EARTH_HPP_
#define SRC_EARTH_HPP_

/**
 * @brief Earth model parameters and coordinate conversion utilities
 * 
 * Provides WGS84 Earth model constants and methods for converting between
 * ECEF and geodetic coordinates, as well as calculating distances using
 * the Vincenty formula.
 */
class Earth {
 public:
  static const double mu;                ///< Earth's gravitational parameter [m³/s²]
  static const double omega_earth_rps;   ///< Earth's rotation rate [rad/s]
  static const double Ra;                ///< WGS84 equatorial radius [m]
  static const double f;                 ///< WGS84 flattening [-]
  static const double Rb;                ///< WGS84 polar radius [m]
  static const double e2;                ///< First eccentricity squared [-]
  static const double ep2;               ///< Second eccentricity squared [-]

  /**
   * @brief Helper function to compute square of a value
   * @param x Input value
   * @return x²
   */
  static inline double pow2(double x) { return x * x; }

  /**
   * @brief Convert ECEF coordinates to geodetic coordinates
   * @param pos_ecef Position in ECEF frame [m] as (x, y, z)
   * @return Geodetic coordinates (latitude [rad], longitude [rad], altitude [m])
   */
  static Eigen::Vector3d ecef2geodetic(Eigen::Vector3d pos_ecef);
  
  /**
   * @brief Convert geodetic coordinates to ECEF coordinates
   * @param geodetic Geodetic coordinates (latitude [rad], longitude [rad], altitude [m])
   * @return Position in ECEF frame [m] as (x, y, z)
   */
  static Eigen::Vector3d geodetic2ecef(Eigen::Vector3d geodetic);
  
  /**
   * @brief Calculate distance and azimuth between two points using Vincenty formula
   * @param observer_LLH Observer position (latitude [rad], longitude [rad], height [m])
   * @param target_LLH Target position (latitude [rad], longitude [rad], height [m])
   * @return Pair of (distance [m], azimuth [rad])
   */
  static std::pair<double, double> distance_vincenty(
      Eigen::Vector3d observer_LLH, Eigen::Vector3d target_LLH);
};

#endif  // SRC_EARTH_HPP_
