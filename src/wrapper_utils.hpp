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
 * @file wrapper_utils.hpp
 * @brief Utility functions for Python bindings
 * 
 * This header provides utility functions including distance calculations
 * and interpolation routines for use in Python via PyBind11.
 */

#include <algorithm>
#include <vector>

#include "wrapper_air.hpp"
#include "wrapper_coordinate.hpp"

namespace py = pybind11;

#ifndef SRC_WRAPPER_UTILS_HPP_
#define SRC_WRAPPER_UTILS_HPP_

/**
 * @brief Calculate great circle distance using Haversine formula
 * @param lon1 Longitude of first point [deg]
 * @param lat1 Latitude of first point [deg]
 * @param lon2 Longitude of second point [deg]
 * @param lat2 Latitude of second point [deg]
 * @param r Radius of sphere [m]
 * @return Great circle distance [m]
 */
double haversine(double lon1, double lat1, double lon2, double lat2, double r) {
  lon1 = lon1 * M_PI / 180.0;
  lat1 = lat1 * M_PI / 180.0;
  lon2 = lon2 * M_PI / 180.0;
  lat2 = lat2 * M_PI / 180.0;

  double dlon = lon2 - lon1;
  double dlat = lat2 - lat1;
  double a =
      pow(sin(dlat / 2), 2) + cos(lat1) * cos(lat2) * pow(sin(dlon / 2), 2);

  return 2 * r * asin(sqrt(a));
}

/**
 * @brief Linear interpolation
 * @param x Value to interpolate at
 * @param xp X values of data points (must be sorted in ascending order)
 * @param yp Y values of data points
 * @return Interpolated value at x
 * 
 * Uses binary search for efficiency. Values outside the range
 * of xp are clamped to the boundary values.
 */
double interp(double x, vecXd xp, vecXd yp) {
  int n = xp.rows();

  if (x < xp(0)) {
    return yp(0);
  }

  if (x > xp(n - 1)) {
    return yp(n - 1);
  }

  // Find the right place in the table by means of binary search
  auto lower = std::lower_bound(xp.data(), xp.data() + n, x);
  int idx = lower - xp.data() - 1;

  double x_lower = xp(idx);
  double x_upper = xp(idx + 1);

  double y_lower = yp(idx);
  double y_upper = yp(idx + 1);

  double alpha = (x - x_lower) / (x_upper - x_lower);

  return y_lower + alpha * (y_upper - y_lower);
}

/**
 * @brief Get wind velocity in NED frame at given altitude
 * @param altitude_m Altitude [m]
 * @param wind_data Wind data table (altitude, u-component, v-component)
 * @return Wind velocity in NED frame [m/s] (north, east, down)
 */
vec3d wind_ned(double altitude_m, matXd wind_data) {
  double wind_u = interp(altitude_m, wind_data.col(0), wind_data.col(1));
  double wind_v = interp(altitude_m, wind_data.col(0), wind_data.col(2));

  return vec3d(wind_u, wind_v, 0.0);
}

/**
 * @brief Calculate total angle of attack
 * @param pos_eci Position in ECI frame [m]
 * @param vel_eci Velocity in ECI frame [m/s]
 * @param quat Attitude quaternion (ECI to body) [w, x, y, z]
 * @param t Time [s]
 * @param wind Wind data table
 * @return Total angle of attack [rad]
 */
double angle_of_attack_all_rad(vec3d pos_eci, vec3d vel_eci, vec4d quat,
                               double t, matXd wind) {
  vec3d thrust_dir_eci = quatrot(conj(quat), vec3d(1.0, 0.0, 0.0));

  vec3d pos_llh = ecef2geodetic(pos_eci[0], pos_eci[1], pos_eci[2]);
  double altitude = geopotential_altitude(pos_llh[2]);

  vec3d vel_ecef = vel_eci2ecef(vel_eci, pos_eci, t);
  vec3d vel_wind_ned = wind_ned(altitude, wind);

  vec3d vel_wind_eci = quatrot(quat_nedg2eci(pos_eci, t), vel_wind_ned);
  vec3d vel_air_eci = ecef2eci(vel_ecef, t) - vel_wind_eci;

  double c_alpha = normalize(vel_air_eci).dot(normalize(thrust_dir_eci));

  if (c_alpha > 1.0) {
    return 0.0;
  } else if (vel_air_eci.norm() < 1e-6) {
    return 0.0;
  } else {
    return acos(c_alpha);
  }
}

/**
 * @brief Calculate total angle of attack for array of states
 * @param pos_eci Position array in ECI frame [m]
 * @param vel_eci Velocity array in ECI frame [m/s]
 * @param quat Attitude quaternion array (ECI to body) [w, x, y, z]
 * @param t Time array [s]
 * @param wind Wind data table
 * @return Total angle of attack array [rad]
 */
vecXd angle_of_attack_all_array_rad(matXd pos_eci, matXd vel_eci, matXd quat,
                                    vecXd t, matXd wind) {
  int n = pos_eci.rows();
  vecXd alpha(n);

  for (int i = 0; i < n; i++) {
    alpha(i) = angle_of_attack_all_rad(pos_eci.row(i), vel_eci.row(i),
                                       quat.row(i), t(i), wind);
  }

  return alpha;
}

/**
 * @brief Calculate angle of attack components (alpha and beta)
 * @param pos_eci Position in ECI frame [m]
 * @param vel_eci Velocity in ECI frame [m/s]
 * @param quat Attitude quaternion (ECI to body) [w, x, y, z]
 * @param t Time [s]
 * @param wind Wind data table
 * @return (alpha_z, alpha_y) angle of attack components [rad]
 */
Eigen::Vector2d angle_of_attack_ab_rad(vec3d pos_eci, vec3d vel_eci, vec4d quat,
                                       double t, matXd wind) {
  vec3d thrust_dir_eci = quatrot(conj(quat), vec3d(1.0, 0.0, 0.0));

  vec3d pos_llh = ecef2geodetic(pos_eci[0], pos_eci[1], pos_eci[2]);
  double altitude = geopotential_altitude(pos_llh[2]);

  vec3d vel_ecef = vel_eci2ecef(vel_eci, pos_eci, t);
  vec3d vel_wind_ned = wind_ned(altitude, wind);

  vec3d vel_wind_eci = quatrot(quat_nedg2eci(pos_eci, t), vel_wind_ned);
  vec3d vel_air_eci = ecef2eci(vel_ecef, t) - vel_wind_eci;

  vec3d vel_air_body = quatrot(quat, vel_air_eci);

  if (vel_air_body[0] < 1e-6) {
    return Eigen::Vector2d(0.0, 0.0);
  } else {
    double alpha_z = atan2(vel_air_body[2], vel_air_body[0]);
    double alpha_y = atan2(vel_air_body[1], vel_air_body[0]);
    return Eigen::Vector2d(alpha_z, alpha_y);
  }
}

/**
 * @brief Calculate angle of attack components for array of states
 * @param pos_eci Position array in ECI frame [m]
 * @param vel_eci Velocity array in ECI frame [m/s]
 * @param quat Attitude quaternion array (ECI to body) [w, x, y, z]
 * @param t Time array [s]
 * @param wind Wind data table
 * @return Array of (alpha_z, alpha_y) angle of attack components [rad]
 */
Eigen::MatrixXd angle_of_attack_ab_array_rad(matXd pos_eci, matXd vel_eci,
                                             matXd quat, vecXd t, matXd wind) {
  int n = pos_eci.rows();
  Eigen::MatrixXd alpha(n, 2);

  for (int i = 0; i < n; i++) {
    alpha.row(i) = angle_of_attack_ab_rad(pos_eci.row(i), vel_eci.row(i),
                                          quat.row(i), t(i), wind);
  }

  return alpha;
}

/**
 * @brief Calculate dynamic pressure at a given state
 * @param pos_eci Position in ECI frame [m]
 * @param vel_eci Velocity in ECI frame [m/s]
 * @param t Time [s]
 * @param wind Wind data table
 * @return Dynamic pressure [Pa]
 */
double dynamic_pressure_pa(vec3d pos_eci, vec3d vel_eci, double t, matXd wind) {
  vec3d pos_llh = ecef2geodetic(pos_eci[0], pos_eci[1], pos_eci[2]);
  double altitude = geopotential_altitude(pos_llh[2]);
  double rho = airdensity_at(altitude);

  vec3d vel_ecef = vel_eci2ecef(vel_eci, pos_eci, t);
  vec3d vel_wind_ned = wind_ned(altitude, wind);
  vec3d vel_wind_eci = quatrot(quat_nedg2eci(pos_eci, t), vel_wind_ned);
  vec3d vel_air_eci = ecef2eci(vel_ecef, t) - vel_wind_eci;

  return 0.5 * rho * vel_air_eci.norm() * vel_air_eci.norm();
}

/**
 * @brief Calculate dynamic pressure for array of states
 * @param pos_eci Position array in ECI frame [m]
 * @param vel_eci Velocity array in ECI frame [m/s]
 * @param t Time array [s]
 * @param wind Wind data table
 * @return Dynamic pressure array [Pa]
 */
Eigen::VectorXd dynamic_pressure_array_pa(matXd pos_eci, matXd vel_eci, vecXd t,
                                          matXd wind) {
  int n = pos_eci.rows();
  Eigen::VectorXd q(n);

  for (int i = 0; i < n; i++) {
    q(i) = dynamic_pressure_pa(pos_eci.row(i), vel_eci.row(i), t(i), wind);
  }

  return q;
}

/**
 * @brief Calculate aerodynamic load (q·α) at a given state
 * @param pos_eci Position in ECI frame [m]
 * @param vel_eci Velocity in ECI frame [m/s]
 * @param quat Attitude quaternion (ECI to body) [w, x, y, z]
 * @param t Time [s]
 * @param wind Wind data table
 * @return Aerodynamic load (dynamic pressure × angle of attack) [Pa·rad]
 */
double q_alpha_pa_rad(vec3d pos_eci, vec3d vel_eci, vec4d quat, double t,
                      matXd wind) {
  double alpha = angle_of_attack_all_rad(pos_eci, vel_eci, quat, t, wind);
  double q = dynamic_pressure_pa(pos_eci, vel_eci, t, wind);
  return q * alpha;
}

/**
 * @brief Calculate aerodynamic load (q·α) for array of states
 * @param pos_eci Position array in ECI frame [m]
 * @param vel_eci Velocity array in ECI frame [m/s]
 * @param quat Attitude quaternion array (ECI to body) [w, x, y, z]
 * @param t Time array [s]
 * @param wind Wind data table
 * @return Aerodynamic load array (dynamic pressure × angle of attack) [Pa·rad]
 */
Eigen::VectorXd q_alpha_array_pa_rad(matXd pos_eci, matXd vel_eci, matXd quat,
                                     vecXd t, matXd wind) {
  int n = pos_eci.rows();
  Eigen::VectorXd q_alpha(n);

  for (int i = 0; i < n; i++) {
    q_alpha(i) =
        q_alpha_pa_rad(pos_eci.row(i), vel_eci.row(i), quat.row(i), t(i), wind);
  }

  return q_alpha;
}

#endif  // SRC_WRAPPER_UTILS_HPP_
