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
 * @file Coordinate.hpp
 * @brief Coordinate transformation utilities for aerospace applications
 * 
 * This header provides comprehensive coordinate transformation functions between
 * various reference frames used in aerospace applications including:
 * - ECEF (Earth-Centered Earth-Fixed)
 * - ECI (Earth-Centered Inertial)
 * - NED (North-East-Down)
 * - Body-fixed frame
 * 
 * Also includes quaternion and Euler angle conversions, and orbital element calculations.
 */

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "Earth.hpp"

#ifndef SRC_COORDINATE_HPP_
#define SRC_COORDINATE_HPP_

/**
 * @brief Coordinate transformation utility class
 * 
 * Provides static methods for transforming coordinates between different
 * reference frames (ECEF, ECI, NED, Body) and handling quaternion operations.
 * Uses Eigen library for linear algebra operations.
 */
class Coordinate {
 public:
  /**
   * @brief Convert ECEF position to ECI position
   * @param xyz_ecef Position vector in ECEF frame [m]
   * @param t Time [s]
   * @return Position vector in ECI frame [m]
   */
  static Eigen::Vector3d ecef2eci(Eigen::Vector3d xyz_ecef, double t);
  
  /**
   * @brief Convert ECI position to ECEF position
   * @param xyz_eci Position vector in ECI frame [m]
   * @param t Time [s]
   * @return Position vector in ECEF frame [m]
   */
  static Eigen::Vector3d eci2ecef(Eigen::Vector3d xyz_eci, double t);
  
  /**
   * @brief Convert ECEF ground velocity to ECI inertialvelocity
   * @param vel_ecef Velocity vector in ECEF frame [m/s]
   * @param pos_ecef Position vector in ECEF frame [m]
   * @param t Time [s]
   * @return Velocity vector in ECI frame [m/s]
   */
  static Eigen::Vector3d vel_ecef2eci(Eigen::Vector3d vel_ecef,
                                      Eigen::Vector3d pos_ecef, double t);
  
  /**
   * @brief Convert ECI inertial velocity to ECEF ground velocity
   * @param vel_eci Velocity vector in ECI frame [m/s]
   * @param pos_eci Position vector in ECI frame [m]
   * @param t Time [s]
   * @return Velocity vector in ECEF frame [m/s]
   */
  static Eigen::Vector3d vel_eci2ecef(Eigen::Vector3d vel_eci,
                                      Eigen::Vector3d pos_eci, double t);
  
  /**
   * @brief Get quaternion for ECI to ECEF transformation
   * @param t Time [s]
   * @return Quaternion representing rotation from ECI to ECEF
   */
  static Eigen::Quaterniond quat_eci2ecef(double t);
  
  /**
   * @brief Get quaternion for ECEF to ECI transformation
   * @param t Time [s]
   * @return Quaternion representing rotation from ECEF to ECI
   */
  static Eigen::Quaterniond quat_ecef2eci(double t);
  
  /**
   * @brief Get quaternion for ECEF to NED transformation
   * @param pos_ecef Position vector in ECEF frame [m]
   * @return Quaternion representing rotation from ECEF to NED
   */
  static Eigen::Quaterniond quat_ecef2ned(Eigen::Vector3d pos_ecef);
  
  /**
   * @brief Get quaternion for NED to ECEF transformation
   * @param pos_ecef Position vector in ECEF frame [m]
   * @return Quaternion representing rotation from NED to ECEF
   */
  static Eigen::Quaterniond quat_ned2ecef(Eigen::Vector3d pos_ecef);
  
  /**
   * @brief Get quaternion for ECI to NED transformation
   * @param pos_eci Position vector in ECI frame [m]
   * @param t Time [s]
   * @return Quaternion representing rotation from ECI to NED
   */
  static Eigen::Quaterniond quat_eci2ned(Eigen::Vector3d pos_eci, double t);
  
  /**
   * @brief Get quaternion for NED to ECI transformation
   * @param pos_eci Position vector in ECI frame [m]
   * @param t Time [s]
   * @return Quaternion representing rotation from NED to ECI
   */
  static Eigen::Quaterniond quat_ned2eci(Eigen::Vector3d pos_eci, double t);
  
  /**
   * @brief Get quaternion for NED to Body transformation
   * @param quat_eci2body Quaternion from ECI to Body frame
   * @param pos_eci Position vector in ECI frame [m]
   * @param t Time [s]
   * @return Quaternion representing rotation from NED to Body
   */
  static Eigen::Quaterniond quat_ned2body(Eigen::Quaterniond quat_eci2body,
                                          Eigen::Vector3d pos_eci, double t);
  
  /**
   * @brief Create quaternion from Euler angles in degrees
   * @param az_deg Azimuth angle [deg]
   * @param el_deg Elevation angle [deg]
   * @param ro_deg Roll angle [deg]
   * @return Quaternion representing the rotation
   */
  static Eigen::Quaterniond quat_from_euler_deg(double az_deg, double el_deg,
                                                double ro_deg);
  
  /**
   * @brief Extract Euler angles from quaternion
   * @param q Input quaternion
   * @return Vector3d containing [azimuth, elevation, roll] in radians
   */
  static Eigen::Vector3d euler_from_quat(Eigen::Quaterniond q);
  
  /**
   * @brief Extract Euler angles from direction cosine matrix
   * @param C Direction cosine matrix (DCM)
   * @return Vector3d containing [azimuth, elevation, roll] in radians
   */
  static Eigen::Vector3d euler_from_dcm(Eigen::Matrix3d C);
  
  /**
   * @brief Convert quaternion to direction cosine matrix
   * @param q Input quaternion
   * @return Direction cosine matrix (DCM)
   */
  static Eigen::Matrix3d dcm_from_quat(Eigen::Quaterniond q);
  
  /**
   * @brief Convert direction cosine matrix to quaternion
   * @param C Direction cosine matrix (DCM)
   * @return Quaternion representing the rotation
   */
  static Eigen::Quaterniond quat_from_dcm(Eigen::Matrix3d C);
  
  /**
   * @brief Create DCM from thrust vector
   * @param thrustvec_eci Thrust vector in ECI frame
   * @param pos_eci Position vector in ECI frame [m]
   * @return Direction cosine matrix aligned with thrust vector
   */
  static Eigen::Matrix3d dcm_from_thrustvector(Eigen::Vector3d thrustvec_eci,
                                               Eigen::Vector3d pos_eci);
  
  /**
   * @brief Create quaternion from thrust vector
   * @param thrustvec_eci Thrust vector in ECI frame
   * @param pos_eci Position vector in ECI frame [m]
   * @return Quaternion aligned with thrust vector
   */
  static Eigen::Quaterniond quat_from_thrustvector(
      Eigen::Vector3d thrustvec_eci, Eigen::Vector3d pos_eci);
  
  /**
   * @brief Calculate orbital elements from position and velocity
   * @param pos_eci Position vector in ECI frame [m]
   * @param vel_eci Velocity vector in ECI frame [m/s]
   * @return 6x1 matrix containing [a, e, i, Omega, omega, nu]
   *         a: semi-major axis [m]
   *         e: eccentricity [-]
   *         i: inclination [rad]
   *         Omega: RAAN (right ascension of ascending node) [rad]
   *         omega: argument of periapsis [rad]
   *         nu: true anomaly [rad]
   */
  static Eigen::Matrix<double, 6, 1> orbital_elements(Eigen::Vector3d pos_eci,
                                                      Eigen::Vector3d vel_eci);
  
  /**
   * @brief Calculate position from orbital elements
   * @param elem 6x1 matrix containing orbital elements [a, e, i, Omega, omega, nu]
   * @return Position vector in ECI frame [m]
   */
  static Eigen::Vector3d pos_from_orbital_elements(
      Eigen::Matrix<double, 6, 1> elem);
  
  /**
   * @brief Calculate velocity from orbital elements
   * @param elem 6x1 matrix containing orbital elements [a, e, i, Omega, omega, nu]
   * @return Velocity vector in ECI frame [m/s]
   */
  static Eigen::Vector3d vel_from_orbital_elements(
      Eigen::Matrix<double, 6, 1> elem);
};

#endif  // SRC_COORDINATE_HPP_
