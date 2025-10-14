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
 * @file wrapper_coordinate.hpp
 * @brief Python binding wrapper functions for coordinate transformations
 * 
 * This header provides wrapper functions for coordinate transformations,
 * quaternion operations, and orbital mechanics calculations to facilitate
 * Python bindings via PyBind11. Functions are designed to work with
 * NumPy arrays through Eigen matrices.
 * 
 * Includes wrappers for:
 * - Quaternion operations (multiplication, conjugate, rotation)
 * - Coordinate transformations (ECEF, ECI, NED, geodetic)
 * - Direction cosine matrix conversions
 * - Orbital element calculations
 * - Gravity calculations
 */

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
#include <utility>

#include "Coordinate.hpp"
#include "Earth.hpp"
#include "gravity.hpp"

namespace py = pybind11;

#ifndef SRC_WRAPPER_COORDINATE_HPP_
#define SRC_WRAPPER_COORDINATE_HPP_

using vec3d = Eigen::Matrix<double, 3, 1>;   ///< 3D vector type
using vec4d = Eigen::Matrix<double, 4, 1>;   ///< 4D vector type (for quaternions)
using vecXd = Eigen::Matrix<double, -1, 1>;  ///< Dynamic size vector type
using matXd = Eigen::Matrix<double, -1, -1, Eigen::RowMajor>;  ///< Dynamic size matrix type
using mat3d = Eigen::Matrix<double, 3, 3, Eigen::RowMajor>;    ///< 3x3 matrix type

// ============================================================================
// Quaternion Operations
// ============================================================================

/**
 * @brief Multiply two quaternions
 * @param q First quaternion [w, x, y, z]
 * @param p Second quaternion [w, x, y, z]
 * @return Product quaternion q*p
 */
vec4d quatmult(vec4d q, vec4d p) {
  vec4d qp;
  qp[0] = q[0] * p[0] - q[1] * p[1] - q[2] * p[2] - q[3] * p[3];
  qp[1] = q[0] * p[1] + q[1] * p[0] + q[2] * p[3] - q[3] * p[2];
  qp[2] = q[0] * p[2] - q[1] * p[3] + q[2] * p[0] + q[3] * p[1];
  qp[3] = q[0] * p[3] + q[1] * p[2] - q[2] * p[1] + q[3] * p[0];
  return qp;
}

/**
 * @brief Calculate conjugate of a quaternion
 * @param q Input quaternion [w, x, y, z]
 * @return Conjugate quaternion [w, -x, -y, -z]
 */
vec4d conj(vec4d q) {
  vec4d q_conj;
  q_conj[0] = q[0];
  q_conj[1] = -q[1];
  q_conj[2] = -q[2];
  q_conj[3] = -q[3];
  return q_conj;
}

/**
 * @brief Normalize a vector
 * @param v Input vector
 * @return Normalized vector (unit vector)
 */
vecXd normalize(vecXd v) { return v / v.norm(); }

/**
 * @brief Rotate a 3D vector by a quaternion
 * @param q Rotation quaternion [w, x, y, z]
 * @param v Vector to rotate
 * @return Rotated vector
 */
vec3d quatrot(vec4d q, vec3d v) {
  vec4d vq;
  vq[0] = 0;
  vq[1] = v[0];
  vq[2] = v[1];
  vq[3] = v[2];
  vec4d rvq = quatmult(conj(q), quatmult(vq, q));
  return rvq.tail(3);
}

// ============================================================================
// Direction Cosine Matrix (DCM) Operations
// ============================================================================

/**
 * @brief Convert quaternion to direction cosine matrix
 * @param q Quaternion [w, x, y, z]
 * @return 3x3 direction cosine matrix
 */
Eigen::Matrix3d dcm_from_quat(vec4d q) {
  Eigen::Matrix3d C;
  C(0, 0) = q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3];
  C(0, 1) = 2 * (q[1] * q[2] + q[0] * q[3]);
  C(0, 2) = 2 * (q[1] * q[3] - q[0] * q[2]);

  C(1, 0) = 2 * (q[1] * q[2] - q[0] * q[3]);
  C(1, 1) = q[0] * q[0] - q[1] * q[1] + q[2] * q[2] - q[3] * q[3];
  C(1, 2) = 2 * (q[2] * q[3] + q[0] * q[1]);

  C(2, 0) = 2 * (q[1] * q[3] + q[0] * q[2]);
  C(2, 1) = 2 * (q[2] * q[3] - q[0] * q[1]);
  C(2, 2) = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3];
  return C;
}

/**
 * @brief Convert direction cosine matrix to quaternion
 * @param C 3x3 direction cosine matrix
 * @return Quaternion [w, x, y, z]
 */
vec4d quat_from_dcm(Eigen::Matrix3d C) {
  vec4d q;
  q[0] = 0.5 * sqrt(1 + C(0, 0) + C(1, 1) + C(2, 2));
  q[1] = (C(1, 2) - C(2, 1)) / (4 * q[0]);
  q[2] = (C(2, 0) - C(0, 2)) / (4 * q[0]);
  q[3] = (C(0, 1) - C(1, 0)) / (4 * q[0]);
  return q;
}

// ============================================================================
// Coordinate Conversions (ECEF, ECI, Geodetic)
// ============================================================================

/**
 * @brief Convert ECEF coordinates to geodetic coordinates
 * @param x ECEF X coordinate [m]
 * @param y ECEF Y coordinate [m]
 * @param z ECEF Z coordinate [m]
 * @return Geodetic coordinates (latitude [deg], longitude [deg], altitude [m])
 */
vec3d ecef2geodetic(double x, double y, double z) {
  vec3d pos_ecef(x, y, z);
  vec3d geodetic = Earth::ecef2geodetic(pos_ecef);
  geodetic[0] = geodetic[0] * 180.0 / M_PI;
  geodetic[1] = geodetic[1] * 180.0 / M_PI;
  return geodetic;
}

/**
 * @brief Convert geodetic coordinates to ECEF coordinates
 * @param lat Latitude [deg]
 * @param lon Longitude [deg]
 * @param alt Altitude [m]
 * @return ECEF position vector [m]
 */
vec3d geodetic2ecef(double lat, double lon, double alt) {
  vec3d geodetic(lat * M_PI / 180.0, lon * M_PI / 180.0, alt);
  return Earth::geodetic2ecef(geodetic);
}

/**
 * @brief Convert ECEF position to ECI position
 * @param xyz_in ECEF position vector [m]
 * @param t Time [s]
 * @return ECI position vector [m]
 */
vec3d ecef2eci(vec3d xyz_in, double t) {
  return Coordinate::ecef2eci(xyz_in, t);
}

/**
 * @brief Convert ECI position to ECEF position
 * @param xyz_in ECI position vector [m]
 * @param t Time [s]
 * @return ECEF position vector [m]
 */
vec3d eci2ecef(vec3d xyz_in, double t) {
  return Coordinate::eci2ecef(xyz_in, t);
}

/**
 * @brief Convert ECEF ground velocity to ECI inertial velocity
 * @param vel_ecef Ground velocity in ECEF frame [m/s]
 * @param pos_ecef Position in ECEF frame [m]
 * @param t Time [s]
 * @return Inertial velocity in ECI frame [m/s]
 */
vec3d vel_ecef2eci(vec3d vel_ecef, vec3d pos_ecef, double t) {
  return Coordinate::vel_ecef2eci(vel_ecef, pos_ecef, t);
}

/**
 * @brief Convert ECI inertial velocity to ECEF ground velocity
 * @param vel_eci Inertial velocity in ECI frame [m/s]
 * @param pos_eci Position in ECI frame [m]
 * @param t Time [s]
 * @return Ground velocity in ECEF frame [m/s]
 */
vec3d vel_eci2ecef(vec3d vel_eci, vec3d pos_eci, double t) {
  return Coordinate::vel_eci2ecef(vel_eci, pos_eci, t);
}

/**
 * @brief Get quaternion for ECI to ECEF transformation
 * @param t Time [s]
 * @return Quaternion [w, x, y, z]
 */
vec4d quat_eci2ecef(double t) {
  Eigen::Quaterniond q = Coordinate::quat_eci2ecef(t);
  return vec4d(q.w(), q.x(), q.y(), q.z());
}

/**
 * @brief Get quaternion for ECEF to ECI transformation
 * @param t Time [s]
 * @return Quaternion [w, x, y, z]
 */
vec4d quat_ecef2eci(double t) {
  Eigen::Quaterniond q = Coordinate::quat_ecef2eci(t);
  return vec4d(q.w(), q.x(), q.y(), q.z());
}

/**
 * @brief Get quaternion for ECEF to NED (geodetic) transformation
 * @param pos_ecef Position in ECEF frame [m]
 * @return Quaternion [w, x, y, z]
 */
vec4d quat_ecef2nedg(vec3d pos_ecef) {
  Eigen::Quaterniond q = Coordinate::quat_ecef2ned(pos_ecef);
  return vec4d(q.w(), q.x(), q.y(), q.z());
}

/**
 * @brief Get quaternion for NED (geodetic) to ECEF transformation
 * @param pos_ecef Position in ECEF frame [m]
 * @return Quaternion [w, x, y, z]
 */
vec4d quat_nedg2ecef(vec3d pos_ecef) {
  Eigen::Quaterniond q = Coordinate::quat_ned2ecef(pos_ecef);
  return vec4d(q.w(), q.x(), q.y(), q.z());
}

/**
 * @brief Get quaternion for ECI to NED (geodetic) transformation
 * @param pos_eci Position in ECI frame [m]
 * @param t Time [s]
 * @return Quaternion [w, x, y, z]
 */
vec4d quat_eci2nedg(vec3d pos_eci, double t) {
  Eigen::Quaterniond q = Coordinate::quat_eci2ned(pos_eci, t);
  return vec4d(q.w(), q.x(), q.y(), q.z());
}

/**
 * @brief Get quaternion for NED (geodetic) to ECI transformation
 * @param pos_eci Position in ECI frame [m]
 * @param t Time [s]
 * @return Quaternion [w, x, y, z]
 */
vec4d quat_nedg2eci(vec3d pos_eci, double t) {
  Eigen::Quaterniond q = Coordinate::quat_ned2eci(pos_eci, t);
  return vec4d(q.w(), q.x(), q.y(), q.z());
}

/**
 * @brief Create quaternion from Euler angles
 * @param az Azimuth angle [deg]
 * @param el Elevation angle [deg]
 * @param ro Roll angle [deg]
 * @return Quaternion [w, x, y, z]
 */
vec4d quat_from_euler(double az, double el, double ro) {
  Eigen::Quaterniond q = Coordinate::quat_from_euler_deg(az, el, ro);
  return vec4d(q.w(), q.x(), q.y(), q.z());
}

/**
 * @brief Calculate gravity acceleration in ECI frame
 * @param pos Position in ECI frame [m]
 * @return Gravity acceleration vector [m/s²]
 */
vec3d gravity(vec3d pos) { return gravityECI(pos); }

/**
 * @brief Get quaternion for NED (geodetic) to body transformation
 * @param quat_eci2body Quaternion from ECI to body frame [w, x, y, z]
 * @param pos_eci Position in ECI frame [m]
 * @param t Time [s]
 * @return Quaternion from NED (geodetic) to body [w, x, y, z]
 */
vec4d quat_nedg2body(vec4d quat_eci2body, vec3d pos_eci, double t) {
  vec4d q = quat_eci2nedg(pos_eci, t);
  return quatmult(conj(q), quat_eci2body);
}

/**
 * @brief Extract Euler angles from quaternion
 * @param q Quaternion [w, x, y, z]
 * @return Euler angles (azimuth, elevation, roll) [deg]
 */
vec3d euler_from_quat(vec4d q) {
  Eigen::Quaterniond q_(q[0], q[1], q[2], q[3]);
  vec3d euler = Coordinate::euler_from_quat(q_);
  return euler * 180.0 / M_PI;
}

/**
 * @brief Extract Euler angles from direction cosine matrix
 * @param C Direction cosine matrix
 * @return Euler angles (azimuth, elevation, roll) [deg]
 */
vec3d euler_from_dcm(Eigen::Matrix3d C) {
  Eigen::Matrix3d C_(C.data());
  vec3d euler = Coordinate::euler_from_dcm(C_);
  return euler * 180.0 / M_PI;
}

/**
 * @brief Create DCM aligned with thrust vector
 * @param pos_eci Position in ECI frame [m]
 * @param thrustvec_eci Thrust vector in ECI frame
 * @return Direction cosine matrix
 */
Eigen::Matrix3d dcm_from_thrustvector(vec3d pos_eci, vec3d thrustvec_eci) {
  Eigen::Matrix3d C = Coordinate::dcm_from_thrustvector(thrustvec_eci, pos_eci);
  return C;
}

/**
 * @brief Convert ECI position to geodetic coordinates
 * @param pos_eci Position in ECI frame [m]
 * @param t Time [s]
 * @return Geodetic coordinates (latitude [deg], longitude [deg], altitude [m])
 */
vec3d eci2geodetic(vec3d pos_eci, double t) {
  vec3d pos_ecef = Coordinate::eci2ecef(pos_eci, t);
  vec3d geodetic = Earth::ecef2geodetic(pos_ecef);
  geodetic[0] = geodetic[0] * 180.0 / M_PI;
  geodetic[1] = geodetic[1] * 180.0 / M_PI;
  return geodetic;
}

// ============================================================================
// Orbital Mechanics
// ============================================================================

/**
 * @brief Calculate orbital elements from position and velocity
 * @param pos_eci Position in ECI frame [m]
 * @param vel_eci Velocity in ECI frame [m/s]
 * @return 6x1 vector containing [a, e, i, Omega, omega, nu]
 *         where angles (i, Omega, omega, nu) are in degrees
 */
Eigen::Matrix<double, 6, 1> orbital_elements(vec3d pos_eci, vec3d vel_eci) {
  Eigen::Matrix<double, 6, 1> elem =
      Coordinate::orbital_elements(pos_eci, vel_eci);
  elem[2] = elem[2] * 180.0 / M_PI;
  elem[3] = elem[3] * 180.0 / M_PI;
  elem[4] = elem[4] * 180.0 / M_PI;
  elem[5] = elem[5] * 180.0 / M_PI;
  return elem;
}

/**
 * @brief Calculate geodesic distance using Vincenty formula
 * @param lat_origin Origin latitude [deg]
 * @param lon_origin Origin longitude [deg]
 * @param lat_target Target latitude [deg]
 * @param lon_target Target longitude [deg]
 * @return Geodesic distance [m]
 */
double distance_vincenty(double lat_origin, double lon_origin,
                         double lat_target, double lon_target) {
  std::pair<double, double> dist_azimuth =
      Earth::distance_vincenty(Eigen::Vector3d(lat_origin * M_PI / 180.0,
                                               lon_origin * M_PI / 180.0, 0.0),
                               Eigen::Vector3d(lat_target * M_PI / 180.0,
                                               lon_target * M_PI / 180.0, 0.0));

  return dist_azimuth.first;
}

/**
 * @brief Calculate angular momentum vector
 * @param pos_eci Position in ECI frame [m]
 * @param vel_eci Velocity in ECI frame [m/s]
 * @return Angular momentum vector [m²/s]
 */
vec3d angular_momentum_vec(vec3d pos_eci, vec3d vel_eci) {
  return pos_eci.cross(vel_eci);
}

/**
 * @brief Calculate angular momentum magnitude
 * @param pos_eci Position in ECI frame [m]
 * @param vel_eci Velocity in ECI frame [m/s]
 * @return Angular momentum magnitude [m²/s]
 */
double angular_momentum(vec3d pos_eci, vec3d vel_eci) {
  return angular_momentum_vec(pos_eci, vel_eci).norm();
}

/**
 * @brief Calculate cosine of orbital inclination
 * @param pos_eci Position in ECI frame [m]
 * @param vel_eci Velocity in ECI frame [m/s]
 * @return Cosine of inclination [-]
 */
double inclination_cosine(vec3d pos_eci, vec3d vel_eci) {
  return angular_momentum_vec(pos_eci, vel_eci)[2] /
         angular_momentum(pos_eci, vel_eci);
}

/**
 * @brief Calculate orbital inclination
 * @param pos_eci Position in ECI frame [m]
 * @param vel_eci Velocity in ECI frame [m/s]
 * @return Inclination [rad]
 */
double inclination_rad(vec3d pos_eci, vec3d vel_eci) {
  return acos(inclination_cosine(pos_eci, vel_eci));
}

/**
 * @brief Calculate Laplace eccentricity vector
 * @param pos_eci Position in ECI frame [m]
 * @param vel_eci Velocity in ECI frame [m/s]
 * @return Laplace vector (eccentricity vector) [-]
 */
vec3d laplace_vector(vec3d pos_eci, vec3d vel_eci) {
  vec3d h = angular_momentum_vec(pos_eci, vel_eci);
  vec3d r = pos_eci;
  vec3d v = vel_eci;
  vec3d e = v.cross(h) / Earth::mu - r / r.norm();
  return e;
}

/**
 * @brief Calculate orbital energy (specific)
 * @param pos_eci Position in ECI frame [m]
 * @param vel_eci Velocity in ECI frame [m/s]
 * @return Specific orbital energy [J/kg]
 */
double orbit_energy(vec3d pos_eci, vec3d vel_eci) {
  double r = pos_eci.norm();
  double v = vel_eci.norm();
  return 0.5 * v * v - Earth::mu / r;
}

/**
 * @brief Calculate angular momentum from apogee and perigee altitudes
 * @param ha Apogee altitude [m]
 * @param hp Perigee altitude [m]
 * @return Angular momentum magnitude [m²/s]
 */
double angular_momentum_from_altitude(double ha, double hp) {
  double ra = Earth::Ra + ha;
  double rp = Earth::Ra + hp;
  double a = (ra + rp) / 2.0;
  double vp = sqrt(Earth::mu * (2.0 / rp - 1.0 / a));
  return rp * vp;
}

/**
 * @brief Calculate orbital energy from apogee and perigee altitudes
 * @param ha Apogee altitude [m]
 * @param hp Perigee altitude [m]
 * @return Specific orbital energy [J/kg]
 */
double orbit_energy_from_altitude(double ha, double hp) {
  double ra = Earth::Ra + ha;
  double rp = Earth::Ra + hp;
  double a = (ra + rp) / 2.0;
  return -Earth::mu / 2.0 / a;
}

#endif  // SRC_WRAPPER_COORDINATE_HPP_
