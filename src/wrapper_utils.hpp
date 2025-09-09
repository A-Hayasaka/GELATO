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

#include <algorithm>
#include <vector>

#include "wrapper_air.hpp"
#include "wrapper_coordinate.hpp"

namespace py = pybind11;

#ifndef SRC_WRAPPER_UTILS_HPP_
#define SRC_WRAPPER_UTILS_HPP_

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

double interp(double x, vecXd xp, vecXd yp) {
  // Linear interpolation
  // x: value to interpolate
  // xp: x values of the data points
  // yp: y values of the data points

  const int n = xp.size();

  if (x <= xp(0)) {
    return yp(0);
  }

  if (x >= xp(n - 1)) {
    return yp(n - 1);
  }

  // Find index i such that xp[i] <= x < xp[i+1]
  auto upper = std::upper_bound(xp.data(), xp.data() + n, x);
  int idx = static_cast<int>(upper - xp.data()) - 1;
  if (idx < 0) idx = 0;
  if (idx >= n - 1) idx = n - 2;

  const double x_lower = xp(idx);
  const double x_upper = xp(idx + 1);
  const double y_lower = yp(idx);
  const double y_upper = yp(idx + 1);

  const double dx = x_upper - x_lower;
  if (std::abs(dx) < 1e-12) {
    // Avoid divide-by-zero if xp has duplicate values
    return y_lower;
  }
  const double alpha = (x - x_lower) / dx;

  return y_lower + alpha * (y_upper - y_lower);
}

vec3d wind_ned(double altitude_m, matXd wind_data) {
  double wind_u = interp(altitude_m, wind_data.col(0), wind_data.col(1));
  double wind_v = interp(altitude_m, wind_data.col(0), wind_data.col(2));

  return vec3d(wind_u, wind_v, 0.0);
}

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

double angle_of_attack_all_dimless(vec3d pos_eci, vec3d vel_eci, vec4d quat,
                                  double t, matXd wind, vecXd units) {
  return angle_of_attack_all_rad(pos_eci * units[0], vel_eci * units[1], quat,
                                 t * units[2], wind) / units[3];
}

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

vecXd angle_of_attack_all_array_dimless(matXd pos_eci_e, matXd vel_eci_e,
                                       matXd quat, vecXd t_e, matXd wind,
                                       vecXd units) {
  int n = pos_eci_e.rows();
  vecXd alpha(n);

  for (int i = 0; i < n; i++) {
    alpha(i) = angle_of_attack_all_dimless(pos_eci_e.row(i), vel_eci_e.row(i),
                                          quat.row(i), t_e(i), wind, units);
  }

  return alpha;
}

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

matXd angle_of_attack_ab_array_rad(matXd pos_eci, matXd vel_eci, matXd quat,
                                   vecXd t, matXd wind) {
  int n = pos_eci.rows();
  matXd alpha(n, 2);

  for (int i = 0; i < n; i++) {
    alpha.row(i) = angle_of_attack_ab_rad(pos_eci.row(i), vel_eci.row(i),
                                          quat.row(i), t(i), wind);
  }

  return alpha;
}

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

double dynamic_pressure_dimless(vec3d pos_eci, vec3d vel_eci, double t,
                       matXd wind, vecXd units) {
  return dynamic_pressure_pa(pos_eci * units[0], vel_eci * units[1],
                        t * units[2], wind) /
         units[3];
}

vecXd dynamic_pressure_array_pa(matXd pos_eci, matXd vel_eci, vecXd t,
                                matXd wind) {
  int n = pos_eci.rows();
  vecXd q(n);

  for (int i = 0; i < n; i++) {
    q(i) = dynamic_pressure_pa(pos_eci.row(i), vel_eci.row(i), t(i), wind);
  }

  return q;
}

vecXd dynamic_pressure_array_dimless(matXd pos_eci_e, matXd vel_eci_e,
                                   vecXd t_e, matXd wind, vecXd units) {
  int n = pos_eci_e.rows();
  vecXd q(n);

  for (int i = 0; i < n; i++) {
    q(i) = dynamic_pressure_dimless(pos_eci_e.row(i), vel_eci_e.row(i), t_e(i),
                                   wind, units);
  }

  return q;
}

double q_alpha_pa_rad(vec3d pos_eci, vec3d vel_eci, vec4d quat, double t,
                      matXd wind) {
  double alpha = angle_of_attack_all_rad(pos_eci, vel_eci, quat, t, wind);
  double q = dynamic_pressure_pa(pos_eci, vel_eci, t, wind);
  return q * alpha;
}

double q_alpha_dimless(vec3d pos_eci, vec3d vel_eci, vec4d quat, double t,
                       matXd wind, vecXd units) {
  return q_alpha_pa_rad(pos_eci * units[0], vel_eci * units[1], quat,
                        t * units[2], wind) /
         units[3];
}

vecXd q_alpha_array_pa_rad(matXd pos_eci, matXd vel_eci, matXd quat, vecXd t,
                           matXd wind) {
  int n = pos_eci.rows();
  vecXd q_alpha(n);

  for (int i = 0; i < n; i++) {
    q_alpha(i) =
        q_alpha_pa_rad(pos_eci.row(i), vel_eci.row(i), quat.row(i), t(i), wind);
  }

  return q_alpha;
}

vecXd q_alpha_array_dimless(matXd pos_eci_e, matXd vel_eci_e, matXd quat,
                            vecXd t_e, matXd wind, vecXd units) {
  int n = pos_eci_e.rows();
  vecXd q_alpha(n);

  for (int i = 0; i < n; i++) {
    q_alpha(i) = q_alpha_dimless(pos_eci_e.row(i), vel_eci_e.row(i),
                                 quat.row(i), t_e(i), wind, units);
  }

  return q_alpha;
}

py::dict angle_of_attack_all_gradient_dimless_core(int n, matXd pos_ki, matXd vel_ki,
                                       matXd quat_ki, vecXd t_ki, matXd wind,
                                       vecXd units, double dx) {
  matXd grad_pos = Eigen::MatrixXd::Zero(n, 3);
  matXd grad_vel = Eigen::MatrixXd::Zero(n, 3);
  matXd grad_quat = Eigen::MatrixXd::Zero(n, 4);
  vecXd grad_to = Eigen::VectorXd::Zero(n);
  vecXd grad_tf = Eigen::VectorXd::Zero(n);

  double to = t_ki(0);
  double tf = 0.0;
  if (n > 1)
    tf = t_ki(t_ki.size() - 1);

  double f_c, t_po, t_pf;

  for (int i = 0; i < n; i++) {
    f_c = angle_of_attack_all_dimless(pos_ki.row(i), vel_ki.row(i), quat_ki.row(i), t_ki(i),
                          wind, units);

    for (int j = 0; j < 3; j++) {
      pos_ki(i, j) += dx;
      grad_pos(i, j) = (angle_of_attack_all_dimless(pos_ki.row(i), vel_ki.row(i),
                                        quat_ki.row(i), t_ki(i), wind, units) -
                        f_c) /
                       dx;
      pos_ki(i, j) -= dx;
    }

    for (int j = 0; j < 3; j++) {
      vel_ki(i, j) += dx;
      grad_vel(i, j) = (angle_of_attack_all_dimless(pos_ki.row(i), vel_ki.row(i),
                                        quat_ki.row(i), t_ki(i), wind, units) -
                        f_c) /
                       dx;
      vel_ki(i, j) -= dx;
    }

    for (int j = 0; j < 4; j++) {
      quat_ki(i, j) += dx;
      grad_quat(i, j) = (angle_of_attack_all_dimless(pos_ki.row(i), vel_ki.row(i),
                                         quat_ki.row(i), t_ki(i), wind, units) -
                         f_c) /
                        dx;
      quat_ki(i, j) -= dx;
    }

    if (n > 1) {
      t_po = tf - (tf - (to + dx)) / (tf - to) * (tf - t_ki(i));
      grad_to(i) = (angle_of_attack_all_dimless(pos_ki.row(i), vel_ki.row(i), quat_ki.row(i),
                                    t_po, wind, units) -
                    f_c) /
                  dx;

      t_pf = to + ((tf + dx) - to) / (tf - to) * (t_ki(i) - to);
      grad_tf(i) = (angle_of_attack_all_dimless(pos_ki.row(i), vel_ki.row(i), quat_ki.row(i),
                                    t_pf, wind, units) -
                    f_c) /
                  dx;
    } else {
      grad_to(i) = (angle_of_attack_all_dimless(pos_ki.row(i), vel_ki.row(i), quat_ki.row(i),
                                    to + dx, wind, units) -
                    f_c) /
                  dx;
    }
  }

  py::dict grad;
  grad["position"] = grad_pos;
  grad["velocity"] = grad_vel;
  grad["quaternion"] = grad_quat;
  grad["to"] = grad_to;
  grad["tf"] = grad_tf;
  return grad;
}


py::dict dynamic_pressure_gradient_dimless_core(int n, matXd pos_ki, matXd vel_ki,
                                       vecXd t_ki, matXd wind,
                                       vecXd units, double dx) {
  matXd grad_pos = Eigen::MatrixXd::Zero(n, 3);
  matXd grad_vel = Eigen::MatrixXd::Zero(n, 3);
  vecXd grad_to = Eigen::VectorXd::Zero(n);
  vecXd grad_tf = Eigen::VectorXd::Zero(n);

  // assumption: when n = 1, t_ki is assumed to be the value of t_o.
  double to = t_ki(0);

  double tf = 0.0;
  if (n > 1)
    tf = t_ki(t_ki.size() - 1);

  double f_c, t_po, t_pf;

  for (int i = 0; i < n; i++) {
    f_c = dynamic_pressure_dimless(pos_ki.row(i), vel_ki.row(i), t_ki(i),
                          wind, units);

    for (int j = 0; j < 3; j++) {
      pos_ki(i, j) += dx;
      grad_pos(i, j) = (dynamic_pressure_dimless(pos_ki.row(i), vel_ki.row(i),
                                       t_ki(i), wind, units) -
                        f_c) /
                       dx;
      pos_ki(i, j) -= dx;
    }

    for (int j = 0; j < 3; j++) {
      vel_ki(i, j) += dx;
      grad_vel(i, j) = (dynamic_pressure_dimless(pos_ki.row(i), vel_ki.row(i),
                                       t_ki(i), wind, units) -
                        f_c) /
                       dx;
      vel_ki(i, j) -= dx;
    }

    if (n > 1) {
      t_po = tf - (tf - (to + dx)) / (tf - to) * (tf - t_ki(i));
      grad_to(i) = (dynamic_pressure_dimless(pos_ki.row(i), vel_ki.row(i),
                                    t_po, wind, units) -
                    f_c) /
                  dx;

      t_pf = to + ((tf + dx) - to) / (tf - to) * (t_ki(i) - to);
      grad_tf(i) = (dynamic_pressure_dimless(pos_ki.row(i), vel_ki.row(i),
                                    t_pf, wind, units) -
                    f_c) /
                  dx;
    } else {
      grad_to(i) = (dynamic_pressure_dimless(pos_ki.row(i), vel_ki.row(i),
                                    to + dx, wind, units) -
                    f_c) /
                  dx;
    }
  }

  py::dict grad;
  grad["position"] = grad_pos;
  grad["velocity"] = grad_vel;
  grad["to"] = grad_to;
  grad["tf"] = grad_tf;
  return grad;
}


py::dict q_alpha_gradient_dimless_core(int n, matXd pos_ki, matXd vel_ki,
                                       matXd quat_ki, vecXd t_ki, matXd wind,
                                       vecXd units, double dx) {
  matXd grad_pos = Eigen::MatrixXd::Zero(n, 3);
  matXd grad_vel = Eigen::MatrixXd::Zero(n, 3);
  matXd grad_quat = Eigen::MatrixXd::Zero(n, 4);
  vecXd grad_to = Eigen::VectorXd::Zero(n);
  vecXd grad_tf = Eigen::VectorXd::Zero(n);

  double to = t_ki(0);
  double tf = 0.0;
  if (n > 1)
    tf = t_ki(t_ki.size() - 1);

  double f_c, t_po, t_pf;

  for (int i = 0; i < n; i++) {
    f_c = q_alpha_dimless(pos_ki.row(i), vel_ki.row(i), quat_ki.row(i), t_ki(i),
                          wind, units);

    for (int j = 0; j < 3; j++) {
      pos_ki(i, j) += dx;
      grad_pos(i, j) = (q_alpha_dimless(pos_ki.row(i), vel_ki.row(i),
                                        quat_ki.row(i), t_ki(i), wind, units) -
                        f_c) /
                       dx;
      pos_ki(i, j) -= dx;
    }

    for (int j = 0; j < 3; j++) {
      vel_ki(i, j) += dx;
      grad_vel(i, j) = (q_alpha_dimless(pos_ki.row(i), vel_ki.row(i),
                                        quat_ki.row(i), t_ki(i), wind, units) -
                        f_c) /
                       dx;
      vel_ki(i, j) -= dx;
    }

    for (int j = 0; j < 4; j++) {
      quat_ki(i, j) += dx;
      grad_quat(i, j) = (q_alpha_dimless(pos_ki.row(i), vel_ki.row(i),
                                         quat_ki.row(i), t_ki(i), wind, units) -
                         f_c) /
                        dx;
      quat_ki(i, j) -= dx;
    }

    if (n > 1) {
      t_po = tf - (tf - (to + dx)) / (tf - to) * (tf - t_ki(i));
      grad_to(i) = (q_alpha_dimless(pos_ki.row(i), vel_ki.row(i), quat_ki.row(i),
                                    t_po, wind, units) -
                    f_c) /
                  dx;

      t_pf = to + ((tf + dx) - to) / (tf - to) * (t_ki(i) - to);
      grad_tf(i) = (q_alpha_dimless(pos_ki.row(i), vel_ki.row(i), quat_ki.row(i),
                                    t_pf, wind, units) -
                    f_c) /
                  dx;
    } else {
      grad_to(i) = (q_alpha_dimless(pos_ki.row(i), vel_ki.row(i), quat_ki.row(i),
                                    to + dx, wind, units) -
                    f_c) /
                  dx;
    }
  }

  py::dict grad;
  grad["position"] = grad_pos;
  grad["velocity"] = grad_vel;
  grad["quaternion"] = grad_quat;
  grad["to"] = grad_to;
  grad["tf"] = grad_tf;
  return grad;
}

double yb_r_dot(vec3d pos_eci, vec4d quat_eci2body) {
  // Returns sine of roll angles.
  vec3d yb_dir_eci = quatrot(conj(quat_eci2body), vec3d(0.0, 1.0, 0.0));
  return yb_dir_eci.dot(normalize(pos_eci));
}

vecXd roll_direction_array(matXd pos, matXd quat) {
  // Returns array of sine of roll angles for each state values.
  int n = pos.rows();
  vecXd yb_rd(n);

  for (int i = 0; i < n; i++) {
    yb_rd(i) = yb_r_dot(pos.row(i), quat.row(i));
  }

  return yb_rd;
}

py::dict roll_direction_array_gradient(
  matXd pos, matXd quat,
  double unit_pos, double dx) {

  int n = pos.rows();
  matXd grad_pos = Eigen::MatrixXd::Zero(n, 3);
  matXd grad_quat = Eigen::MatrixXd::Zero(n, 4);

  double f_c;

  for (int i = 0; i < n; i++) {
    f_c = yb_r_dot(pos.row(i), quat.row(i));

    for (int j = 0; j < 3; j++) {
      pos(i, j) += dx;
      grad_pos(i, j) = (yb_r_dot(pos.row(i), quat.row(i)) -
                        f_c) /
                       dx;
      pos(i, j) -= dx;
    }

    for (int j = 0; j < 4; j++) {
      quat(i, j) += dx;
      grad_quat(i, j) = (yb_r_dot(pos.row(i), quat.row(i)) -
                         f_c) /
                        dx;
      quat(i, j) -= dx;
    }
  }

  py::dict grad;
  grad["position"] = grad_pos;
  grad["quaternion"] = grad_quat;
  return grad;
}

#endif  // SRC_WRAPPER_UTILS_HPP_
