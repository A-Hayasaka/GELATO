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

#include "wrapper_air.hpp"
#include "wrapper_coordinate.hpp"
#include "wrapper_utils.hpp"

vec3d dynamics_velocity(double mass_e, vec3d pos_eci_e, vec3d vel_eci_e,
                        vec4d quat_eci2body, double t, vecXd param,
                        matXd wind_table, matXd CA_table, vecXd units) {
  double mass = mass_e * units[0];
  vec3d pos_eci = pos_eci_e * units[1];
  vec3d vel_eci = vel_eci_e * units[2];

  double thrust_vac = param[0];
  double air_area = param[2];
  double nozzle_area = param[4];

  vec3d aeroforce_eci(0.0, 0.0, 0.0);
  double thrust = thrust_vac;
  bool use_air = (air_area > 0.0);


  if (use_air) {
    vec3d pos_llh = ecef2geodetic(pos_eci[0], pos_eci[1], pos_eci[2]);
    double altitude = geopotential_altitude(pos_llh[2]);
    double rho = airdensity_at(altitude);
    double p = airpressure_at(altitude);

    vec3d vel_ecef = vel_eci2ecef(vel_eci, pos_eci, t);
    vec3d vel_wind_ned = wind_ned(altitude, wind_table);

    vec3d vel_wind_eci =
        quatrot(quat_nedg2eci(pos_eci, t), vel_wind_ned);
    vec3d vel_air_eci = ecef2eci(vel_ecef, t) - vel_wind_eci;
    double mach_number = vel_air_eci.norm() / speed_of_sound(altitude);

    double ca = interp(mach_number, CA_table.col(0), CA_table.col(1));

    aeroforce_eci =
        0.5 * rho * air_area * ca * vel_air_eci.norm() * -vel_air_eci;

    thrust = thrust_vac - nozzle_area * p;
  } else {
    aeroforce_eci.setZero();
    thrust = thrust_vac;
  }
  vec3d thrustdir_eci =
      quatrot(conj(quat_eci2body), vec3d(1.0, 0.0, 0.0));
  vec3d thrust_eci = thrust * thrustdir_eci;
  vec3d gravity_eci = gravity(pos_eci);
  vec3d acc_eci = (thrust_eci + aeroforce_eci) / mass + gravity_eci;

  return acc_eci / units[2];
}

matXd dynamics_velocity_array(vecXd mass_e, matXd pos_eci_e, matXd vel_eci_e,
                        matXd quat_eci2body, vecXd t, vecXd param,
                        matXd wind_table, matXd CA_table, vecXd units) {

  matXd acc_eci_e = matXd::Zero(pos_eci_e.rows(), 3);
  for (int i = 0; i < pos_eci_e.rows(); i++) {
    acc_eci_e.row(i) = dynamics_velocity(mass_e(i), pos_eci_e.row(i), vel_eci_e.row(i),
                        quat_eci2body.row(i), t(i), param,
                        wind_table, CA_table, units);
  }

  return acc_eci_e;
}

vec4d dynamics_quaternion(vec4d quat_eci2body, vec3d u_e, double unit_u) {

  vec3d omega_rps_body = u_e * unit_u * M_PI / 180.0;

  vec4d d_quat = vec4d::Zero();
  d_quat[0] = 0.5 * (-quat_eci2body[1] * omega_rps_body[0] -
                      quat_eci2body[2] * omega_rps_body[1] -
                      quat_eci2body[3] * omega_rps_body[2]);
  d_quat[1] = 0.5 * ( quat_eci2body[0] * omega_rps_body[0] -
                      quat_eci2body[3] * omega_rps_body[1] +
                      quat_eci2body[2] * omega_rps_body[2]);
  d_quat[2] = 0.5 * ( quat_eci2body[3] * omega_rps_body[0] +
                      quat_eci2body[0] * omega_rps_body[1] -
                      quat_eci2body[1] * omega_rps_body[2]);
  d_quat[3] = 0.5 * (-quat_eci2body[2] * omega_rps_body[0] +
                      quat_eci2body[1] * omega_rps_body[1] +
                      quat_eci2body[0] * omega_rps_body[2]);

  return d_quat;
}


matXd dynamics_quaternion_array(matXd quat_eci2body, matXd u_e, double unit_u) {

  matXd d_quat = matXd::Zero(quat_eci2body.rows(), 4);

  for (int i = 0; i < quat_eci2body.rows(); i++) {
    d_quat.row(i) = dynamics_quaternion(quat_eci2body.row(i), u_e.row(i), unit_u);
  }

  return d_quat;
}

py::dict dynamics_velocity_rh_gradient(
  double mass, vec3d pos, vec3d vel, vec4d quat, double t,
  vecXd param, matXd wind_table, matXd CA_table, vecXd units,
  double to, double tf, double unit_time, double dx
) {


  vec3d f_c = dynamics_velocity(
    mass, pos, vel, quat, t, param, wind_table, CA_table, units
  );

  // mass
  vec3d grad_mass = vec3d::Zero();
  mass += dx;
  vec3d f_p_mass = dynamics_velocity(
    mass, pos, vel, quat, t, param, wind_table, CA_table, units
  );
  mass -= dx;
  grad_mass = -(f_p_mass - f_c) / dx * (tf - to) * unit_time / 2.0;

  // position
  mat3d grad_position = mat3d::Zero();
  for (int k = 0; k < 3; k++) {
    pos(k) += dx;
    vec3d f_p_pos = dynamics_velocity(
      mass, pos, vel, quat, t, param, wind_table, CA_table, units
    );
    pos(k) -= dx;
    grad_position.col(k) = -(f_p_pos - f_c) / dx * (tf - to) * unit_time / 2.0;
  }

  // velocity: changes only affect aerodynamic forces
  mat3d grad_velocity = mat3d::Zero();
  if (param[2] > 0.0) {
    for (int k = 0; k < 3; k++) {
      vel(k) += dx;
      vec3d f_p_vel = dynamics_velocity(
        mass, pos, vel, quat, t, param, wind_table, CA_table, units
      );
      vel(k) -= dx;
      grad_velocity.col(k) = -(f_p_vel - f_c) / dx * (tf - to) * unit_time / 2.0;
    }
  }

  // quaternion
  matXd grad_quaternion = matXd::Zero(3, 4);
  for (int k = 0; k < 4; k++) {
    quat(k) += dx;
    vec3d f_p_quat = dynamics_velocity(
      mass, pos, vel, quat, t, param, wind_table, CA_table, units
    );
    quat(k) -= dx;
    grad_quaternion.col(k) = -(f_p_quat - f_c) / dx * (tf - to) * unit_time / 2.0;
  }

  // to, tf: changes only affect aerodynamic forces
  vec3d grad_to = vec3d::Zero();
  vec3d grad_tf = vec3d::Zero();
  if (param[2] > 0.0) {
    double to_p = to + dx;
    double t_p1 = to_p + t / (tf - to) * (tf - to_p);
    vec3d f_p_to = dynamics_velocity(
      mass, pos, vel, quat, t_p1, param, wind_table, CA_table, units
    );
    grad_to = -(f_p_to * (tf - to_p) - f_c * (tf - to)) / dx * unit_time / 2.0;

    double tf_p = tf + dx;
    double t_p2 = to + t / (tf - to) * (tf_p - to);
    vec3d f_p_tf = dynamics_velocity(
      mass, pos, vel, quat, t_p2, param, wind_table, CA_table, units
    );
    grad_tf = -(f_p_tf * (tf_p - to) - f_c * (tf - to)) / dx * unit_time / 2.0;
  } else {
    grad_to = f_c * unit_time / 2.0;
    grad_tf = -grad_to;
  }

  py::dict grad;
  grad["mass"] = grad_mass;
  grad["position"] = grad_position;
  grad["velocity"] = grad_velocity;
  grad["quaternion"] = grad_quaternion;
  grad["to"] = grad_to;
  grad["tf"] = grad_tf;

  return grad;
}

py::dict dynamics_quaternion_rh_gradient(
  vec4d quat, vec3d u, double unit_u,
  double to, double tf, double unit_time, double dx
) {

  vec4d f_c = dynamics_quaternion(quat, u, unit_u);

  vec3d omega_rps_body = u * unit_u * M_PI / 180.0;

  // quaternion
  matXd grad_quaternion = matXd::Zero(4, 4);
  grad_quaternion << 0.0, -0.5 * omega_rps_body[0], -0.5 * omega_rps_body[1], -0.5 * omega_rps_body[2],
                     0.5 * omega_rps_body[0], 0.0, 0.5 * omega_rps_body[2], -0.5 * omega_rps_body[1],
                     0.5 * omega_rps_body[1], -0.5 * omega_rps_body[2], 0.0, 0.5 * omega_rps_body[0],
                     0.5 * omega_rps_body[2], 0.5 * omega_rps_body[1], -0.5 * omega_rps_body[0], 0.0;
  grad_quaternion = -grad_quaternion * (tf - to) * unit_time / 2.0;

  // u (angular velocity)
  matXd grad_u = matXd::Zero(4, 3);
  grad_u << -0.5 * quat[1], -0.5 * quat[2], -0.5 * quat[3],
             0.5 * quat[0], -0.5 * quat[3], 0.5 * quat[2],
             0.5 * quat[3], 0.5 * quat[0], -0.5 * quat[1],
             -0.5 * quat[2], 0.5 * quat[1], 0.5 * quat[0];
  grad_u = -grad_u * unit_u * M_PI / 180.0 * (tf - to) * unit_time / 2.0;

  //to, tf
  vec4d grad_to = f_c * unit_time / 2.0;
  vec4d grad_tf = -grad_to;

  py::dict grad;
  grad["quaternion"] = grad_quaternion;
  grad["u"] = grad_u;
  grad["to"] = grad_to;
  grad["tf"] = grad_tf;

  return grad;
}

PYBIND11_MODULE(dynamics_c, m) {
  m.def("dynamics_velocity_array", &dynamics_velocity_array,
        "velocity with aerodynamic forces (array)");
  m.def("dynamics_quaternion_array", &dynamics_quaternion_array,
        "quaternion dynamics array");
  m.def("dynamics_velocity_rh_gradient", &dynamics_velocity_rh_gradient,
        "gradient of equality_dynamics_velocity RHS components");
  m.def("dynamics_quaternion_rh_gradient", &dynamics_quaternion_rh_gradient,
        "gradient of equality_dynamics_quaternion RHS components");
}
