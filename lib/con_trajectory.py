#
# The MIT License
#
# Copyright (c) 2022 Interstellar Technologies Inc.
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files
# (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#

"""
Trajectory Constraints Module
==============================

Module for defining constraint conditions for the entire trajectory.

This module defines constraints that span the entire flight path (angular velocity limits,
mass constraints, etc.) and provides them as inequality constraints for optimization problems.

Main Features:
    * Constraints on roll, pitch, and yaw angular velocities
    * Non-negative constraints on propellant mass
    * Limits on attitude change rate
    * Smoothness constraints on control inputs

Constraint Functions:
    inequality_turn_rate: Inequality constraints on angular velocity
    inequality_propellant_mass: Inequality constraints on propellant mass
"""

# constraints_b.py
# constraints about stage mass and turn rate


import sys
import numpy as np
from .utils_c import *
from .coordinate_c import conj, quatrot, normalize


def yb_r_dot(pos_eci, quat_eci2body):
    """
    Calculate sine of roll angle for a single state.
    
    Computes the dot product between the body y-axis direction and the normalized 
    position vector, which gives the sine of the roll angle in the trajectory plane.
    
    Args:
        pos_eci (numpy.ndarray): Position vector in ECI frame [m] (3,)
        quat_eci2body (numpy.ndarray): Quaternion from ECI to body frame (4,)
    
    Returns:
        float: Sine of roll angle (ranges from -1 to 1)
    """
    yb_dir_eci = quatrot(conj(quat_eci2body), np.array([0.0, 1.0, 0.0]))
    return yb_dir_eci.dot(normalize(pos_eci))


def roll_direction_array(pos, quat):
    """
    Calculate sine of roll angles for all trajectory points.
    
    Computes roll angle sines for an array of position and quaternion states.
    
    Args:
        pos (numpy.ndarray): Array of position vectors in ECI frame [m], shape (N,3)
        quat (numpy.ndarray): Array of quaternions from ECI to body frame, shape (N,4)
    
    Returns:
        numpy.ndarray: Array of sine of roll angles for each state (N,)
    """
    return np.array([yb_r_dot(pos[i], quat[i]) for i in range(len(pos))])


def inequality_mass(xdict, pdict, unitdict, condition):
    """
    Inequality constraint enforcing non-negative propellant mass.
    
    Ensures that the remaining propellant mass in each stage never becomes negative
    during the trajectory by constraining mass to remain above the dry mass.
    
    Args:
        xdict (dict): Dictionary containing state variables with 'mass'
        pdict (dict): Dictionary with 'RocketStage' parameters including:
            - 'mass_dry': Dry mass of each stage [kg]
            - 'ignition_at', 'jettison_at': Stage event names
        unitdict (dict): Dictionary of unit scaling factors
        condition (dict): Configuration dictionary
    
    Returns:
        numpy.ndarray: Array of mass constraint violations (should be non-negative).
            Negative values indicate constraint violation.
    """

    con = []

    mass_ = xdict["mass"]
    for index, stage in pdict["RocketStage"].items():

        # read index number
        section_ig = [
            i
            for i, value in enumerate(pdict["params"])
            if value["name"] == stage["ignition_at"]
        ][0]
        section_co = [
            i
            for i, value in enumerate(pdict["params"])
            if value["name"] == stage["cutoff_at"]
        ][0]

        mass_ig = mass_[pdict["ps_params"].index_start_x(section_ig)]
        mass_co = mass_[pdict["ps_params"].index_start_x(section_co)]

        d_mass = stage["mass_propellant"]
        if stage["dropMass"] is not None:
            d_mass += sum([item["mass"] for item in stage["dropMass"].values()])
        con.append(-mass_ig + mass_co + d_mass / unitdict["mass"])

    return con


def inequality_jac_mass(xdict, pdict, unitdict, condition):
    """
    Compute Jacobian matrix of inequality_mass constraint.
    
    Calculates the partial derivatives of propellant mass constraints with respect 
    to the mass state variable.
    
    Args:
        xdict (dict): Dictionary containing state variables
        pdict (dict): Dictionary with 'RocketStage' and section parameters
        unitdict (dict): Dictionary of unit scaling factors
        condition (dict): Configuration dictionary
    
    Returns:
        dict: Jacobian matrix in COO sparse format with key 'mass', containing
            'coo' (row, col, data) and 'shape' (nRow, nCol).
    """

    jac = {}

    data = []
    row = []
    col = []

    counter = 0
    for index, stage in pdict["RocketStage"].items():
        section_ig = [
            i
            for i, value in enumerate(pdict["params"])
            if value["name"] == stage["ignition_at"]
        ][0]
        section_co = [
            i
            for i, value in enumerate(pdict["params"])
            if value["name"] == stage["cutoff_at"]
        ][0]
        data.extend([-1.0, 1.0])
        row.extend([counter, counter])
        col.extend(
            [
                pdict["ps_params"].index_start_x(section_ig),
                pdict["ps_params"].index_start_x(section_co),
            ]
        )
        counter += 1

    jac["mass"] = {
        "coo": [
            np.array(row, dtype="i4"),
            np.array(col, dtype="i4"),
            np.array(data, dtype="f8"),
        ],
        "shape": (counter, len(xdict["mass"])),
    }
    return jac


def inequality_kickturn(xdict, pdict, unitdict, condition):
    """
    Inequality constraint for minimum pitch rate during kick-turn maneuver.
    
    Enforces that pitch angular velocity remains above minimum during kick-turn 
    sections to ensure proper trajectory shaping.
    
    Args:
        xdict (dict): Dictionary containing state variables with 'u' (angular velocity)
        pdict (dict): Dictionary with section parameters indicating 'kick' attitude mode
        unitdict (dict): Dictionary of unit scaling factors with 'u'
        condition (dict): Configuration dictionary
    
    Returns:
        numpy.ndarray: Array of angular velocity constraint violations (negative values
            indicate constraint satisfaction).
    """

    con = []
    unit_u = unitdict["u"]
    u_ = xdict["u"].reshape(-1, 3) * unit_u
    num_sections = pdict["num_sections"]

    for i in range(num_sections - 1):
        # kick turn
        if "kick" in pdict["params"][i]["attitude"]:
            ua, ub, xa, xb, n = pdict["ps_params"].get_index(i)
            u_i_ = u_[ua:ub]
            con.append(-u_i_[:, 1])
            # con.append(u_i_[:,1]+0.36)

    return np.concatenate(con, axis=None)


def inequality_jac_kickturn(xdict, pdict, unitdict, condition):
    """
    Compute Jacobian matrix of inequality_kickturn constraint.
    
    Calculates the partial derivatives of kick-turn rate constraints with respect 
    to angular velocity.
    
    Args:
        xdict (dict): Dictionary containing state variables
        pdict (dict): Dictionary with section parameters
        unitdict (dict): Dictionary of unit scaling factors
        condition (dict): Configuration dictionary
    
    Returns:
        dict: Jacobian matrix in COO sparse format with key 'u', containing
            'coo' (row, col, data) and 'shape' (nRow, nCol).
    """

    jac = {}
    num_sections = pdict["num_sections"]

    data = []
    row = []
    col = []

    nRow = 0
    for i in range(num_sections - 1):

        # kick turn
        if "kick" in pdict["params"][i]["attitude"]:
            ua, ub, xa, xb, n = pdict["ps_params"].get_index(i)
            row.extend(range(nRow, nRow + n))
            col.extend(range(ua * 3 + 1, ub * 3 + 1, 3))
            data.extend([-1.0] * n)
            nRow += n

    jac["u"] = {
        "coo": [
            np.array(row, dtype="i4"),
            np.array(col, dtype="i4"),
            np.array(data, dtype="f8"),
        ],
        "shape": (nRow, len(xdict["u"])),
    }

    return jac


def equality_6DoF_rate(xdict, pdict, unitdict, condition):
    """
    Equality constraint relating roll rate to trajectory geometry.
    
    Enforces the relationship between the body-frame roll rate and the rate of change 
    of roll angle relative to the trajectory plane. This ensures consistency between 
    angular velocity commands and actual trajectory rotation.
    
    Args:
        xdict (dict): Dictionary containing state variables:
            - 'position': Position vectors (N,3)
            - 'quaternion': Quaternions (N,4)
            - 'u': Angular velocities (N,3) [rad/s]
        pdict (dict): Dictionary with section parameters and 'ps_params'
        unitdict (dict): Dictionary of unit scaling factors
        condition (dict): Configuration dictionary
    
    Returns:
        numpy.ndarray: Array of roll rate constraint violations (should be zero at optimum).
    """

    con = []

    unit_pos = unitdict["position"]

    pos_ = xdict["position"].reshape(-1, 3)
    quat_ = xdict["quaternion"].reshape(-1, 4)

    u_ = xdict["u"].reshape(-1, 3)

    num_sections = pdict["num_sections"]

    for i in range(num_sections):
        ua, ub, xa, xb, n = pdict["ps_params"].get_index(i)
        pos_i_ = pos_[xa:xb]
        quat_i_ = quat_[xa:xb]
        u_i_ = u_[ua:ub]

        # rate constraint

        att = pdict["params"][i]["attitude"]

        # attitude hold : angular velocity is zero
        if att in ["hold", "vertical"]:
            con.append(u_i_)

        # kick-turn : pitch rate constant, roll/yaw rate is zero
        elif att == "kick-turn" or att == "pitch":
            con.append(u_i_[:, 0])
            con.append(u_i_[:, 2])
            con.append(u_i_[1:, 1] - u_i_[0, 1])

        # pitch-yaw : pitch/yaw constant, roll ANGLE is zero
        elif att == "pitch-yaw":
            con.append(u_i_[1:, 1] - u_i_[0, 1])
            con.append(u_i_[1:, 2] - u_i_[0, 2])
            con.append(roll_direction_array(pos_i_[1:] * unit_pos, quat_i_[1:]))

        # same-rate : pitch/yaw is the same as previous section, roll ANGLE is zero
        elif att == "same-rate":
            uf_prev = u_[ua - 1]
            con.append(u_i_[:, 1] - uf_prev[1])
            con.append(u_i_[:, 2] - uf_prev[2])
            con.append(roll_direction_array(pos_i_[1:] * unit_pos, quat_i_[1:]))

        # zero-lift-turn or free : roll hold
        elif att == "zero-lift-turn" or att == "free":
            con.append(u_i_[:, 0])

        else:
            print("ERROR: UNKNOWN ATTITUDE OPTION! ({})".format(att))
            sys.exit()

    return np.concatenate(con, axis=None)


def equality_length_6DoF_rate(xdict, pdict, unitdict, condition):
    """length of equality_6DoF_rate"""

    res = 0
    num_sections = pdict["num_sections"]

    for i in range(num_sections):
        ua, ub, xa, xb, n = pdict["ps_params"].get_index(i)
        # rate constraint

        att = pdict["params"][i]["attitude"]

        # attitude hold : angular velocity is zero
        if att in ["hold", "vertical"]:
            res += 3 * n

        # kick-turn : pitch rate constant, roll/yaw rate is zero
        elif att == "kick-turn" or att == "pitch":
            res += 3 * n - 1

        # pitch-yaw : pitch/yaw constant, roll ANGLE is zero
        elif att == "pitch-yaw":
            res += 3 * n - 2

        # same-rate : pitch/yaw is the same as previous section, roll ANGLE is zero
        elif att == "same-rate":
            res += 3 * n

        # zero-lift-turn or free : roll hold
        elif att == "zero-lift-turn" or att == "free":
            res += n

        else:
            print("ERROR: UNKNOWN ATTITUDE OPTION! ({})".format(att))
            sys.exit()

    return res


def roll_direction_array_gradient(pos, quat, unit_pos, dx):
    """
    Calculate gradient of roll direction vector array.
    
    Computes the partial derivatives of roll direction vectors with respect to 
    position and quaternion for an array of trajectory points using finite difference.
    
    Args:
        pos (numpy.ndarray): Normalized position array in ECI frame (n, 3)
        quat (numpy.ndarray): Attitude quaternion array (n, 4)
        unit_pos (float): Position unit scaling factor [m]
        dx (float): Finite difference step size
    
    Returns:
        dict: Dictionary containing gradient arrays:
            - 'position': ∂roll_dir/∂pos (n, 3)
            - 'quaternion': ∂roll_dir/∂quat (n, 4)
    """
    grad = {
        "position": np.zeros((len(pos), 3)),
        "quaternion": np.zeros((len(pos), 4)),
    }

    f_c = roll_direction_array(pos * unit_pos, quat)

    for j in range(3):
        pos[:, j] += dx
        f_p = roll_direction_array(pos * unit_pos, quat)
        pos[:, j] -= dx
        grad["position"][:, j] = (f_p - f_c) / dx

    for j in range(4):
        quat[:, j] += dx
        f_p = roll_direction_array(pos * unit_pos, quat)
        quat[:, j] -= dx
        grad["quaternion"][:, j] = (f_p - f_c) / dx

    return grad


def equality_jac_6DoF_rate(xdict, pdict, unitdict, condition):
    """Compute Jacobian of 6-DOF rate equality constraints.
    
    Calculates partial derivatives of attitude rate continuity constraints with respect to 
    position, quaternion, and control inputs (angular rates).
    
    Args:
        xdict (dict): Dictionary containing state variables ('position', 'quaternion', 'u')
        pdict (dict): Dictionary with problem parameters and 'dx' (finite difference step)
        unitdict (dict): Dictionary of unit scaling factors
        condition (dict): Configuration dictionary with attitude control specifications
    
    Returns:
        dict: Jacobian matrices in COO sparse format for 'position', 'quaternion', and 'u'
    """

    jac = {}
    dx = pdict["dx"]

    unit_pos = unitdict["position"]

    pos_ = xdict["position"].reshape(-1, 3)
    quat_ = xdict["quaternion"].reshape(-1, 4)

    num_sections = pdict["num_sections"]

    nRow = equality_length_6DoF_rate(xdict, pdict, unitdict, condition)
    jac["position"] = {"coo": [[], [], []], "shape": (nRow, pdict["M"] * 3)}
    jac["quaternion"] = {"coo": [[], [], []], "shape": (nRow, pdict["M"] * 4)}
    jac["u"] = {"coo": [[], [], []], "shape": (nRow, pdict["N"] * 3)}

    iRow = 0

    for i in range(num_sections):
        ua, ub, xa, xb, n = pdict["ps_params"].get_index(i)
        pos_i_ = pos_[xa:xb]
        quat_i_ = quat_[xa:xb]

        # rate constraint

        att = pdict["params"][i]["attitude"]
        # attitude hold : angular velocity is zero
        if att in ["hold", "vertical"]:
            jac["u"]["coo"][0].extend(list(range(iRow, iRow + n * 3)))
            jac["u"]["coo"][1].extend(list(range(ua * 3, (ua + n) * 3)))
            jac["u"]["coo"][2].extend([1.0] * (n * 3))
            iRow += n * 3

        # kick-turn : pitch rate constant, roll/yaw rate is zero
        elif att == "kick-turn" or att == "pitch":
            jac["u"]["coo"][0].extend(list(range(iRow, iRow + n)))
            jac["u"]["coo"][1].extend(list(range(ua * 3, (ua + n) * 3, 3)))
            jac["u"]["coo"][2].extend([1.0] * n)
            iRow += n
            jac["u"]["coo"][0].extend(list(range(iRow, iRow + n)))
            jac["u"]["coo"][1].extend(list(range(ua * 3 + 2, (ua + n) * 3 + 2, 3)))
            jac["u"]["coo"][2].extend([1.0] * n)
            iRow += n
            jac["u"]["coo"][0].extend(list(range(iRow, iRow + n - 1)))
            jac["u"]["coo"][1].extend([ua * 3 + 1] * (n - 1))
            jac["u"]["coo"][2].extend([-1.0] * (n - 1))
            jac["u"]["coo"][0].extend(list(range(iRow, iRow + n - 1)))
            jac["u"]["coo"][1].extend(
                list(range((ua + 1) * 3 + 1, (ua + n) * 3 + 1, 3))
            )
            jac["u"]["coo"][2].extend([1.0] * (n - 1))
            iRow += n - 1

        # pitch-yaw : pitch/yaw constant, roll ANGLE is zero
        elif att == "pitch-yaw":
            jac["u"]["coo"][0].extend(list(range(iRow, iRow + n - 1)))
            jac["u"]["coo"][1].extend([ua * 3 + 1] * (n - 1))
            jac["u"]["coo"][2].extend([-1.0] * (n - 1))
            jac["u"]["coo"][0].extend(list(range(iRow, iRow + n - 1)))
            jac["u"]["coo"][1].extend(
                list(range((ua + 1) * 3 + 1, (ua + n) * 3 + 1, 3))
            )
            jac["u"]["coo"][2].extend([1.0] * (n - 1))
            iRow += n - 1
            jac["u"]["coo"][0].extend(list(range(iRow, iRow + n - 1)))
            jac["u"]["coo"][1].extend([ua * 3 + 2] * (n - 1))
            jac["u"]["coo"][2].extend([-1.0] * (n - 1))
            jac["u"]["coo"][0].extend(list(range(iRow, iRow + n - 1)))
            jac["u"]["coo"][1].extend(
                list(range((ua + 1) * 3 + 2, (ua + n) * 3 + 2, 3))
            )
            jac["u"]["coo"][2].extend([1.0] * (n - 1))
            iRow += n - 1

            dfdx = roll_direction_array_gradient(pos_i_[1:], quat_i_[1:], unit_pos, dx)
            for j in range(3):
                jac["position"]["coo"][0].extend(list(range(iRow, iRow + n)))
                jac["position"]["coo"][1].extend(
                    list(range((xa + 1) * 3 + j, (xa + n + 1) * 3 + j, 3))
                )
                jac["position"]["coo"][2].extend(dfdx["position"][:, j])
            for j in range(4):
                jac["quaternion"]["coo"][0].extend(list(range(iRow, iRow + n)))
                jac["quaternion"]["coo"][1].extend(
                    list(range((xa + 1) * 4 + j, (xa + n + 1) * 4 + j, 4))
                )
                jac["quaternion"]["coo"][2].extend(dfdx["quaternion"][:, j])
            iRow += n

        # same-rate : pitch/yaw is the same as previous section, roll ANGLE is zero
        elif att == "same-rate":
            jac["u"]["coo"][0].extend(list(range(iRow, iRow + n)))
            jac["u"]["coo"][1].extend([ua * 3 - 2] * n)
            jac["u"]["coo"][2].extend([-1.0] * n)
            jac["u"]["coo"][0].extend(list(range(iRow, iRow + n)))
            jac["u"]["coo"][1].extend(list(range(ua * 3 + 1, (ua + n) * 3 + 1, 3)))
            jac["u"]["coo"][2].extend([1.0] * n)
            iRow += n
            jac["u"]["coo"][0].extend(list(range(iRow, iRow + n)))
            jac["u"]["coo"][1].extend([ua * 3 - 1] * n)
            jac["u"]["coo"][2].extend([-1.0] * n)
            jac["u"]["coo"][0].extend(list(range(iRow, iRow + n)))
            jac["u"]["coo"][1].extend(list(range(ua * 3 + 2, (ua + n) * 3 + 2, 3)))
            jac["u"]["coo"][2].extend([1.0] * n)
            iRow += n

            dfdx = roll_direction_array_gradient(pos_i_[1:], quat_i_[1:], unit_pos, dx)
            for j in range(3):
                jac["position"]["coo"][0].extend(list(range(iRow, iRow + n)))
                jac["position"]["coo"][1].extend(
                    list(range((xa + 1) * 3 + j, (xa + n + 1) * 3 + j, 3))
                )
                jac["position"]["coo"][2].extend(dfdx["position"][:, j])
            for j in range(4):
                jac["quaternion"]["coo"][0].extend(list(range(iRow, iRow + n)))
                jac["quaternion"]["coo"][1].extend(
                    list(range((xa + 1) * 4 + j, (xa + n + 1) * 4 + j, 4))
                )
                jac["quaternion"]["coo"][2].extend(dfdx["quaternion"][:, j])
            iRow += n

        # zero-lift-turn or free : roll hold
        elif att == "zero-lift-turn" or att == "free":
            jac["u"]["coo"][0].extend(list(range(iRow, iRow + n)))
            jac["u"]["coo"][1].extend(list(range(ua * 3, (ua + n) * 3, 3)))
            jac["u"]["coo"][2].extend([1.0] * n)
            iRow += n

        else:
            print("ERROR: UNKNOWN ATTITUDE OPTION! ({})".format(att))
            sys.exit()

    for key in jac.keys():
        jac[key]["coo"][0] = np.array(jac[key]["coo"][0], dtype="i4")
        jac[key]["coo"][1] = np.array(jac[key]["coo"][1], dtype="i4")
        jac[key]["coo"][2] = np.array(jac[key]["coo"][2], dtype="f8")

    return jac
