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
Initial, Terminal, and Knot Constraints Module
===============================================

Module for defining boundary conditions and constraints at knot points for optimization problems.

This module defines initial conditions, terminal conditions, and continuity/connection
conditions at knot points between phases for trajectory optimization problems.

Main Features:
    * Equality constraints for initial conditions (position, velocity, attitude, mass)
    * Equality constraints for terminal conditions (orbital elements, velocity, altitude, etc.)
    * Continuity constraints for state variables at knot points
    * Mass change constraints during stage separation
    * Constraints on orbital elements (semi-major axis, eccentricity, inclination, etc.)

Constraint Functions:
    equality_initial: Equality constraints for initial conditions
    equality_terminal: Equality constraints for terminal conditions
    equality_knot: Continuity constraints at knot points
"""

import numpy as np
from math import cos, radians
from .coordinate_c import (
    angular_momentum,
    orbit_energy,
    inclination_rad,
    angular_momentum_from_altitude,
    orbit_energy_from_altitude,
)


# constraints_a.py
# constraints about initial, knotting and terminal conditions


def equality_init(xdict, pdict, unitdict, condition):
    """Equality constraint for initial conditions at trajectory start.
    
    Enforces that the initial state (mass, position, velocity, quaternion) matches 
    the specified initial conditions from the configuration.
    
    Args:
        xdict (dict): Dictionary containing state variables:
            - 'mass': Array of mass values
            - 'position': Flattened array of position vectors (reshape to (N,3))
            - 'velocity': Flattened array of velocity vectors (reshape to (N,3))
            - 'quaternion': Flattened array of quaternions (reshape to (N,4))
        pdict (dict): Dictionary containing optimization parameters
        unitdict (dict): Dictionary of unit scaling factors for normalization
        condition (dict): Configuration dictionary with 'init' and 'OptimizationMode'
    
    Returns:
        numpy.ndarray: Concatenated array of constraint violations (should be zero at optimum).
            Length depends on OptimizationMode (10 or 11 constraints).
    """

    con = []
    mass_ = xdict["mass"]
    pos_ = xdict["position"].reshape(-1, 3)
    vel_ = xdict["velocity"].reshape(-1, 3)
    quat_ = xdict["quaternion"].reshape(-1, 4)

    # initial condition
    if condition["OptimizationMode"] != "Payload":
        con.append(mass_[0] - condition["init"]["mass"] / unitdict["mass"])
    con.append(pos_[0] - condition["init"]["position"] / unitdict["position"])
    con.append(vel_[0] - condition["init"]["velocity"] / unitdict["velocity"])
    con.append(quat_[0] - condition["init"]["quaternion"])

    return np.concatenate(con, axis=None)


def equality_jac_init(xdict, pdict, unitdict, condition):
    """Compute Jacobian matrix of equality_init constraint.
    
    Calculates the partial derivatives of initial condition constraints with respect 
    to all state variables. Returns sparse Jacobian in COO format.
    
    Args:
        xdict (dict): Dictionary containing state variables
        pdict (dict): Dictionary containing optimization parameters with 'M' (total nodes)
        unitdict (dict): Dictionary of unit scaling factors
        condition (dict): Configuration dictionary with 'OptimizationMode'
    
    Returns:
        dict: Jacobian matrices in COO sparse format for each state variable:
            - Keys: 'position', 'velocity', 'quaternion', and optionally 'mass'
            - Values: Dict with 'coo' (row, col, data) and 'shape' (nRow, nCol)
    """

    jac = {}

    if condition["OptimizationMode"] == "Payload":
        jac["position"] = {
            "coo": [
                np.arange(0, 3, dtype="i4"),
                np.arange(0, 3, dtype="i4"),
                np.ones(3),
            ],
            "shape": (10, pdict["M"] * 3),
        }
        jac["velocity"] = {
            "coo": [
                np.arange(3, 6, dtype="i4"),
                np.arange(0, 3, dtype="i4"),
                np.ones(3),
            ],
            "shape": (10, pdict["M"] * 3),
        }
        jac["quaternion"] = {
            "coo": [
                np.arange(6, 10, dtype="i4"),
                np.arange(0, 4, dtype="i4"),
                np.ones(4),
            ],
            "shape": (10, pdict["M"] * 4),
        }

    else:
        jac["mass"] = {
            "coo": [np.zeros(1, dtype="i4"), np.zeros(1, dtype="i4"), np.ones(1)],
            "shape": (11, pdict["M"]),
        }
        jac["position"] = {
            "coo": [
                np.arange(1, 4, dtype="i4"),
                np.arange(0, 3, dtype="i4"),
                np.ones(3),
            ],
            "shape": (11, pdict["M"] * 3),
        }
        jac["velocity"] = {
            "coo": [
                np.arange(4, 7, dtype="i4"),
                np.arange(0, 3, dtype="i4"),
                np.ones(3),
            ],
            "shape": (11, pdict["M"] * 3),
        }
        jac["quaternion"] = {
            "coo": [
                np.arange(7, 11, dtype="i4"),
                np.arange(0, 4, dtype="i4"),
                np.ones(4),
            ],
            "shape": (11, pdict["M"] * 4),
        }

    return jac


def equality_time(xdict, pdict, unitdict, condition):
    """Equality constraint fixing event times at specified knot points.
    
    Enforces that event times match predefined reference values when specified
    in the configuration. Forces initial time and relative timing between events.
    
    Args:
        xdict (dict): Dictionary containing state variables including 't' (event times)
        pdict (dict): Dictionary with 'params' (section parameters) and 'event_index'
        unitdict (dict): Dictionary of unit scaling factors including 't' (time unit)
        condition (dict): Configuration dictionary
    
    Returns:
        numpy.ndarray: Constraint values (equals 0 when satisfied)
    """

    con = []
    unit_t = unitdict["t"]

    t_ = xdict["t"]

    num_sections = pdict["num_sections"]

    # force to fix initial time
    con.append(t_[0] - pdict["params"][0]["time"] / unit_t)
    for i in range(1, num_sections + 1):
        if pdict["params"][i]["time_ref"] in pdict["event_index"].keys():
            i_ref = pdict["event_index"][pdict["params"][i]["time_ref"]]
            con.append(
                t_[i]
                - t_[i_ref]
                - (pdict["params"][i]["time"] - pdict["params"][i_ref]["time"]) / unit_t
            )

    return np.concatenate(con, axis=None)


def equality_jac_time(xdict, pdict, unitdict, condition):
    """Compute Jacobian matrix of equality_time constraint.
    
    Calculates the partial derivatives of time constraints with respect to 
    event times (knot points).
    
    Args:
        xdict (dict): Dictionary containing state variables including 't' (event times)
        pdict (dict): Dictionary containing parameters with 'event_index' and 'num_sections'
        unitdict (dict): Dictionary of unit scaling factors
        condition (dict): Configuration dictionary
    
    Returns:
        dict: Jacobian matrix in COO sparse format with key 't', containing
            'coo' (row, col, data) and 'shape' (nRow, nCol).
    """

    jac = {}

    data = [1.0]
    row = [0]
    col = [0]

    iRow = 1
    for i in range(1, pdict["num_sections"] + 1):
        if pdict["params"][i]["time_ref"] in pdict["event_index"].keys():
            i_ref = pdict["event_index"][pdict["params"][i]["time_ref"]]
            data.extend([1.0, -1.0])
            row.extend([iRow, iRow])
            col.extend([i, i_ref])
            iRow += 1

    jac["t"] = {
        "coo": [np.array(row, dtype="i4"), np.array(col, dtype="i4"), np.array(data)],
        "shape": (iRow, len(xdict["t"])),
    }

    return jac


def equality_knot_LGR(xdict, pdict, unitdict, condition):
    """Equality constraint for state continuity at knot points (LGR method).
    
    Enforces continuity of state variables (mass, position, velocity, quaternion) 
    at the boundaries between trajectory sections for Legendre-Gauss-Radau (LGR) 
    pseudospectral method. Handles mass discontinuities at stage separations.
    
    Args:
        xdict (dict): Dictionary containing state variables:
            - 'mass': Array of mass values
            - 'position', 'velocity': Flattened arrays reshaped to (N,3)
            - 'quaternion': Flattened array reshaped to (N,4)
        pdict (dict): Dictionary with 'ps_params' for section indexing, 'RocketStage' for staging info
        unitdict (dict): Dictionary of unit scaling factors
        condition (dict): Configuration dictionary
    
    Returns:
        numpy.ndarray or None: Concatenated constraint violations at all knot points,
            or None if there are no knot points to constrain.
    """

    con = []

    mass_ = xdict["mass"]
    pos_ = xdict["position"].reshape(-1, 3)
    vel_ = xdict["velocity"].reshape(-1, 3)
    quat_ = xdict["quaternion"].reshape(-1, 4)

    num_sections = pdict["num_sections"]

    param = np.zeros(5)

    section_sep_list = []
    for key, stage in pdict["RocketStage"].items():
        if stage["separation_at"] is not None:
            section_ig = [
                i
                for i, value in enumerate(pdict["params"])
                if value["name"] == stage["ignition_at"]
            ][0]
            section_sep = [
                i
                for i, value in enumerate(pdict["params"])
                if value["name"] == stage["separation_at"]
            ][0]
            section_sep_list.append(section_sep)

            # mass after separation
            mass_stage = (
                stage["mass_dry"]
                + stage["mass_propellant"]
                + sum([item["mass"] for item in stage["dropMass"].values()])
            )
            index_ig = pdict["ps_params"].index_start_x(section_ig)
            index_sep = pdict["ps_params"].index_start_x(section_sep)
            con.append(
                mass_[index_ig] - mass_[index_sep] - mass_stage / unitdict["mass"]
            )

    for i in range(1, num_sections):
        xa = pdict["ps_params"].index_start_x(i)

        param[0] = pdict["params"][i]["thrust"]
        param[1] = pdict["params"][i]["massflow"]
        param[2] = pdict["params"][i]["reference_area"]
        param[4] = pdict["params"][i]["nozzle_area"]

        # knotting constraints
        mass_init_ = mass_[xa]
        mass_prev_ = mass_[xa - 1]
        if not (i in section_sep_list):
            con.append(
                mass_init_
                - mass_prev_
                + pdict["params"][i]["mass_jettison"] / unitdict["mass"]
            )

        pos_init_ = pos_[xa]
        pos_prev_ = pos_[xa - 1]
        con.append(pos_init_ - pos_prev_)

        vel_init_ = vel_[xa]
        vel_prev_ = vel_[xa - 1]
        con.append(vel_init_ - vel_prev_)

        quat_init_ = quat_[xa]
        quat_prev_ = quat_[xa - 1]
        con.append(quat_init_ - quat_prev_)

    return np.concatenate(con, axis=None)


def equality_jac_knot_LGR(xdict, pdict, unitdict, condition):
    """Compute Jacobian matrix of equality_knot_LGR constraint.
    
    Calculates the partial derivatives of knot point continuity constraints with 
    respect to state variables at section boundaries. Returns sparse Jacobian in COO format.
    
    Args:
        xdict (dict): Dictionary containing state variables
        pdict (dict): Dictionary with section parameters, 'ps_params', and 'RocketStage'
        unitdict (dict): Dictionary of unit scaling factors
        condition (dict): Configuration dictionary
    
    Returns:
        dict or None: Jacobian matrices in COO sparse format for 'mass', 'position',
            'velocity', 'quaternion', or None if no knot constraints exist.
    """

    jac = {}

    num_sections = pdict["num_sections"]

    f_center = equality_knot_LGR(xdict, pdict, unitdict, condition)
    nRow = len(f_center)

    jac["mass"] = {"coo": [[], [], []], "shape": (nRow, pdict["M"])}
    jac["position"] = {"coo": [[], [], []], "shape": (nRow, pdict["M"] * 3)}
    jac["velocity"] = {"coo": [[], [], []], "shape": (nRow, pdict["M"] * 3)}
    jac["quaternion"] = {"coo": [[], [], []], "shape": (nRow, pdict["M"] * 4)}

    iRow = 0

    section_sep_list = []
    for key, stage in pdict["RocketStage"].items():
        if stage["separation_at"] is not None:
            section_ig = [
                i
                for i, value in enumerate(pdict["params"])
                if value["name"] == stage["ignition_at"]
            ][0]
            section_sep = [
                i
                for i, value in enumerate(pdict["params"])
                if value["name"] == stage["separation_at"]
            ][0]
            section_sep_list.append(section_sep)

            # mass after separation
            index_ig = pdict["ps_params"].index_start_x(section_ig)
            index_sep = pdict["ps_params"].index_start_x(section_sep)
            jac["mass"]["coo"][0].extend([iRow, iRow])
            jac["mass"]["coo"][1].extend([index_ig, index_sep])
            jac["mass"]["coo"][2].extend([1.0, -1.0])
            iRow += 1

    for i in range(1, num_sections):
        xa = pdict["ps_params"].index_start_x(i)

        if not (i in section_sep_list):
            jac["mass"]["coo"][0].extend([iRow, iRow])
            jac["mass"]["coo"][1].extend([xa - 1, xa])
            jac["mass"]["coo"][2].extend([-1.0, 1.0])
            iRow += 1

        jac["position"]["coo"][0].extend(list(range(iRow, iRow + 3)))
        jac["position"]["coo"][1].extend(list(range((xa - 1) * 3, (xa) * 3)))
        jac["position"]["coo"][2].extend([-1.0] * 3)
        jac["position"]["coo"][0].extend(list(range(iRow, iRow + 3)))
        jac["position"]["coo"][1].extend(list(range((xa) * 3, (xa + 1) * 3)))
        jac["position"]["coo"][2].extend([1.0] * 3)
        iRow += 3

        jac["velocity"]["coo"][0].extend(list(range(iRow, iRow + 3)))
        jac["velocity"]["coo"][1].extend(list(range((xa - 1) * 3, (xa) * 3)))
        jac["velocity"]["coo"][2].extend([-1.0] * 3)
        jac["velocity"]["coo"][0].extend(list(range(iRow, iRow + 3)))
        jac["velocity"]["coo"][1].extend(list(range((xa) * 3, (xa + 1) * 3)))
        jac["velocity"]["coo"][2].extend([1.0] * 3)
        iRow += 3

        jac["quaternion"]["coo"][0].extend(list(range(iRow, iRow + 4)))
        jac["quaternion"]["coo"][1].extend(list(range((xa - 1) * 4, (xa) * 4)))
        jac["quaternion"]["coo"][2].extend([-1.0] * 4)
        jac["quaternion"]["coo"][0].extend(list(range(iRow, iRow + 4)))
        jac["quaternion"]["coo"][1].extend(list(range((xa) * 4, (xa + 1) * 4)))
        jac["quaternion"]["coo"][2].extend([1.0] * 4)
        iRow += 4

    for key in jac.keys():
        jac[key]["coo"][0] = np.array(jac[key]["coo"][0], dtype="i4")
        jac[key]["coo"][1] = np.array(jac[key]["coo"][1], dtype="i4")
        jac[key]["coo"][2] = np.array(jac[key]["coo"][2], dtype="f8")

    return jac


def equality_6DoF_LGR_terminal(xdict, pdict, unitdict, condition):
    """Equality constraint for terminal (final) orbital conditions.
    
    Enforces terminal constraints such as target perigee/apogee altitude, inclination,
    right ascension of ascending node (RAAN), and argument of perigee based on the
    final state of the trajectory.
    
    Args:
        xdict (dict): Dictionary containing state variables with final position and velocity
        pdict (dict): Dictionary containing optimization parameters
        unitdict (dict): Dictionary of unit scaling factors for position and velocity
        condition (dict): Configuration dictionary with terminal conditions:
            - 'altitude_perigee', 'altitude_apogee': Target orbit altitudes [m]
            - 'inclination', 'arg_perigee', 'RAAN': Orbital elements [degrees or rad]
    
    Returns:
        numpy.ndarray or None: Array of terminal constraint violations, or None if 
            no terminal constraints are specified.
    """

    con = []

    unit_pos = unitdict["position"]
    unit_vel = unitdict["velocity"]

    # terminal conditions

    pos_f = xdict["position"][-3:] * unit_pos
    vel_f = xdict["velocity"][-3:] * unit_vel

    GMe = 3.986004418e14
    if (
        condition["altitude_perigee"] is not None
        and condition["altitude_apogee"] is not None
    ):
        c_target = angular_momentum_from_altitude(
            condition["altitude_perigee"],
            condition["altitude_apogee"],
        )
        e_target = orbit_energy_from_altitude(
            condition["altitude_perigee"],
            condition["altitude_apogee"],
        )
    else:
        c_target = condition["radius"] * condition["vel_tangential_geocentric"]
        vf_target = condition["vel_tangential_geocentric"] / cos(
            radians(condition["flightpath_vel_inertial_geocentric"])
        )
        e_target = vf_target**2 / 2.0 - GMe / condition["radius"]

    c = angular_momentum(pos_f, vel_f)
    e = orbit_energy(pos_f, vel_f)
    con.append((e / e_target) - 1.0)  # orbit energy
    con.append((c / c_target) - 1.0)  # angular momentum

    if condition["inclination"] is not None:
        inc = inclination_rad(pos_f, vel_f)
        inc_target = radians(condition["inclination"])
        con.append(inc - inc_target)

    return np.concatenate(con, axis=None)


def equality_jac_6DoF_LGR_terminal(xdict, pdict, unitdict, condition):
    """Compute Jacobian matrix of equality_6DoF_LGR_terminal constraint.
    
    Calculates the partial derivatives of terminal orbital condition constraints 
    with respect to final position and velocity using finite differences.
    
    Args:
        xdict (dict): Dictionary containing state variables
        pdict (dict): Dictionary with 'dx' (finite difference step size) and 'M' (total nodes)
        unitdict (dict): Dictionary of unit scaling factors
        condition (dict): Configuration with terminal orbital constraints
    
    Returns:
        dict or None: Jacobian matrices in COO sparse format for 'position' and 
            'velocity', or None if no terminal constraints exist.
    """

    jac = {}
    dx = pdict["dx"]

    f_center = equality_6DoF_LGR_terminal(xdict, pdict, unitdict, condition)

    nRow = 2
    if condition["inclination"] is not None:
        nRow += 1

    nCol = pdict["M"] * 3
    jac["position"] = {"coo": [[], [], []], "shape": (nRow, nCol)}
    jac["velocity"] = {"coo": [[], [], []], "shape": (nRow, nCol)}

    for key in ["position", "velocity"]:

        for j in range(nCol - 3, nCol):
            xdict[key][j] += dx
            f_p = equality_6DoF_LGR_terminal(xdict, pdict, unitdict, condition)
            xdict[key][j] -= dx
            jac[key]["coo"][0].extend(list(range(nRow)))
            jac[key]["coo"][1].extend([j] * nRow)
            jac[key]["coo"][2].extend(((f_p - f_center) / dx).tolist())

        jac[key]["coo"][0] = np.array(jac[key]["coo"][0], dtype="i4")
        jac[key]["coo"][1] = np.array(jac[key]["coo"][1], dtype="i4")
        jac[key]["coo"][2] = np.array(jac[key]["coo"][2], dtype="f8")

    return jac


def inequality_time(xdict, pdict, unitdict, condition):
    """Inequality constraint ensuring forward time progression at knot points.
    
    Enforces that event times are monotonically increasing (later events occur after earlier ones)
    for sections where time is not explicitly fixed by the configuration.
    
    Args:
        xdict (dict): Dictionary containing state variables including 't' (event times)
        pdict (dict): Dictionary with 'params' (section parameters) and 'event_index'
        unitdict (dict): Dictionary of unit scaling factors
        condition (dict): Configuration dictionary
    
    Returns:
        numpy.ndarray: Constraint values (>= 0 when satisfied, ensuring t[i+1] >= t[i])
    """

    con = []
    t_normal = xdict["t"]

    for i in range(pdict["num_sections"]):
        if not (
            pdict["params"][i]["time_ref"] in pdict["event_index"].keys()
            and pdict["params"][i + 1]["time_ref"] in pdict["event_index"].keys()
        ):
            con.append(t_normal[i + 1] - t_normal[i])

    return np.array(con)


def inequality_jac_time(xdict, pdict, unitdict, condition):
    """Jacobian of the time ordering inequality constraint.
    
    Computes the derivatives of inequality_time with respect to event times.
    
    Args:
        xdict (dict): Dictionary containing state variables including 't' (event times)
        pdict (dict): Dictionary with 'params' (section parameters) and 'event_index'
        unitdict (dict): Dictionary of unit scaling factors
        condition (dict): Configuration dictionary
    
    Returns:
        dict: Jacobian in COO sparse format with keys 'data', 'row', 'col'
    """

    jac = {}

    data = []
    row = []
    col = []

    counter = 0
    for i in range(pdict["num_sections"]):
        if not (
            pdict["params"][i]["time_ref"] in pdict["event_index"].keys()
            and pdict["params"][i + 1]["time_ref"] in pdict["event_index"].keys()
        ):
            data.extend([-1.0, 1.0])
            row.extend([counter, counter])
            col.extend([i, i + 1])
            counter += 1

    jac["t"] = {
        "coo": [
            np.array(row, dtype="i4"),
            np.array(col, dtype="i4"),
            np.array(data, dtype="f8"),
        ],
        "shape": (counter, len(xdict["t"])),
    }
    return jac
