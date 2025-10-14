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
 * @file pybind_dynamics.cpp
 * @brief PyBind11 module definition for dynamics calculations
 * 
 * This file creates the Python module 'dynamics_c' which exposes C++
 * dynamics calculation functions to Python for trajectory optimization.
 * The module includes velocity and acceleration calculations considering:
 * - Thrust forces
 * - Aerodynamic forces
 * - Gravity
 * - Wind effects
 */

#include "wrapper_dynamics.hpp"

PYBIND11_MODULE(dynamics_c, m) {
  m.def("dynamics_velocity", &dynamics_velocity,
        "velocity with aerodynamic forces");
  m.def("dynamics_velocity_NoAir", &dynamics_velocity_NoAir,
        "velocity without aerodynamic forces");
  m.def("dynamics_quaternion", &dynamics_quaternion, "quaternion dynamics");
}
