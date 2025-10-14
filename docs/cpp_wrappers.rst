Python Binding Wrappers
========================

These modules provide wrapper functions to facilitate Python bindings via PyBind11.
The wrappers convert between C++ types (Eigen vectors) and Python-friendly types 
(scalars and arrays) to make the C++ functionality easily accessible from Python.

Atmospheric Property Wrappers
------------------------------

Wrapper functions for atmospheric calculations exposed to Python.

.. doxygenfile:: wrapper_air.hpp
   :project: GELATO_CPP

Coordinate Transformation Wrappers
-----------------------------------

Wrapper functions for coordinate transformations, quaternion operations, and orbital mechanics.

.. doxygenfile:: wrapper_coordinate.hpp
   :project: GELATO_CPP
   :sections: briefdescription detaileddescription func

Dynamics Calculation Wrappers
------------------------------

Wrapper functions for vehicle dynamics calculations including velocity and quaternion dynamics.

.. doxygenfile:: wrapper_dynamics.hpp
   :project: GELATO_CPP

Utility Function Wrappers
--------------------------

Wrapper functions for utility calculations including distance and interpolation.

.. doxygenfile:: wrapper_utils.hpp
   :project: GELATO_CPP
   :sections: briefdescription detaileddescription func
