C++ API Reference
==================

This section documents the C++ implementation and Python bindings that provide 
high-performance computation for GELATO's trajectory optimization.

.. toctree::
   :maxdepth: 2
   
   cpp_core_api
   cpp_wrappers
   pybind_modules

Overview
--------

GELATO uses C++ for performance-critical calculations and exposes them to Python 
through PyBind11 bindings. The documentation is organized into three sections:

**C++ Core API Reference**
   Core C++ classes and functions implementing coordinate transformations, 
   atmospheric models, gravity calculations, and IIP computations.

**Python Binding Wrappers**
   Wrapper functions that convert between C++ types (Eigen vectors) and 
   Python-friendly types to facilitate PyBind11 bindings.

**PyBind11 Module Definitions**
   Documentation of the Python modules (``*_c``) that are automatically 
   generated from C++ code: ``coordinate_c``, ``dynamics_c``, ``IIP_c``, 
   ``USStandardAtmosphere_c``, and ``utils_c``.
