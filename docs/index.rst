GELATO Documentation
====================

GELATO (Generic Launch Trajectory Optimizer) is an open source tool for launch trajectory optimization, written in Python. It solves trajectory optimization problems using the Legendre-Gauss-Radau pseudospectral method, which is stable and robust for highly complex nonlinear problems.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   tutorial
   usage
   examples
   modules
   api
   cpp_api

Key Features
------------

* **Pseudospectral Method Optimization**: High-precision trajectory optimization using Legendre-Gauss-Radau pseudospectral method
* **Flexible Constraint Settings**: Support for user-defined constraint conditions
* **Atmospheric Model**: Support for US Standard Atmosphere 1976
* **High-Speed Computation**: Acceleration through C++ bindings (PyBind11)
* **KML Output**: Generate KML files for visualization in Google Earth

Indices and Tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
