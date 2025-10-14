C++ Core API Reference
=======================

This section documents the core C++ classes and modules that implement the performance-critical components of GELATO.

Coordinate Transformations
---------------------------

The Coordinate module provides functions for converting between different coordinate systems
(ECEF, ECI, NED) and handling quaternion operations.

.. doxygenclass:: Coordinate
   :members:
   :undoc-members:

Atmospheric Model
-----------------

The Air module implements the U.S. Standard Atmosphere model for calculating atmospheric properties
at different altitudes.

.. doxygenfile:: Air.hpp
   :sections: briefdescription detaileddescription

.. doxygenclass:: Air
   :members:
   :undoc-members:

Earth Model
-----------

The Earth module provides Earth-related constants and models (e.g., WGS84 parameters).

.. doxygenclass:: Earth
   :members:
   :undoc-members:

Gravity Calculations
---------------------

The gravity module implements gravity model calculations.

.. doxygenfile:: gravity.hpp

Instantaneous Impact Point (IIP)
---------------------------------

The IIP module calculates the instantaneous impact point for trajectory analysis.

.. doxygenfile:: iip.hpp
