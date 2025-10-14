PyBind11 Module Definitions
============================

Overview
--------

GELATO uses PyBind11 to expose C++ functionality to Python. The following Python modules 
are automatically generated from C++ code and provide high-performance implementations of 
physics calculations, coordinate transformations, and utilities.

Available Modules
-----------------

coordinate_c
^^^^^^^^^^^^

Coordinate transformations, quaternion operations, and orbital mechanics.

**Key Functions:**

- ``ecef2geodetic``, ``geodetic2ecef`` - ECEF/geodetic conversions
- ``ecef2eci``, ``eci2ecef`` - ECEF/ECI conversions  
- ``quatmult``, ``quatrot``, ``normalize`` - Quaternion operations
- ``quat_*`` - Quaternion transformations between frames
- ``orbital_elements`` - Orbital element calculations
- ``distance_vincenty`` - Geodesic distance calculations

dynamics_c
^^^^^^^^^^

Vehicle dynamics calculations including thrust forces, aerodynamic forces, gravity, and wind effects.

**Key Functions:**

- ``dynamics_velocity`` - Velocity dynamics with aerodynamic forces
  
  Computes acceleration considering thrust, aerodynamic drag, gravity, and wind.
  Includes atmospheric properties calculation and Mach number-dependent drag coefficient.

- ``dynamics_velocity_NoAir`` - Velocity dynamics without aerodynamic forces
  
  Simplified dynamics calculation for vacuum flight or when aerodynamic effects are negligible.
  Computes acceleration from thrust and gravity only.

- ``dynamics_quaternion`` - Quaternion rate of change
  
  Computes the time derivative of attitude quaternion from angular velocity (body rates).

IIP_c
^^^^^

Instantaneous Impact Point (IIP) calculations for trajectory analysis.

**Key Functions:**

- ``posLLH_IIP_FAA`` - Calculate IIP using FAA iterative method

USStandardAtmosphere_c
^^^^^^^^^^^^^^^^^^^^^^^

U.S. Standard Atmosphere 1976 model for atmospheric properties.

**Key Functions:**

- ``airdensity_at`` - Air density at given altitude
- ``airpressure_at`` - Air pressure at given altitude
- ``airtemperature_at`` - Air temperature at given altitude
- ``speed_of_sound`` - Speed of sound at given altitude

utils_c
^^^^^^^

Utility functions for wind models, angle of attack, and aerodynamic calculations.

**Key Functions:**

- ``wind_ned`` - Wind velocity in NED frame
- ``angle_of_attack_all_rad`` - Total angle of attack calculation
- ``dynamic_pressure_pa`` - Dynamic pressure calculation
- ``q_alpha_pa_rad`` - Q-alpha (dynamic pressure × angle of attack)

Implementation Details
----------------------

The modules are defined in ``pybind_*.cpp`` files in the ``src/`` directory:

- ``src/pybind_coordinate.cpp`` → ``coordinate_c``
- ``src/pybind_dynamics.cpp`` → ``dynamics_c``
- ``src/pybind_IIP.cpp`` → ``IIP_c``
- ``src/pybind_USStandardAtmosphere.cpp`` → ``USStandardAtmosphere_c``
- ``src/pybind_utils.cpp`` → ``utils_c``

These modules use the wrapper functions from ``wrapper_*.hpp`` files to provide 
Python-friendly interfaces to the C++ implementations.

Usage Example
-------------

.. code-block:: python

   from lib.coordinate_c import ecef2geodetic, distance_vincenty
   from lib.USStandardAtmosphere_c import airdensity_at
   
   # Convert ECEF coordinates to latitude/longitude/altitude
   lat, lon, alt = ecef2geodetic(x, y, z)
   
   # Calculate distance between two points
   distance = distance_vincenty(lat1, lon1, lat2, lon2)
   
   # Get air density at altitude
   rho = airdensity_at(altitude_m)
