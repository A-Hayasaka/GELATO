Module Structure
================

GELATO consists of the following modules organized by functionality:

Main Scripts
------------

* **Trajectory_Optimization**: Main optimization solver
* **initialize**: Initial trajectory generation and simulation
* **output_result**: Result output and data processing
* **make_kml**: KML file generation for Google Earth visualization

Library Modules
---------------

Constraint Modules
^^^^^^^^^^^^^^^^^^

Modules implementing optimization constraints:

* **con_aero**: Aerodynamic constraints (dynamic pressure, angle of attack, Q-alpha)
* **con_dynamics**: Dynamics constraints (mass, position, velocity, quaternion)
* **con_init_terminal_knot**: Initial/terminal/knot point constraints
* **con_trajectory**: Trajectory-specific constraints
* **con_user**: User-defined constraints
* **con_waypoint**: Waypoint constraints

Physics and Mathematics Modules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Core physics and mathematical models:

* **dynamics**: Rocket equations of motion (6-DOF dynamics)
* **coordinate**: Coordinate transformations and quaternion operations
* **USStandardAtmosphere**: US Standard Atmosphere 1976 model
* **IIP**: Instantaneous Impact Point (IIP) calculation
* **downrange**: Downrange distance calculation using Vincenty formula

Optimization Modules
^^^^^^^^^^^^^^^^^^^^

Pseudospectral method and optimization tools:

* **PSfunctions**: Pseudospectral method functions (LGL, LG, LGR)
* **SectionParameters**: Multi-phase optimization parameter management
* **cost_gradient**: Cost function and gradient calculation
* **jac_fd**: Jacobian computation using finite differences

Utilities
^^^^^^^^^

Helper functions and utilities:

* **utils**: Aerodynamic calculations, wind models, angle of attack

For detailed API documentation, see :doc:`api`.
