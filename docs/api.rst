API Reference
=============

This section provides detailed API reference for all GELATO modules, their functions, and classes.

Main Scripts
------------

.. toctree::
   :maxdepth: 1

   api/Trajectory_Optimization
   api/initialize
   api/output_result
   api/make_kml

Constraint Modules
------------------

.. toctree::
   :maxdepth: 1

   api/con_aero
   api/con_dynamics
   api/con_init_terminal_knot
   api/con_trajectory
   api/con_user
   api/con_waypoint

Optimization
------------

.. toctree::
   :maxdepth: 1

   api/PSfunctions
   api/SectionParameters
   api/cost_gradient
   api/jac_fd

Legacy Code (Reference Only)
-----------------------------

.. warning::
   **Deprecated Python Implementations**
   
   All modules in this section have been moved to ``lib/_LEGACYCODE/`` and are 
   **no longer used** in the codebase. They have been replaced by optimized C++ 
   extensions for performance:
   
   - ``lib._LEGACYCODE.dynamics`` → Use ``lib.dynamics_c``
   - ``lib._LEGACYCODE.coordinate`` → Use ``lib.coordinate_c``
   - ``lib._LEGACYCODE.USStandardAtmosphere`` → Use ``lib.USStandardAtmosphere_c``
   - ``lib._LEGACYCODE.IIP`` → Use ``lib.IIP_c``
   - ``lib._LEGACYCODE.downrange`` → Use ``lib.coordinate_c`` (distance_vincenty)
   - ``lib._LEGACYCODE.utils`` → Use ``lib.utils_c``
   
   These modules are kept for reference purposes only to understand the original 
   Python implementation.

.. toctree::
   :maxdepth: 1

   api/dynamics
   api/coordinate
   api/USStandardAtmosphere
   api/IIP
   api/downrange
   api/utils
