Examples
========

This section provides detailed examples to help you understand how to use GELATO for trajectory optimization.

Basic Example
-------------

The ``example/`` directory contains a complete working example of a two-stage launch vehicle optimization.

Running the Example
^^^^^^^^^^^^^^^^^^^

To run the basic example:

.. code-block:: bash

   cd example
   python ../Trajectory_Optimization.py example-settings.json

This will optimize a trajectory for a two-stage rocket launching to a 200 km circular orbit.

Example Files Overview
----------------------

The example directory contains the following files:

* ``example-settings.json``: Main settings file defining vehicle and optimization parameters
* ``example-events.csv``: Flight events timeline (liftoff, staging, engine cutoff, etc.)
* ``example-CA.csv``: Axial force coefficient data as function of Mach number
* ``example-wind_average.csv``: Wind profile data
* ``example-trajectory_init.csv``: Initial trajectory guess (optional)
* ``user_constraints.py``: User-defined constraint functions
* ``output/``: Directory for output files

Settings File (example-settings.json)
--------------------------------------

The settings file defines all parameters for the optimization problem.

Vehicle Configuration
^^^^^^^^^^^^^^^^^^^^^

The rocket configuration is defined in the ``RocketStage`` section:

.. code-block:: json

   {
     "RocketStage": {
       "1": {
         "mass_dry": 1361.0,
         "mass_propellant": 21500.0,
         "dropMass": {},
         "Isp_vac": 304.0,
         "reference_area": 2.21,
         "ignition_at": "LIFTOFF",
         "cutoff_at": "MECO",
         "separation_at": "SEP1"
       },
       "2": {
         "mass_dry": 544.0,
         "mass_propellant": 4037.0,
         "dropMass": {
           "fairing": {
             "mass": 140.0,
             "separation_at": "FAIRING"
           }
         },
         "Isp_vac": 317.0,
         "reference_area": 0.0,
         "ignition_at": "SEIG",
         "cutoff_at": "SECO",
         "separation_at": "SEP2"
       }
     },
     "mass_payload": 0.0
   }

**Stage Parameters:**

* ``mass_dry``: Dry mass of the stage (without propellant) [kg]
* ``mass_propellant``: Propellant mass [kg]
* ``dropMass``: Jettisoned components (e.g., fairing) with mass and separation event
* ``Isp_vac``: Specific impulse in vacuum [s]
* ``reference_area``: Reference area for aerodynamic calculations [m²]
* ``ignition_at``, ``cutoff_at``, ``separation_at``: Event names for engine and stage operations

Launch Conditions
^^^^^^^^^^^^^^^^^

Define the launch site and initial flight azimuth:

.. code-block:: json

   {
     "LaunchCondition": {
       "lon": 143.45659,
       "lat": 42.50587,
       "altitude": 50.0,
       "flight_azimuth_init": 90.0
     }
   }

* ``lon``, ``lat``: Launch site coordinates [degrees]
* ``altitude``: Launch site altitude [m]
* ``flight_azimuth_init``: Initial flight azimuth angle [degrees] (0=North, 90=East)

Terminal Conditions
^^^^^^^^^^^^^^^^^^^

Define the target orbit:

.. code-block:: json

   {
     "TerminalCondition": {
       "altitude_perigee": 200000,
       "altitude_apogee": 200000,
       "inclination": null,
       "radius": 6578137,
       "vel_tangential_geocentric": 7784.3,
       "flightpath_vel_inertial_geocentric": 0.0
     }
   }

* ``altitude_perigee``, ``altitude_apogee``: Orbit altitudes [m]
* ``inclination``: Orbital inclination [degrees] (null = determined by launch azimuth)
* ``radius``: Target orbital radius at terminal point [m]
* ``vel_tangential_geocentric``: Tangential velocity component [m/s]
* ``flightpath_vel_inertial_geocentric``: Radial velocity component [m/s] (0 for circular orbit)

Flight Constraints
^^^^^^^^^^^^^^^^^^

Define safety and operational constraints:

.. code-block:: json

   {
     "FlightConstraint": {
       "AOA_max": {
         "MECO": {
           "value": 10.0,
           "range": "initial"
         }
       },
       "dynamic_pressure_max": {},
       "Q_alpha_max": {
         "ZEROLIFT_START": {
           "value": 30000.0,
           "range": "all"
         }
       },
       "waypoint": {
         "FAIRING": {
           "altitude": {
             "exact": 100000.0
           },
           "lon_IIP": {
             "min": 145.0
           }
         }
       },
       "antenna": {
         "ANT1": {
           "lon": 143.45659,
           "lat": 42.50587,
           "altitude": 50.0,
           "elevation_min": {
             "SECO": 0.0
           }
         }
       }
     }
   }

**Constraint Types:**

* ``AOA_max``: Maximum angle of attack [degrees]
* ``dynamic_pressure_max``: Maximum dynamic pressure [Pa]
* ``Q_alpha_max``: Maximum Q-alpha (dynamic pressure × angle of attack) [Pa·deg]
* ``waypoint``: Specific state requirements at certain events
* ``antenna``: Ground station visibility constraints

Optimizer Settings
^^^^^^^^^^^^^^^^^^

Configure the IPOPT solver:

.. code-block:: json

   {
     "IPOPT_": {
       "linear_solver": "mumps",
       "tol": 1.0e-6,
       "acceptable_tol": 1.0e-4,
       "max_iter": 2000
     }
   }

Events File (example-events.csv)
---------------------------------

The events file defines the flight timeline and attitude control segments.

CSV Format
^^^^^^^^^^

.. csv-table:: Events File Structure
   :header: "Column", "Description", "Example"
   :widths: 20, 50, 30

   "name", "Event name (must be unique)", "LIFTOFF"
   "time", "Event time [s] (can be empty for optimization)", "10"
   "time_ref", "Reference event for relative timing", "LIFTOFF"
   "rocketStage", "Active stage number", "1"
   "engineOn", "Engine status (True/False)", "True"
   "thrust", "Thrust magnitude [N]", "420000"
   "nozzle_area", "Nozzle exit area [m²]", "0.68"
   "attitude", "Attitude control mode", "vertical"
   "pitchrate_init", "Initial pitch rate [deg/s]", "-1.0"
   "yawrate_init", "Initial yaw rate [deg/s]", "0.0"
   "num_nodes", "Number of collocation nodes", "5"

Attitude Control Modes
^^^^^^^^^^^^^^^^^^^^^^

* ``vertical``: Maintain vertical attitude
* ``kick-turn``: Initial pitch maneuver with specified pitch rate
* ``zero-lift-turn``: Gravity turn with zero angle of attack
* ``pitch-yaw``: Direct pitch and yaw rate control (optimization variables)
* ``same-rate``: Continue with the same pitch/yaw rates as previous segment
* ``hold``: Maintain current attitude

Example Event Sequence
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   LIFTOFF (t=0s)        → Vertical ascent, Stage 1 ignition
   KICKTURN (t=10s)      → Begin pitch maneuver
   ZEROLIFT_START (t=20s) → Start gravity turn
   ZEROLIFT_END (t=90s)   → End gravity turn, begin rate-controlled ascent
   MECO (t=169s)         → Main Engine Cutoff (optimized)
   SEP1 (t=174s)         → Stage 1 separation
   SEIG (t=179s)         → Stage 2 Engine Ignition
   FAIRING (t=189s)      → Fairing jettison at 100 km altitude
   SECO (t=597s)         → Second stage Engine Cutoff (optimized)
   SEP2 (t=627s)         → Stage 2 separation
   SIMEND (t=630s)       → End of simulation

User Constraints (user_constraints.py)
---------------------------------------

You can define custom constraint functions in ``user_constraints.py``.

Example Template
^^^^^^^^^^^^^^^^

.. code-block:: python

   import numpy as np
   from lib.coordinate import quat2dcm, pos_ecef2geodetic

   def user_con_event(phase, Variables, params):
       """
       User-defined constraints at event boundaries.
       
       Args:
           phase: Phase number
           Variables: Dictionary containing state variables
           params: Parameters dictionary
           
       Returns:
           List of constraint values (>= 0 satisfied)
       """
       con_list = []
       
       # Example: Minimum altitude at FAIRING event
       if params['event_name'] == 'FAIRING':
           altitude = params['altitude']
           altitude_min = 90000  # meters
           con_list.append(altitude - altitude_min)
       
       return con_list

   def user_con_path(phase, t, Variables, params):
       """
       User-defined path constraints throughout trajectory.
       
       Args:
           phase: Phase number
           t: Time
           Variables: Dictionary containing state variables
           params: Parameters dictionary
           
       Returns:
           List of constraint values (>= 0 satisfied)
       """
       con_list = []
       
       # Example: Maximum heating rate constraint
       rho = params['rho']  # Atmospheric density
       v = np.linalg.norm(Variables['vel_NED'])
       q_heating = 0.5 * rho * v**3 / 1e6  # MW/m²
       q_max = 1.0  # MW/m²
       con_list.append(q_max - q_heating)
       
       return con_list

Available Parameters
^^^^^^^^^^^^^^^^^^^^

The ``params`` dictionary provides access to:

* ``rho``: Atmospheric density [kg/m³]
* ``altitude``: Altitude above sea level [m]
* ``event_name``: Current event name (for event constraints)
* ``mass``: Current vehicle mass [kg]
* ``quat``: Attitude quaternion
* Aerodynamic coefficients, dynamic pressure, etc.

Aerodynamic Data (example-CA.csv)
----------------------------------

Axial force coefficient as a function of Mach number:

.. code-block:: text

   Mach,CA
   0.0,0.35
   0.5,0.36
   0.8,0.38
   1.0,0.45
   1.2,0.50
   2.0,0.48
   5.0,0.42

GELATO interpolates this data during trajectory calculation.

Output Files
------------

After successful optimization, the following files are generated in the output directory:

Trajectory Result CSV
^^^^^^^^^^^^^^^^^^^^^

``example-trajectoryResult.csv`` contains the complete trajectory data:

* Time [s]
* Position (latitude, longitude, altitude)
* Velocity (ECEF and NED frames)
* Attitude (quaternion, Euler angles)
* Mass [kg]
* Thrust [N]
* Aerodynamic forces and moments
* Dynamic pressure, Mach number, angle of attack
* And more...

KML Visualization
^^^^^^^^^^^^^^^^^

``example.kml`` can be opened in Google Earth to visualize:

* 3D ground track
* Altitude profile
* Velocity vectors
* Event markers

Optimization Log
^^^^^^^^^^^^^^^^

``example-log.txt`` contains:

* Optimizer convergence information
* Constraint violation history
* Final objective function value
* Computation time

Advanced Examples
-----------------

Multi-Pass Optimization
^^^^^^^^^^^^^^^^^^^^^^^

For complex trajectories, you can perform multi-pass optimization:

1. **Pass 1**: Relax constraints, get feasible solution
   
   .. code-block:: bash
   
      python Trajectory_Optimization.py example-settings.json

2. **Pass 2**: Use Pass 1 result as initial guess, tighten constraints
   
   Update settings to use previous result:
   
   .. code-block:: json
   
      {
        "Initial trajectory file": "example-trajectoryResult.csv"
      }

Parameter Studies
^^^^^^^^^^^^^^^^^

To study the effect of different parameters, create multiple settings files
and use batch processing:

.. code-block:: bash

   # Create folder with multiple settings files
   mkdir parameter_study
   cp example-settings.json parameter_study/case1.json
   cp example-settings.json parameter_study/case2.json
   # Edit case2.json to change parameters
   
   # Run batch optimization
   ./run_batch.sh parameter_study

Tips for Setting Up Your Own Problem
-------------------------------------

1. **Start with the example**: Copy and modify the example files rather than starting from scratch

2. **Initial trajectory**: Provide a good initial guess using ``example-trajectory_init.csv`` or from a previous optimization

3. **Number of nodes**: Start with fewer collocation nodes (5-8) for faster iterations, then increase for final optimization

4. **Constraint relaxation**: Begin with relaxed constraints and gradually tighten them

5. **Event times**: Leave event times empty in the CSV file to allow the optimizer to find optimal timing

6. **Check feasibility**: Ensure your terminal conditions are physically achievable with the given rocket performance

7. **Debugging**: Use ``IPOPT_`` options to increase output verbosity:

   .. code-block:: json
   
      {
        "IPOPT_": {
          "print_level": 5,
          "output_file": "ipopt_debug.txt"
        }
      }

