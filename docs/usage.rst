Usage
=====

Basic Usage
-----------

GELATO requires the following input files to solve trajectory optimization problems:

1. **Settings File** (JSON format)
2. **Event File** (CSV format)
3. **User Constraints File** (Python format, optional)

Basic Execution
^^^^^^^^^^^^^^^

.. code-block:: bash

   cd /path/to/GELATO
   python Trajectory_Optimization.py [settings_file_name]

Example:

.. code-block:: bash

   cd example
   python ../Trajectory_Optimization.py example-settings.json

Input Files Description
-----------------------

Settings File (JSON)
^^^^^^^^^^^^^^^^^^^^

The settings file describes the basic parameters of the optimization problem. Main items:

* **rocket**: Rocket physical parameters (mass, thrust, specific impulse, etc.)
* **aero**: Aerodynamic parameters (reference area, drag coefficient, etc.)
* **optimize**: Optimization settings (objective function, variable ranges, etc.)
* **output**: Output file settings

Example settings file:

.. code-block:: json

   {
     "name": "example-rocket",
     "stage": 2,
     "rocket": {
       "mass_burnout": [1500, 500],
       "mass_propellant": [3500, 1500],
       "thrust": [50000, 20000],
       "Isp": [280, 320]
     },
     "optimize": {
       "method": "ipopt",
       "maxiter": 1000
     }
   }

Event File (CSV)
^^^^^^^^^^^^^^^^

The event file describes flight events (stage separation, attitude control changes, etc.) in chronological order.

CSV columns:

* **time**: Event occurrence time [s]
* **type**: Event type (staging, pitch_program, etc.)
* **parameter**: Event parameters

User Constraints File (Python)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can add user-defined constraint conditions. Create a ``user_constraints.py`` file and
define constraint functions.

.. code-block:: python

   def user_constraint(phase, t, x, u, param):
       """
       User-defined constraint
       
       Args:
           phase: Phase number
           t: Time
           x: State variables
           u: Control variables
           param: Parameters
           
       Returns:
           Constraint value (satisfies constraint when >= 0)
       """
       # Example: Dynamic pressure constraint
       q = 0.5 * rho * v**2
       q_max = 50000  # Pa
       return q_max - q

Batch Processing
----------------

To process multiple settings files at once, use the batch script.

.. code-block:: bash

   ./run_batch.sh [input_folder]

This script executes optimization for all ``.json`` files in the specified folder.

AWS S3 Integration
^^^^^^^^^^^^^^^^^^

You can also use input files stored on AWS S3:

.. code-block:: bash

   ./run_batch.sh s3://your-bucket/input-folder/

Output Files
------------

When optimization completes, the following files are generated:

* **trajectoryResult.csv**: Trajectory data (time, position, velocity, attitude, etc.)
* **[name].kml**: Trajectory visualization file for Google Earth
* **optimization_log.txt**: Optimization process log

Trajectory Visualization
^^^^^^^^^^^^^^^^^^^^^^^^

Open the KML file in Google Earth to display the optimized trajectory in 3D.

You can also plot using Python scripts:

.. code-block:: bash

   python tools/plot_output.py trajectoryResult.csv

Troubleshooting
---------------

When Optimization Does Not Converge
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Improve initial guess values
* Relax constraint conditions
* Increase the number of optimization iterations
* Adjust solver (IPOPT/SNOPT) parameters

When Computation is Slow
^^^^^^^^^^^^^^^^^^^^^^^^^

* Reduce the number of knot points
* Verify that C++ extension modules are built correctly
* Analyze computation time with profiler

Detailed Examples
-----------------

The ``example/`` folder contains a complete set of input files.
Refer to these to configure your own optimization problems.
