make_kml
========

KML generation module for converting trajectory to Google Earth format.

.. note::
   This module requires trajectory result CSV files to be present in the output directory.
   It is designed to be run as a standalone script after trajectory optimization.

Module Overview
---------------

The ``make_kml`` module converts trajectory optimization results from CSV format
to KML (Keyhole Markup Language) files for visualization in Google Earth.

Key Features:
    * Generate 3D trajectory paths with altitude information
    * Create event markers for stage separation and other key points
    * Visualize Instantaneous Impact Point (IIP) trajectories
    * Support for custom styling and colors

Usage:
    Run the script after trajectory optimization completes::

        python make_kml.py

    The generated KML file can be opened in Google Earth for 3D visualization.
