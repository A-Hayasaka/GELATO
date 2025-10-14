# %%
"""
KML Generation Module
=====================

Module for converting trajectory optimization results to KML files for Google Earth.

This module loads trajectory data from CSV files and converts them to KML format
that can be visualized in Google Earth.

Main Features:
    * Generate trajectory linestrings (lines)
    * Generate markers for event points (stage separation, etc.)
    * Visualize IIP (Instantaneous Impact Point) trajectory
    * Add altitude information for 3D display

Functions:
    km_point_save: Save event points to KML
    kml_folder_save: Save trajectory folder to KML
"""

import pandas as pd

# load csv file
out = pd.read_csv("./output/RP_ORv2A_SYS_0133_B-trajectoryResult.csv")


event_time = out.dropna(subset=["event"]).set_index("event").loc[:, "time"].to_dict()


def km_point_save(f, event, lat, lon, alt):
    f.write("\t\t\t<Placemark>\n")
    f.write(f"\t\t\t\t<name>{event}</name>\n")
    f.write("\t\t\t\t<Point>\n")
    f.write("\t\t\t\t\t<altitudeMode>absolute</altitudeMode>\n")
    f.write(f"\t\t\t\t\t<coordinates>{lon},{lat},{alt}</coordinates>\n")
    f.write("\t\t\t\t</Point>\n")
    f.write("\t\t\t</Placemark>\n")


def kml_folder_save(f, name, df):
    f.write("\t\t<Folder>\n")
    f.write(f"\t\t\t<name>{name}</name>\n")
    f.write("\t\t\t<Placemark>\n")
    f.write(f"\t\t\t\t<name>{name}</name>\n")
    f.write("\t\t\t\t<LineString>\n")
    f.write("\t\t\t\t\t<altitudeMode>absolute</altitudeMode>\n")
    f.write("\t\t\t\t\t<coordinates>")
    for index, row in df.iterrows():
        if name == "PPI":
            lat = row["lat"]
            lon = row["lon"]
            alt = row["altitude"]
        elif name == "IIP":
            lat = row["lat_IIP"]
            lon = row["lon_IIP"]
            alt = 0.0
        f.write(f" {lon},{lat},{alt}")
    f.write("</coordinates>\n")
    f.write("\t\t\t\t</LineString>\n")
    f.write("\t\t\t</Placemark>\n")
    for event in [
        "LIFTOFF",
        "TKT",
        "Q10S",
        "TGTS",
        "H8KM",
        "TGTE",
        "MEPCO",
        "Q11E",
        "MECO",
        "MSEP",
        "SELI1",
        "Q21S",
        "FRGDR",
        "IIP1",
        "SECO1",
        "SCSEP",
        "EOS",
    ]:
        if name == "PPI":
            lat, lon, alt = df[df["time"] == event_time[event]][
                ["lat", "lon", "altitude"]
            ].values[0]
        elif name == "IIP":
            lat, lon = df[df["time"] == event_time[event]][
                ["lat_IIP", "lon_IIP"]
            ].values[0]
            alt = 0.0
        km_point_save(f, event, lat, lon, alt)

    f.write("\t\t</Folder>\n")


kml_file = "./output/RP_ORv2A_SYS_0133_B-trajectoryResult.kml"

with open(kml_file, "w") as f:
    f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
    f.write(
        '<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2">'
    )
    f.write("\t<Document>\n")
    kml_folder_save(f, "PPI", out)
    kml_folder_save(f, "IIP", out)
    f.write("\t</Document>\n")
    f.write("</kml>\n")
