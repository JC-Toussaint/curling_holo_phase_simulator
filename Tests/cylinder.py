#!/usr/bin/python3

from meshMaker import Cylinder

Cylinder(
    radius = 100,  # radius [nm]
    thickness = 200,  # length [nm]
    mesh_size = 4,
    surfName = "surface",
    volName = "cyl_volume").make("cyl.msh")

