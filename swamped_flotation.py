#!/usr/bin/env python3

import numpy as np
import trimesh
import sys
import math
import argparse
from itertools import *
from scipy import signal
from quantiphy import Quantity

import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Calculate bulkhead locations")
parser.add_argument('mesh', metavar='file', help="A mesh in a supported format (STL, OBJ, etc.)")
parser.add_argument('--freeboard', '-f', metavar='fb', help="Swamped freeboard (with units)")
parser.add_argument('--angle', '-a', metavar='theta', default=90, help="Bulkhead angle (degrees)")

args = parser.parse_args()

mesh = trimesh.load_mesh(args.mesh)
mesh.fix_normals()

epsilon = 1e-6

loa = mesh.bounds[1][0]
beam = mesh.bounds[1][1] * 2
depth = mesh.bounds[1][2]

section = mesh.section([1, 0, 0], [loa/2, 0, 0])
depth_amidships = section.bounds[1][2]

swamped_freeboard = depth_amidships * 0.25

if args.freeboard:
    q = Quantity(args.freeboard)
    q = Quantity(q, scale='m')
    swamped_freeboard = q.real

water_box = trimesh.primitives.Box(extents=[loa, beam, swamped_freeboard],
            transform=trimesh.transformations.compose_matrix(translate=[loa/2, 0, depth_amidships - swamped_freeboard / 2]))
water = mesh.intersection(water_box)
water.fix_normals()
water.fill_holes()

print("water: {}".format(water.volume))

air_box = trimesh.primitives.Box(extents=[loa, beam, depth_amidships - swamped_freeboard],
            transform=trimesh.transformations.compose_matrix(translate=[loa/2, 0, (depth_amidships - swamped_freeboard) / 2 + epsilon]))
air = mesh.intersection(air_box)
air.fix_normals()
air.fill_holes()

air.export("air.stl")

print("air: {}".format(air.volume))

bulkhead_location = loa / (air.volume / water.volume) / 2 
error = 2

while error > 1.01 or error < 0.99: 

    bulkhead_box = trimesh.primitives.Box(
        extents=[10, 10, 10],
        transform=trimesh.transformations.compose_matrix(translate=[-5 + bulkhead_location/2.0, 0, depth/2.0])
        )
    bulkhead = air.intersection(bulkhead_box)
    
    
    #bulkhead.visual.face_colors = [0, 0, 0, 127]
    #air.visual.face_colors = [127, 0, 0, 127]
    #bulkhead_box.visual.face_colors = [0, 127, 0, 127]
    #scene = trimesh.Scene([air, bulkhead, bulkhead_box])
    #scene.show()
    
    fwd_bulkhead_upper = water.intersection(bulkhead_box)

    error = ((water.volume - fwd_bulkhead_upper.volume * 2)/ 2) / bulkhead.volume

    print("bulkhead: {}".format(bulkhead.volume))

    bulkhead_location = bulkhead_location * (1 + error) / 2

print("fwd bulkhead location: {}".format(bulkhead_location))

aft_bulkhead_location = bulkhead_location
error = 2

while error > 1.01 or error < 0.99: 

    aft_bulkhead_box = trimesh.primitives.Box(
        extents=[10, 10, 10],
        transform=trimesh.transformations.compose_matrix(translate=[5 + loa - aft_bulkhead_location/2, 0, depth/2])
        )
    aft_bulkhead = air.intersection(aft_bulkhead_box)
    print("bulkhead: {}".format(aft_bulkhead.volume))
    
    aft_bulkhead_upper = water.intersection(aft_bulkhead_box)
    
    error = ((water.volume - fwd_bulkhead_upper.volume - aft_bulkhead_upper.volume) / 2) / aft_bulkhead.volume
    aft_bulkhead_location = aft_bulkhead_location * (1 + error) / 2
    
print("aft bulkhead location: {}".format(aft_bulkhead_location))
