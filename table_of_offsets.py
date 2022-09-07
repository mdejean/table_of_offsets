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

parser = argparse.ArgumentParser(description="Generate a table of offsets for a hull from a mesh")

parser.add_argument('mesh', metavar='file', help="A mesh in a supported format (STL, OBJ, etc.)")
parser.add_argument('--stations', '-s', metavar='sta', nargs='+', default=[8], help="Number of stations, station spacing, or a list of station locations (with units)")
parser.add_argument('--waterlines', '-w', metavar='wl', nargs='+', default=[6], help="Waterline spacing or a list of waterline locations")
parser.add_argument('--buttocks', '-b', metavar='butt', nargs='+', default=[4], help="Buttock spacing or a list of buttock locations")
parser.add_argument('--format', '-f', default='fie', choices=['fie', 'fie+', 'mm'], help="Unit format to use for the table")
parser.add_argument('--chines', action='store_true', help="Attempt to automatically locate chines")
parser.add_argument('--show-chines', action='store_true', help="Show chine-locating graph")

args = parser.parse_args()

def get_offset(mesh, station, butt, wl, accept_angle=None, ap=False):
    #a coordinate outside the mesh
    outside = mesh.scale
    
    origin = [station, butt, wl]
    direction = [0, 0, 0]
    normal = [1, 0, 0] #normal of the plane the offset is in, used for angle
    
    #for i in range(3):
    #    if origin[i] is None:
    #        direction[i] = -1
    #        origin[i] = outside
    #        break
    if butt is None: #in from +y
        direction[1] = -1
        origin[1] = outside
    elif wl is None: #up from base
        direction[2] = 1
        origin[2] = -outside
    elif station is None:
        if ap:
            direction[0] = -1
            origin[0] = outside
        else:
            direction[0] = 1
            origin[0] = -outside
        if butt == 0:
            normal = [0, 0, 1]
        else:
            normal = [0, -1, 0]

    locations, index_ray, index_tri = mesh.ray.intersects_location(ray_origins=[origin], ray_directions=[direction])
    if len(locations):
        if accept_angle:
            reflection_angle = \
                np.degrees(
                    (
                        np.arctan2(
                            np.dot(mesh.face_normals[index_tri][0], direction),
                            np.dot(
                                normal, 
                                np.cross(mesh.face_normals[index_tri][0], direction)
                            )
                        ) 
                        + 2 * np.pi
                    ) % (2 * np.pi)
                )
            
            if reflection_angle < accept_angle[0] or reflection_angle > accept_angle[1]:
                raise ValueError
        if butt is None:
            return max(locations[:,1])
        elif wl is None:
            return min(locations[:,2])
        elif station is None:
            if ap:
                return max(locations[:,0])
            else:
                return min(locations[:,0])
    else:
        return None

mesh = trimesh.load_mesh(args.mesh)
mesh.fix_normals()


#todo: reorient the mesh
# x = widest dimension = 0 @ bow (the pointier end)
# y = 2nd widest = port/starboard (the second widest dimension) 0 @ centerline (autodetect half/full?)
# z = narrow dimension = 0 @ keel (also the pointier end)

#mesh.rezero()
if mesh.extents[1] < mesh.extents[1]:
    pass

loa = mesh.bounds[1][0]
beam = mesh.bounds[1][1] * 2
depth = mesh.bounds[1][2]

epsilon = 1e-5

def parse_section(a, dimension):
    if len(a) == 1:
        q = Quantity(a[0])
        if not q.units:
            #number of stations equally spaced between FP and AP
            sections = np.linspace(0, dimension, round(q.real) + 2)[1:-1] 
        else:
            #spacing between stations, centered between FP and AP
            q = Quantity(q, scale='m')
            sections = np.arange((dimension % q.real)/2, dimension, q.real)
    else:
        sections = []
        for s in a:
            q = Quantity(s, scale='m')
            sections.append(q.real)
    return sections

stations = parse_section(args.stations, loa)
butts = parse_section(args.buttocks, beam/2)
wls = parse_section(args.waterlines, depth)

no_chines = 0

half_breadths = {'fp': {}, 'ap': {}}
heights = {'fp': {}, 'ap': {}}
fp_dist = {'wl': {}, 'butt': {}}
ap_dist = {'wl': {}, 'butt': {}}

section = mesh.section([1, 0, 0], [epsilon, 0, 0])
heights['fp']['sheer'] = section.bounds[1][2]
half_breadths['fp']['sheer'] = section.bounds[1][1]
fp_dist['wl']['sheer'] = 0

#should we consider reverse transoms?
section = mesh.section([1, 0, 0], [loa-epsilon, 0, 0])
heights['ap']['sheer'] = section.bounds[1][2]
half_breadths['ap']['sheer'] = section.bounds[1][1]
ap_dist['wl']['sheer'] = loa

fp_dist['wl']['keel'] = loa/2
h = epsilon
while True:
    h += epsilon
    try:
        dist = get_offset(mesh, None, 0, h, accept_angle=(45, 360-45))
        if dist:
            fp_dist['wl']['keel'] = dist
            break
    except ValueError:
        pass
heights['fp']['keel'] = h
half_breadths['fp']['keel'] = get_offset(mesh, fp_dist['wl']['keel'], None, h) 

ap_dist['wl']['keel'] = loa/2
h = epsilon
while True:
    h += epsilon
    try:
        dist = get_offset(mesh, None, 0, h, accept_angle=(45, 360-45), ap=True)
        if dist:
            ap_dist['wl']['keel'] = dist
            break
    except ValueError:
        pass
heights['ap']['keel'] = h
half_breadths['ap']['keel'] = get_offset(mesh, ap_dist['wl']['keel'], None, h) 

for station in stations:
    half_breadths[station] = {}
    heights[station] = {}
    
    for wl in wls:
        half_breadths[station][wl] = get_offset(mesh, station, None, wl)
    
    for butt in butts:
        heights[station][butt] = get_offset(mesh, station, butt, None)
    
    section = mesh.section([1, 0, 0], [station, 0, 0])
    
    heights[station]['keel'] = section.bounds[0][2]
    half_breadths[station]['keel'] = get_offset(mesh, station, None, heights[station]['keel'] + epsilon)
    
    # ignore the 'top' of the boat, even if it's not flat
    #  by moving down from the top of the model if the ray hits the mesh
    #  at the candidate sheer at less than a 45 degree angle
    #  this will handle sprung decks, but not wheelhouses or anything complicated
    h = section.bounds[1][2]
    while True:
        h -= epsilon
        try:
            half_breadths[station]['sheer'] = get_offset(mesh, station, None, h, accept_angle=(45, 360-45))
            break
        except ValueError:
            pass
    
    heights[station]['sheer']  = h
    
    if args.chines:
        ys = section.discrete[0][:, 1]
        zs = section.discrete[0][:, 2]
        
        dy_dt = np.gradient(ys)
        dz_dt = np.gradient(zs)
        curvature = np.gradient(np.divide(dz_dt, dy_dt))
        
        #peaks = signal.find_peaks_cwt(np.square(curvature), np.arange(20, 25), min_snr=1)
        peaks, props = signal.find_peaks(np.square(curvature)/np.sum(np.square(curvature)), prominence=0.01)
        if args.show_chines:
            plt.subplot(211)
            plt.scatter(np.arange(len(curvature)), np.square(curvature), c = ['red' if k in peaks else 'blue' for k in range(len(curvature))])
            plt.subplot(212)
            plt.scatter(ys, zs, marker = '+', c = ['red' if k in peaks else 'blue' for k in range(len(ys))])
            plt.show()
        
        chine = 1
        
        for peak in peaks:
            #only consider the positive half of the model
            if ys[peak] < 0:
                continue
            
            #ignore the sheer and/or keel chines
            if np.abs(zs[peak] - heights[station]['sheer']) < epsilon or \
                np.abs(zs[peak] - heights[station]['keel']) < epsilon:
                continue
            
            heights[station]['chine %d' % chine] = zs[peak]
            half_breadths[station]['chine %d' % chine] = ys[peak]
            
            chine += 1
        no_chines = max(no_chines, chine - 1)

for wl in wls:
    section = mesh.section([0, 0, 1], [0, 0, wl])
    fp_dist['wl'][wl] = section.bounds[0][0]
    ap_dist['wl'][wl] = section.bounds[1][0]
    
    half_breadths['fp'][wl] = get_offset(mesh, fp_dist['wl'][wl]+epsilon, None, wl)
    half_breadths['ap'][wl] = get_offset(mesh, ap_dist['wl'][wl]-epsilon, None, wl)

for butt in butts:
    section = mesh.section([0, 1, 0], [0, butt, 0])
    fp_dist['butt'][butt] = section.bounds[0][0]
    ap_dist['butt'][butt] = section.bounds[1][0]
    
    heights['fp'][butt] = get_offset(mesh, fp_dist['butt'][butt]+epsilon, butt, None)
    heights['ap'][butt] = get_offset(mesh, ap_dist['butt'][butt]-epsilon, butt, None)




fmt = args.format

def format_offset(a, fmt):
    if a is None:
        return ''
    if fmt == 'fie' or fmt == 'fie+':
        eighths = a / 0.3048 * 12 * 8
        plus = ''
        if fmt == 'fie+':
            if eighths % 1 > 0.5 and eighths % 1 < 0.75:
                plus = '-'
            if eighths % 1 < 0.5 and eighths % 1 > 0.25:
                plus = '+'
        eighths = round(eighths)
        inches = eighths // 8
        eighths = eighths % 8
        feet = inches // 12
        inches = inches % 12
        
        return '%d-%d-%d%s' % (feet, inches, eighths, plus)
    if fmt == 'mm':
        return '%d' % round(a * 1000)

def format(a, fmt):
    if fmt == 'fie':
        eighths = round(a / 0.3048 * 12 * 8)
        inches = eighths // 8
        eighths = eighths % 8
        div = math.gcd(8, eighths)
        num = eighths/div
        denom = 8/div
        return "%d%s\"" % (inches, "-%d/%d" % (num, denom) if num else '') 
    elif fmt == 'fie+':
        s = round(a / 0.3048 * 12 * 16)
        inches = s // 16
        s = s % 16
        div = math.gcd(16, s)
        num = s/div
        denom = 16/div
        return "%d%s\"" % (inches, "-%d/%d" % (num, denom) if num else '') 
    if fmt == 'mm':
        return Quantity(a, scale='m').render()

table = np.empty(
    ( 
      # header   sheer  wls        chines      keel sheer butts        chines      keel
        2      + 1    + len(wls) + no_chines + 1  + 1   + len(butts) + no_chines + 1,
      # header   fp  stations        ap   fp_dist   ap_dist
        2      + 1 + len(stations) + 1  + 1       + 1
    ),
    dtype='U25'
    )

col = 1

table[0, 1] = 'Station'
table[1, 1] = 'Distance from FP'

col += 1
table[0, col] = 'FP'
table[1, col] = format(0, fmt)

for i, station in enumerate(stations):
    col += 1
    table[0, col] = '%d' % (i + 1)
    table[1, col] = format(station, fmt)
    
    table[2, col] = format_offset(half_breadths[station]['sheer'], fmt)
    table[1 + 1 + len(wls) + no_chines + 1, col] = format_offset(half_breadths[station]['keel'], fmt)
    table[1 + 1 + len(wls) + no_chines + 1 + 1, col] = format_offset(heights[station]['sheer'], fmt)
    table[1 + 1 + len(wls) + no_chines + 1 + 1 + len(butts) + no_chines + 1, col] = format_offset(heights[station]['keel'], fmt)

col += 1
table[0, col] = 'AP'
table[1, col] = format(loa, fmt)

col += 1
table[0, col] = 'Position of'
table[1, col] = 'FP'
col += 1
table[1, col] = 'AP'

row = 2

col = 0
table[row, col] = 'Half-breadths'
col += 1
table[row, col] = 'Sheer'
col += 1
table[row, col] = format_offset(half_breadths['fp']['sheer'], fmt)
for k, station in enumerate(stations):
    col += 1
    table[row, col] = format_offset(half_breadths[station]['sheer'], fmt)
col += 1
table[row, col] = format_offset(half_breadths['ap']['sheer'], fmt)
col += 1
table[row, col] = format_offset(fp_dist['wl']['sheer'], fmt)
col += 1
table[row, col] = format_offset(loa - ap_dist['wl']['sheer'], fmt)

for wl in reversed(wls):
    row += 1
    col = 1
    table[row, col] = '%s above base' % format(wl, fmt)
    
    col += 1
    table[row, col] = format_offset(half_breadths['fp'][wl], fmt)
    
    for station in stations:
        col += 1
        table[row, col] = format_offset(half_breadths[station][wl], fmt)
    
    col += 1
    table[row, col] = format_offset(half_breadths['ap'][wl], fmt)
    
    col += 1
    table[row, col] = format_offset(fp_dist['wl'][wl], fmt)
    col += 1
    table[row, col] = format_offset(loa - ap_dist['wl'][wl], fmt)

#todo: chine at ap/fp
for i in range(no_chines):
    row += 1
    col = 1
    table[row, col] = 'Chine %d' % (i+1)
    
    col += 1
    #fp
    
    for station in stations:
        col += 1
        table[row, col] = format_offset(half_breadths[station].get('chine %d' % (i+1)), fmt)


row += 1
col = 1
table[row, col] = 'Keel'
col += 1
table[row, col] = format_offset(half_breadths['fp']['keel'], fmt)
for k, station in enumerate(stations):
    col += 1
    table[row, col] = format_offset(half_breadths[station]['keel'], fmt)
col += 1
table[row, col] = format_offset(half_breadths['ap']['keel'], fmt)
col += 1
table[row, col] = format_offset(fp_dist['wl']['keel'], fmt)
col += 1
table[row, col] = format_offset(loa - ap_dist['wl']['keel'], fmt)

row += 1
col = 0
table[row, col] = 'Heights above base'
col += 1
table[row, col] = 'Sheer'
col += 1
table[row, col] = format_offset(heights['fp']['sheer'], fmt)
for k, station in enumerate(stations):
    col += 1
    table[row, col] = format_offset(heights[station]['sheer'], fmt)
col += 1
table[row, col] = format_offset(heights['ap']['sheer'], fmt)
col += 1
table[row, col] = format_offset(fp_dist['wl']['sheer'], fmt)
col += 1
table[row, col] = format_offset(loa - ap_dist['wl']['sheer'], fmt)

for butt in reversed(butts):
    row += 1
    col = 1
    table[row, col] = '%s buttock' % format(butt, fmt)
    
    col += 1
    table[row, col] = format_offset(heights['fp'][butt], fmt)
    
    for k, station in enumerate(stations):
        col += 1
        table[row, col] = format_offset(heights[station][butt], fmt)
    
    col += 1
    table[row, col] = format_offset(heights['ap'][butt], fmt)
    
    col += 1
    table[row, col] = format_offset(fp_dist['butt'][butt], fmt)
    col += 1
    table[row, col] = format_offset(loa - ap_dist['butt'][butt], fmt)

for i in range(no_chines):
    row += 1
    col = 1
    table[row, col] = 'Chine %d' % (i+1)
    
    col += 1
    #fp
    
    for k, station in enumerate(stations):
        col += 1
        table[row, col] = format_offset(heights[station].get('chine %d' % (i+1)), fmt)

row += 1
col = 1
table[row, col] = 'Keel'
col += 1
table[row, col] = format_offset(heights['fp']['keel'], fmt)
for station in stations:
    col += 1
    table[row, col] = format_offset(heights[station]['keel'], fmt)
col += 1
table[row, col] = format_offset(heights['ap']['keel'], fmt)
col += 1
table[row, col] = format_offset(fp_dist['wl']['keel'], fmt)
col += 1
table[row, col] = format_offset(loa - ap_dist['wl']['keel'], fmt)

print("\n".join([",".join(list(t)) for t in table]))