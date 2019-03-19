#!/usr/bin/env python3
# Sam Lundeen
# Term Project, functions
import math


# Functions
def size(hsn):  # basic geometry : (length, spacing, width, thickness, # of fins) [m]
    if hsn == 0:
        g_b = (0.0089, 0.0018, 0.0078, 0.0007, 180)
    if hsn == 1:
        geo_fp = (ra(hsn), 0.0508, 0.00225806)  # (Ra, Lc, A) for flat plate
        return geo_fp
    if 2 <= hsn <= 15:
        g_b = (1, 0.13125, 2, .1, 7)
    if hsn > 15:
        return print('heat sink number not recognized (geo.size)')
    # calculated geometry: (perimeter, l_c, a_space, a_cross-s, a_base) [m, m^2], fixed tuple
    g_c = (2 * g_b[2] + 2 * g_b[3], 2 * g_b[1], g_b[0] * g_b[1], g_b[2] * g_b[3], g_b[1] * g_b[2] * (g_b[4] + 1))
    a_fin = sa(hsn)                                             # fin area, [m^2]
    a_tot = g_b[4] * a_fin + g_c[4]                             # fin area + base area, [m^2]
    d_h = 4 * g_c[2] / (2 * g_b[0] + g_b[1])                    # hydraulic diameter, [m] only used for f
    sr = ra(hsn)
    # (hydraulic diameter, surface roughness, perimeter, lc, cross-section, length, fin area, # of fins, total area)
    geo = (d_h, sr, g_c[0], g_c[1], g_c[3], g_b[0], a_fin, g_b[4], g_c[4], a_tot)
    return geo


def ra(hsn):  # returns the value of the surface roughness for a given heat sink number
    ra = [0.00015, 0.000005, 0.000005, 0.0000707484, 0.000031586, 0.000017961, 0.000031586, 0.000031586, 0.000017961,
          0.0000707484, 0.000031586, 0.000017961, 0.000031586, 0.000031586, 0.000017961]
    return ra[hsn]


def sa(hsn):  # returns the exposed surface area of each heat sink
    sa = [0.00015676, fp(), 4.4, 4.745053246, 4.570185148, 4.504747482, 4.570185148, 4.570185148, 4.504747482, 5.6691,
          4.80404, 4.57612, 4.80404, 4.80404, 4.57612]
    return sa[hsn]


def fp():
    return ValueError('Must use flat plate correlations for HSN 1')
