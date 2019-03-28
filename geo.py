#!/usr/bin/env python3
# Sam Lundeen
# Term Project, functions
import math


# Functions
def size(hsn):  # basic geometry : (length, spacing, width, thickness, # of fins) [m]
    if hsn == 0:
        g_b = (0.0089, 0.0018, 0.0078, 0.0007, 180)                     # geometry of test heat sink
    if hsn == 1:
        geo_fp = (ra(hsn), 0.0508, 0.00225806, 1.7145 * 10 ** -6)                          # (Ra, Lc [m], A [m^2], u_A [m^2]) for flat plate
        return geo_fp
    if 2 <= hsn <= 15:
        h = 1 * 0.0254                                                  # fin height (length), [m]
        s = 0.13125 * 0.0254                                            # fin spacing, [m]
        w = 2 * 0.0254                                                  # fin width, [m]
        t = .1 * 0.0254                                                 # fin thickness, [m]
        n = 7                                                           # number of fins, [unit-less]
    if hsn > 15:
        return print('heat sink number not recognized (geo.size)')

    u_b_st = 2.54 * 10 ** -5                                            # standard caliper bias, [m]
    perimeter = 2 * w + 2 * t                                           # space perimeter, [m]
    u_b_perimeter = math.sqrt(2 * (2 * u_b_st ) ** 2)                   # perimeter bias unc, [m]
    l_c = 2 * s                                                         # characteristic length, [m]
    u_b_l_c = 2 * u_b_st                                                # l_c bias unc, [m]
    a_space = h * s                                                     # area of the space, [m^2]
    u_b_a_space = math.sqrt(((s * u_b_st) ** 2) + ((h * u_b_st) ** 2))  # a_space bias unc, [m^2]
    a_xc = w * t                                                        # fin cross-sectional area, [m^2]
    u_b_a_xc = math.sqrt(((w * u_b_st) ** 2) + ((t * u_b_st) ** 2))     # a_xc bias unc, [m^2]
    a_base = s * w * (n + 1)                                            # base area, [m^2]
    u_b_a_base = math.sqrt(((w * (n + 1) * u_b_st) ** 2) + ((s * (n + 1) * u_b_st) ** 2))  # a_base bias unc, [m^2]
    a_fin = sa(hsn)                                                     # fin area, [m^2]
    a_tot = n * a_fin + a_base                                          # fin area + base area, [m^2]
    u_b_a_tot = math.sqrt(((w * (n + 1) * u_b_st) ** 2) + ((s * (n + 1) * u_b_st) ** 2))  # a_tot bias unc, [m^2]
    d_h = 4 * a_space / (2 * h + s)                                     # hydraulic diameter, [m] only used for f
    u_b_d_h = math.sqrt((4 * u_b_a_space / (2 * h + s)) ** 2 + 2 * (4 * a_space * u_b_st / ((2 * h + s) ** 2)) ** 2)  # d_h bias unc, [m]
    sr = ra(hsn)                                                        # surface roughness, [m]
    geo = (d_h, sr, perimeter, l_c, a_xc, h, a_fin, n, a_base, a_tot)
    u_geo = (u_b_d_h, u_b_perimeter, u_b_l_c, u_b_a_xc, u_b_st, u_b_a_base, u_b_a_tot)
    return geo, u_geo


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


if __name__ == '__main__':
    g = size(2)
