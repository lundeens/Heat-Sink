#!/usr/bin/env python3


import csv
import math

def hs_analysis(temp_base, temp_ambient, n_fins):
    # Fluid & Solid thermal properties
    k_al = 150  # W/m*k
    k_air = 0.02719  # W/m*k
    density = 1.04  # kg/m^3
    c_p = 1.007  # kJ/kg*k
    visc_d = 190.27 * 10 ** -7  # mu, N*s/m^2
    visc_k = 17.1 * 10 ** -6  # nu, m^2/s
    diff_therm = 24.28 * 10 ** -6  # alpha, m^2/s
    pr = 0.7  # Prandtl Number

    # Heat Sink Geometry
    diam_hyd = 0.003269  # hydraulic diameter, m

    # Thermal Analysis
    v = list(range(15, 50, 5))  # velocity range
    re = []
    for i in v:
        re.append(i * diam_hyd / visc_k)
    print(re)
    return

def darcy()

if __name__ == '__main__':
    hs_analysis(393, 311, 180)
