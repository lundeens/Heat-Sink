#!/usr/bin/env python3
# Sam Lundeen
# Term Project

import math


def hs_analysis(temp_base, temp_ambient, n_fins):
    # Fluid & Solid thermal properties
    temp_film = (temp_base + temp_ambient) / 2
    k_al = 150  # W/m*k
    k_air = 0.02719  # W/m*k, evaluated @ temp_film
    density = 1.04  # kg/m^3
    c_p = 1.007  # kJ/kg*k
    visc_d = 190.27 * 10 ** -7  # mu, N*s/m^2
    visc_k = 17.1 * 10 ** -6  # nu, m^2/s
    diff_therm = 24.28 * 10 ** -6  # alpha, m^2/s
    pr = 0.7  # Prandtl Number

    # Heat Sink Geometry
    diam_hyd = 0.003269  # hydraulic diameter, m
    surf_rough = 0.00015  # Ra / epsilon, m

    # Thermal Analysis
    v = list(range(15, 50, 5))  # velocity range
    re = reynolds(v, density, diam_hyd, visc_d)  # Reynolds Number
    f = darcy(re, surf_rough, diam_hyd)  # Darcy Friction Factor
    nu = nusselt(f, re, pr)  # Nusselt Number, internal pipe flow correlation
    coeff_hx = hxc(nu, k_air, diam_hyd)  # convective heat transfer coefficient, h, W/m^2*k
    return


def reynolds(v, density, diam_hyd, visc_d):
    re = []
    for i in range(0, len(v)):
        re.append(density * v[i] * diam_hyd / visc_d)
    return re


def darcy(re, surf_rough, diam_hyd):
    f = []
    for i in range(0, len(re)):
        guess = 0.1  # arbitrary guess to begin the calculation
        new_guess = (1 / (-2 * math.log10((surf_rough / 2.7 * diam_hyd)) + (2.51 / re[i] * math.sqrt(guess)))) ** 2  # new guess
        while abs(guess - new_guess) > 0.001 * guess:
            guess = new_guess
            new_guess = (1 / (-2 * math.log10((surf_rough / (3.7 * diam_hyd)) + (2.51 / (re[i] * math.sqrt(guess)))))) ** 2
        f.append(new_guess)
    return f


def nusselt(f, re, pr):
    nu = []
    for i in range(0, len(f)):
        nu.append((f[i] / 8) * (re[i] - 1000) * pr / (1 + 12.7 * math.sqrt(f[i] / 8) * ((pr ** (2 / 3)) - 1)))
    return nu


def hxc(nu, k_air, diam_hyd):
    hxc = []
    for i in range(0, len(nu)):
        hxc.append(nu[i] * k_air / diam_hyd)
    return hxc


if __name__ == '__main__':
    hs_analysis(393, 311, 180)
