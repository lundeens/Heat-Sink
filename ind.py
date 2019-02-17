#!/usr/bin/env python3
# Sam Lundeen
# Term Project, functions
import math


# Functions
def heatsink(hsn):  # returns the value of the surface roughness for a given heat sink number
    hs = [0.00015, 0.000005, 0.000005, 0.0000707484, 0.000031586, 0.000017961, 0.000031586, 0.000031586, 0.000017961,
          0.0000707484, 0.000031586, 0.000017961, 0.000031586, 0.000031586, 0.000017961]
    return hs[hsn]


def reynolds(v, density, diam_hyd, visc_d):  # finds the reynolds number by calling reynolds()
    re = []
    for i in range(0, len(v)):
        re.append(density * v[i] * diam_hyd / visc_d)
    return re


def darcy(re, surf_rough, diam_hyd):  # finds the darcy friction factor
    f = []
    for i in range(0, len(re)):
        guess = 0.1  # arbitrary guess to begin the calculation
        new_guess = (1 / (-2 * math.log10((surf_rough / 2.7 * diam_hyd)) + (2.51 / re[i] * math.sqrt(guess)))) ** 2
        while abs(guess - new_guess) > 0.001 * guess:
            guess = new_guess
            new_guess = (1 / (-2 * math.log10((surf_rough / (3.7 * diam_hyd)) + (2.51 / (re[i] * math.sqrt(guess)))))) ** 2
        f.append(new_guess)
    return f


def nusselt(f, re, pr):  # finds the nusselt number
    nu = []
    for i in range(0, len(f)):
        nu.append((f[i] / 8) * (re[i] - 1000) * pr / (1 + 12.7 * math.sqrt(f[i] / 8) * ((pr ** (2 / 3)) - 1)))
    return nu


def hxc(nu, k_air, diam_hyd):  # finds the convective heat transfer coefficient
    hxc = []
    for i in range(0, len(nu)):
        hxc.append(nu[i] * k_air / diam_hyd)
    return hxc


def little_m(h, p, k, a_c):  # finds the fin temperature profile constant
    m = []
    for i in range(0, len(h)):
        m.append(math.sqrt(h[i] * p / (k * a_c)))
    return m


def big_m(h, p, k, a_c, dtb):  # finds the fin heat rate constant
    m = []
    for i in range(0, len(h)):
        m.append(math.sqrt(h[i] * p * k * a_c) * dtb)
    return m


def q_f(em, m, l, h, k):  # finds the fin heat rate
    qf = []
    for i in range(0, len(em)):
        qf.append(em[i] * (math.sinh(m[i] * l) + (h[i] / (m[i] * k)) * math.cosh(m[i] * l)) / (math.cosh(m[i] * l) + (h[i] / (m[i] * k)) * math.sinh(m[i] * l)))
    return qf


def effectiveness(qf, h, a_c, dtb):  # finds the single fin effectiveness
    eff = []
    for i in range(0, len(qf)):
        eff.append(qf[i] / (h[i] * a_c * dtb))
    return eff


def efficiency(qf, h, a_f, dtb):
    eff = []
    for i in range(0, len(qf)):
        eff.append(qf[i] / (h[i] * a_f * dtb))
    return eff


def rfin(dtb, qf):
    r = []
    for i in range(0, len(qf)):
        r.append(dtb / qf[i])
    return r


def rbase(h, a_c):
    r = []
    for i in range(0, len(h)):
        r.append(1/(h[i] * a_c))
    return r


def q_tot(n, eta, h, a_f, dtb, a_b):
    q = []
    for i in range(0, len(eta)):
        q.append(n * eta[i] * h[i] * a_f * dtb + h[i] * a_b * dtb)
    return q


def eta_o(qt, h, a_tot, dtb):
    e = []
    for i in range(0, len(qt)):
        e.append(qt[i] / (h[i] * a_tot * dtb))
    return e


def roverall(dtb, qt):
    r = []
    for i in range(0, len(qt)):
        r.append(dtb / qt[i])
    return r
