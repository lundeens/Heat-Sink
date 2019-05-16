#!/usr/bin/env python3
# Sam Lundeen
# Term Project, functions
import math


# Functions
def kal(t):  # returns the thermal conductivity of aluminum for a given temperature, Al 6061
    return 101.42 + 0.7709*t - ((7.731*10**-5)*t**2) - ((7.514*10**-6)*t**3) + ((1.934*10**-8)*t**4) - ((1.973*10**-11)*t**5) + ((7.406*10**-15)*t**6)


def kcu(t):
    return 399.5 - 0.05618*t + ((3.518*10**-5)*t**2) - ((7.531*10**-8)*t**3) + ((3.338*10**-11)*t**4)


def reynolds(v, density, diam_hyd, visc_d):  # finds the reynolds number by calling reynolds()
    return density * v * diam_hyd / visc_d


def darcy(re, surf_rough, diam_hyd):  # finds the darcy friction factor
    guess = 0.1  # arbitrary guess to begin the calculation
    new_guess = (1 / (-2 * math.log10((surf_rough / 2.7 * diam_hyd)) + (2.51 / re * math.sqrt(guess)))) ** 2
    while abs(guess - new_guess) > 0.001 * guess:
        guess = new_guess
        new_guess = (1 / (-2 * math.log10((surf_rough / (3.7 * diam_hyd)) + (2.51 / (re * math.sqrt(guess)))))) ** 2
    return new_guess


def nusselt(f, re, pr):  # finds the nusselt number
    return (f / 8) * (re - 1000) * pr / (1 + 12.7 * math.sqrt(f / 8) * ((pr ** (2 / 3)) - 1))


def nusselt2(re, pr):  # finds the nusselt number without the friction factor
    return 0.023 * (re ** 0.8) * (pr ** 0.33)


def hxc(nu, k_air, diam_hyd):  # finds the convective heat transfer coefficient
    return nu * k_air / diam_hyd


def little_m(h, p, k, a_c):  # finds the fin temperature profile constant
    return math.sqrt(h * p / (k * a_c))


def big_m(h, p, k, a_c, dtb):  # finds the fin heat rate constant
    return math.sqrt(h * p * k * a_c) * dtb


def q_f(em, m, l, h, k):  # finds the fin heat rate
    return em * (math.sinh(m * l) + (h / (m * k)) * math.cosh(m * l)) / (math.cosh(m * l) + (h / (m * k)) * math.sinh(m * l))


def effectiveness(qf, h, a_c, dtb):  # finds the single fin effectiveness
    return qf / (h * a_c * dtb)


def efficiency(qf, h, a_f, dtb):
    return qf / (h * a_f * dtb)


def rfin(dtb, qf):
    return dtb / qf


def rbase(h, a_c):
    return 1/(h * a_c)


def q_tot(n, eta, h, a_f, dtb, a_b):
    return n * eta * h * a_f * dtb + h * a_b * dtb


def eta_o(qt, h, a_tot, dtb):
    return qt / (h * a_tot * dtb)


def roverall(dtb, qt):
    return dtb / qt


if __name__ == '__main__':
    print(kal(390))
    print(kcu(420))
