#!/usr/bin/env python3
# Sam Lundeen
# Term Project, analysis

from CoolProp.CoolProp import PropsSI
import ind


def hs_analysis(temp_base, temp_amb, p, hsn):
    # base, ambient temp [k], atmospheric pressure [mb] (https://w1.weather.gov/obhistory/KCVO.html), heat sink #
    # fluid & solid thermal properties
    dt = temp_base - temp_amb                                   # temp difference between base and ambient, [k]
    temp_film = (temp_base + temp_amb) / 2                      # film temperature, [k]
    k_al = ind.kal(temp_base)                                   # thermal conductivity of aluminum, [W/m*k]
    k_air = PropsSI('L', 'T', temp_film, 'P', p*100, 'Air')     # thermal conductivity of air @ temp_film, [W/m*k]
    dens = PropsSI('D', 'T', temp_amb, 'P', p*100, 'Air')       # air density, [kg/m^3]
    visc = PropsSI('V', 'T', temp_amb, 'P', p*100, 'Air')       # air dynamic viscosity mu, [N*s/m^2]
    pr = PropsSI('Prandtl', 'T', temp_amb, 'P', p*100, 'Air')   # Prandtl Number, [unit-less]

    # heat sink geometry
    g_b = (0.0089, 0.0018, 0.0078, 0.0007, 180)  # basic geometry : (length, spacing, width, thickness, # of fins) [m]
    # calculated geometry: (perimeter, l_c, a_space, a_cross-s, a_base) [m, m^2], fixed tuple, l_c NEEDS UPDATING, 2*s
    g_c = (2 * g_b[2] + 2 * g_b[3], g_b[0] + g_b[3] / 2, g_b[0] * g_b[1], g_b[2] * g_b[3], g_b[1] * g_b[2] * (g_b[4] + 1))
    a_fin = 2 * g_b[0] * g_b[2] + 2 * g_b[0] * g_b[3] + g_c[3]  # fin area, [m^2] NEEDS UPDATING, based on real surface
    a_tot = g_b[4] * a_fin + g_c[4]                             # fin area + base area, [m^2]
    d_h = 4 * g_c[2] / (2 * g_b[0] + g_b[1])                    # hydraulic diameter, [m] NEEDS UPDATING, parallel plate, change to L-C
    ra = ind.heatsink(hsn)                                      # surface roughness Ra, [micron] NEEDS ERROR MESSAGE

    # Thermal Analysis
    v = list(range(15, 50, 5))                                  # velocity range, [m/s] (V_1, V_N + 5, increment)
    re = ind.reynolds(v, dens, d_h, visc)                       # Reynolds Number, [unit-less]
    f = ind.darcy(re, ra, d_h)                                  # Darcy Friction Factor, [unit-less]
    nu = ind.nusselt(f, re, pr)                                 # Nusselt Number, internal pipe, [unit-less]
    coeff_hx = ind.hxc(nu, k_air, d_h)                          # convective heat transfer coefficient, [h, W/m^2*k]
    m = ind.little_m(coeff_hx, g_c[0], k_al, g_c[3])            # temperature profile constant, [m^-1]
    em = ind.big_m(coeff_hx, g_c[0], k_al, g_c[3], dt)          # fin heat rate constant, [W]
    qf = ind.q_f(em, m, g_b[0], coeff_hx, k_al)                 # single fin heat transfer rate, [W]
    eff = ind.effectiveness(qf, coeff_hx, g_c[3], dt)           # fin effectiveness, [unit-less]
    eta = ind.efficiency(qf, coeff_hx, a_fin, dt)               # fin efficiency, [unit-less]
    res_fin = ind.rfin(dt, qf)                                  # single fin resistance, [k/W]
    res_base = ind.rbase(coeff_hx, g_c[3])                      # base resistance, [k/W]
    qtot = ind.q_tot(g_b[4], eta, coeff_hx, a_fin, dt, g_c[4])  # total heat transfer, [W]
    etao = ind.eta_o(qtot, coeff_hx, a_tot, dt)                 # overall efficiency, [unit-less]
    res_o = ind.roverall(dt, qtot)                              # overall resistance, [k/W]
    return etao


if __name__ == '__main__':
    e = hs_analysis(393, 311, 1016.2, 0)
    print(e)

