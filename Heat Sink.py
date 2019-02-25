#!/usr/bin/env python3
# Sam Lundeen
# Term Project, analysis

from CoolProp.CoolProp import PropsSI
import geo
import ind


def hs_analysis(temp_base, temp_amb, p, hsn):
    # base, ambient temp [k], atmospheric pressure [mb] (https://w1.weather.gov/obhistory/KCVO.html), heat sink #
    # fluid & solid thermal properties
    dt = temp_base - temp_amb                                   # temp difference between base and ambient, [k]
    temp_film = (temp_base + temp_amb) / 2                      # film temperature, [k]
    k_al = ind.kal(temp_base)                                   # thermal conductivity of aluminum, [W/m*k]
    k_air = PropsSI('L', 'T', temp_film, 'P', p*100, 'Air')     # thermal conductivity of air @ temp_film, [W/m*k]
    dens = PropsSI('D', 'T', temp_film, 'P', p*100, 'Air')      # air density, [kg/m^3]
    visc = PropsSI('V', 'T', temp_film, 'P', p*100, 'Air')      # air dynamic viscosity mu, [N*s/m^2]
    pr = PropsSI('Prandtl', 'T', temp_film, 'P', p*100, 'Air')  # Prandtl Number, [unit-less]

    # heat sink geometry
    g = geo.size(hsn)                                           # surface roughness Ra, [micron]

    # Thermal Analysis
    v = list(range(15, 50, 5))                                  # velocity range, [m/s] NEEDS UPDATING
    re = ind.reynolds(v, dens, g[3], visc)                      # Reynolds Number (l_c), [unit-less] @ T_film
    f = ind.darcy(re, g[1], g[0])                               # Darcy Friction Factor, [unit-less]
    nu = ind.nusselt(f, re, pr)                                 # Nusselt Number, internal pipe, [unit-less] @ T_fim
    coeff_hx = ind.hxc(nu, k_air, g[3])                         # convective heat transfer coefficient, [h, W/m^2*k]
    m = ind.little_m(coeff_hx, g[2], k_al, g[4])                # temperature profile constant, [m^-1]
    em = ind.big_m(coeff_hx, g[2], k_al, g[4], dt)              # fin heat rate constant, [W]
    qf = ind.q_f(em, m, g[5], coeff_hx, k_al)                   # single fin heat transfer rate, [W]
    eff = ind.effectiveness(qf, coeff_hx, g[4], dt)             # fin effectiveness, [unit-less]
    eta = ind.efficiency(qf, coeff_hx, g[6], dt)                # fin efficiency, [unit-less]
    res_fin = ind.rfin(dt, qf)                                  # single fin resistance, [k/W]
    res_base = ind.rbase(coeff_hx, g[4])                        # base resistance, [k/W]
    qtot = ind.q_tot(g[7], eta, coeff_hx, g[6], dt, g[8])       # total heat transfer, [W]
    etao = ind.eta_o(qtot, coeff_hx, g[9], dt)                  # overall efficiency, [unit-less]
    res_o = ind.roverall(dt, qtot)                              # overall resistance, [k/W]
    return (etao, eta)


if __name__ == '__main__':

    e = hs_analysis(393, 311, 1016.2, 0)
    print(e[0])
    print(e[1])

