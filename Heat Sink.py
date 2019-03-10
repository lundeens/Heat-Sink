#!/usr/bin/env python3
# Sam Lundeen
# Term Project, analysis

from CoolProp.CoolProp import PropsSI
import pandas
import geo
import ind
import math
from statistics import mean
import matplotlib.pyplot as plt
import numpy


def hs_analysis(ve, p, hsn, data):
    # import data
    data = pandas.read_csv(data)                                 # read the function input
    cu_top_ = data['Channel 2'].values                           # store the temperature in a list
    cu_midd = data['Channel 3'].values
    cu_bott = data['Channel 4'].values
    cu_left = data['Channel 5'].values
    cu_back = data['Channel 6'].values
    cu_righ = data['Channel 7'].values
    cu_fron = data['Channel 8'].values
    de_bott = data['Channel 9'].values
    hs_base = data['Channel 10'].values
    ambient = data['Channel 11'].values
    cu_base = []
    for i in range(0, len(cu_back)):
        cu_base.append((cu_left[i] + cu_back[i] + cu_righ[i] + cu_bott[i]) / 4)
    t = []
    for i in range(1, len(cu_top_) - 301):                      # find the steady state values of the data
        if abs(cu_top_[i] - cu_top_[i + 300]) <= 0.01:
            t.append(i)
    temp_cu_top_ = mean(cu_top_[t[0]:-1])
    temp_cu_midd = mean(cu_midd[t[0]:-1])
    temp_cu_bott = mean(cu_bott[t[0]:-1])
    temp_cu = mean(cu_base[t[0]:-1])
    temp_delb = mean(de_bott[t[0]:-1])
    temp_base = mean(hs_base[t[0]:-1])
    temp_amb = mean(ambient[t[0]:-1])
    print(temp_cu, temp_delb, temp_base, temp_amb)

    # base, ambient temp [k], atmospheric pressure [mb] (https://w1.weather.gov/obhistory/KCVO.html), heat sink #
    # fluid & solid thermal properties
    dt = temp_base - temp_amb                                    # temp difference between base and ambient, [k]
    temp_film = (temp_base + temp_amb) / 2                       # film temperature, [k]
    k_al = ind.kal(temp_base)                                    # thermal conductivity of aluminum @ base_temp, [W/m*k]
    k_air = PropsSI('L', 'T', temp_film, 'P', p*100, 'Air')      # thermal conductivity of air @ temp_film, [W/m*k]
    k_del = 0.4                                                  # thermal conductivity of delrin, [W/m*k]
    dens = PropsSI('D', 'T', temp_amb, 'P', p*100, 'Air')        # air density @ temp_amb, [kg/m^3]
    visc = PropsSI('V', 'T', temp_amb, 'P', p*100, 'Air')        # air dynamic viscosity mu @ temp_amb, [N*s/m^2]
    pr = PropsSI('Prandtl', 'T', temp_film, 'P', p*100, 'Air')   # Prandtl Number @ temp_film, [unit-less]

    # Heat Sink Geometry
    g = geo.size(hsn)                                            # surface roughness Ra, [micron]

    # Thermal Loss
    q = 56.4515  # power supplied to the system for 2.5 W/cm^2 [W]
    l_b = (0.00258064 / 0.0127) * k_del * (temp_cu - temp_delb)  # thermal loss through bottom delrin [W]
    print(l_b)
    re_d = dens * ve * 0.1016 / visc  # Reynolds Number at width scale [unit-less]
    s = 2 * math.pi * 0.0381 / (0.93 * math.log(4 / 2) - 0.05)  # shape factor for conduction model, [m]
    n = k_air * 0.158 * (re_d ** (2 / 3)) * (pr ** (1 / 3)) * 0.0154838 / (s * k_del * 0.1016)  # loss coefficient, [unit-less]
    print(k_air, re_d, pr, s, k_del, n)
    dn = 1
    while abs(dn) > 0.0000001:  # loop to find proper k and Pr
        temp_surf = (temp_cu + n * temp_amb) / (1 + n)  # delrin side surface temperature, [k]
        temp_filmd = (temp_surf + temp_amb) / 2  # film temperature on the delrin sides, [k]
        k_aird = PropsSI('L', 'T', temp_filmd, 'P', p * 100, 'Air')  # thermal conductivity of air @ temp_filmd, [W/m*k]
        prd = PropsSI('Prandtl', 'T', temp_filmd, 'P', p * 100, 'Air')  # Prandtl Number @ temp_filmd, [unit-less]
        nd = k_aird * 0.158 * (re_d ** (2 / 3)) * (prd ** (1 / 3)) * 0.0154838 / (s * k_del * 0.1016)
        dn = nd - n
        n = nd
    print(temp_surf, n)
    l_s = s * k_del * (temp_cu - temp_surf)
    q_leftover = q - l_b - l_s
    print(q_leftover)
    ss = 2 * math.pi * (0.0000254) / (0.93 * math.log(4 / 2) - 0.05)  # shape factor for conduction model, [m]
    ns = k_aird * 0.158 * (re_d ** (2 / 3)) * (prd ** (1 / 3)) * 0.00001032 / (
                ss * k_del * 0.1016)  # loss coefficient, [unit-less]
    print(ns)
    temp_gcu = [temp_cu]
    temp_gsu = [temp_surf]
    for i in range(2625):
        dn = 1
        while abs(dn) > 0.0000001:  # loop to find proper k and Pr by refining ns
            temp_surf = (temp_gcu[i] + ns * temp_amb) / (1 + ns)  # delrin side surface temperature, [k]
            # print(temp_surf)
            temp_filmd = (temp_surf + temp_amb) / 2  # film temperature on the delrin sides, [k]
            k_aird = PropsSI('L', 'T', temp_filmd, 'P', p * 100,
                             'Air')  # thermal conductivity of air @ temp_filmd, [W/m*k]
            prd = PropsSI('Prandtl', 'T', temp_filmd, 'P', p * 100, 'Air')  # Prandtl Number @ temp_filmd, [unit-less]
            nd = k_aird * 0.158 * (re_d ** (2 / 3)) * (prd ** (1 / 3)) * 0.00001032 / (ss * k_del * 0.1016)
            dn = nd - ns
            ns = nd
        # print(ns)
        temp_gsu.append(temp_surf)
        less = s * k_del * (temp_gcu[i] - temp_gsu[i])
        # print(less)
        q_leftover -= less
        temp_gcu.append(temp_gcu[i] - 2.54 * 10 ** -5 * q_leftover / (ind.kcu(temp_gcu[i]) * 0.00258064))
    # print(q_leftover)
    print(temp_gsu)
    print(temp_gcu)

    # Thermal Analysis
    v = list(range(15, 50, 5))                                   # velocity range, [m/s] NEEDS UPDATING
    re = ind.reynolds(v, dens, g[3], visc)                       # Reynolds Number (l_c), [unit-less] @ T_film
    f = ind.darcy(re, g[1], g[0])                                # Darcy Friction Factor, [unit-less]
    nu = ind.nusselt(f, re, pr)                                  # Nusselt Number, internal pipe, [unit-less] @ T_fim
    coeff_hx = ind.hxc(nu, k_air, g[3])                          # convective heat transfer coefficient, [h, W/m^2*k]
    m = ind.little_m(coeff_hx, g[2], k_al, g[4])                 # temperature profile constant, [m^-1]
    em = ind.big_m(coeff_hx, g[2], k_al, g[4], dt)               # fin heat rate constant, [W]
    qf = ind.q_f(em, m, g[5], coeff_hx, k_al)                    # single fin heat transfer rate, [W]
    eff = ind.effectiveness(qf, coeff_hx, g[4], dt)              # fin effectiveness, [unit-less]
    eta = ind.efficiency(qf, coeff_hx, g[6], dt)                 # fin efficiency, [unit-less]
    res_fin = ind.rfin(dt, qf)                                   # single fin resistance, [k/W]
    res_base = ind.rbase(coeff_hx, g[4])                         # base resistance, [k/W]
    qtot = ind.q_tot(g[7], eta, coeff_hx, g[6], dt, g[8])        # total heat transfer, [W]
    etao = ind.eta_o(qtot, coeff_hx, g[9], dt)                   # overall efficiency, [unit-less]
    res_o = ind.roverall(dt, qtot)                               # overall resistance, [k/W]

    # plots
    time = numpy.arange(0, len(cu_top_), 1)  # x range, 0 to 4 pi
    plt.plot(time, cu_top_, 'r', time, cu_midd, 'b', time, cu_bott, 'g', time, cu_base, 'k', time, hs_base, 'c')  # plot the data
    plt.grid(True, alpha=0.25)  # turn on the grid with light lines
    plt.legend(('Top', 'Middle', 'Bottom', 'Copper Base', 'Heat Sink Base'), loc='best')
    plt.xlabel('Time (s)')  # label the plot areas
    plt.ylabel('Temperature (k)')
    plt.title('Temperature of the Copper Block Over Time')
    plt.show()
    return etao, eta, res_o


if __name__ == '__main__':
    #  (ve, temp_dels, p, hsn, data)
    e = hs_analysis(8.9408, 1016.2, 2, '0_20_02-28-2019_1129_Lundeen.csv')




