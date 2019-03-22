#!/usr/bin/env python3
# Sam Lundeen
# Term Project, analysis

from CoolProp.CoolProp import PropsSI
import pandas
import geo
import ind
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
        cu_base.append((cu_left[i] + cu_back[i] + cu_righ[i] + cu_fron[i]) / 4)
    t = []
    for i in range(1, len(cu_top_) - 301):                      # find the steady state values of the data
        if abs(cu_top_[i] - cu_top_[i + 300]) <= 0.01:
            t.append(i)
    temp_cu_top_ = mean(cu_top_[t[0]:-1])                       # top temperature in copper, [k]
    temp_cu_midd = mean(cu_midd[t[0]:-1])                       # middle temperature in copper, [k]
    temp_cu_bott = mean(cu_bott[t[0]:-1])                       # bottom temperature in copper, [k]
    temp_cu = mean(cu_base[t[0]:-1])                            # base temperature of copper, [k]
    temp_delb = mean(de_bott[t[0]:-1])                          # bottom of delrin temperature, [k]
    temp_base = mean(hs_base[t[0]:-1])                          # heat sink base temperature, [k]
    temp_amb = mean(ambient[t[0]:-1])                           # ambient temperature, [k]
    #print(temp_cu, temp_delb, temp_base, temp_amb)

    # Temperature plots
    time = numpy.arange(0, len(cu_top_), 1)  # x range, 0 to 4 pi
    plt.plot(time, cu_top_, 'r', time, cu_midd, 'b', time, cu_bott, 'g', time, cu_base, 'k', time, hs_base, 'c')  # plot the data
    plt.grid(True, alpha=0.25)  # turn on the grid with light lines
    plt.legend(('Top', 'Middle', 'Bottom', 'Copper Base', 'Heat Sink Base'), loc='best')
    plt.xlabel('Time (s)')  # label the plot areas
    plt.ylabel('Temperature (k)')
    plt.title('Temperature of the Copper Block Over Time')
    plt.show()

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
    # l_b = (0.00258064 / 0.0127) * k_del * (temp_cu - temp_delb)    # thermal loss through bottom delrin [W]
    # loss out the bottom
    temp_film_bo = (temp_delb + temp_amb) / 2                        # film temperature at bottom, [k]
    k_air_bo = PropsSI('L', 'T', temp_film_bo, 'P', p*100, 'Air')    # thermal conductivity of air, [W/m*k]
    rho_bo = PropsSI('D', 'T', temp_film_bo, 'P', p*100, 'Air')      # density of air, [kg/m^3]
    cp_bo = PropsSI('CP0MASS', 'T', temp_film_bo, 'P', p*100, 'Air') # specific heat of air, [J/kg*k]
    alpha_bo = k_air_bo / (rho_bo * cp_bo)                           # thermal diffusivity, [m^2/s]
    kvisc_bo = PropsSI('V', 'T', temp_film_bo, 'P', p*100, 'Air') / rho_bo  # kinematic viscosity of air, [m^2/s]
    ra = 9.81 * (1 / temp_film_bo) * (temp_delb - temp_amb) * (0.0508 ** 3) / (alpha_bo * kvisc_bo)  # Rayleigh number [unit-less]
    h_bo = 0.52 * k_air_bo * (ra ** (1 / 5)) / 0.0508                # convective heat transfer coefficient, [W/m^2*k]
    l_bo = 0.00258064 * (temp_cu - temp_amb) / ((1 / h_bo) + (0.0127 / k_del))  # loss through bottom delrin [W]
    #print(l_bo)
    # loss out the base
    k_air_ba = PropsSI('L', 'T', temp_film_bo, 'P', p*100, 'Air')     # thermal conductivity of air, [W/m*k]
    pr_ba = PropsSI('Prandtl', 'T', temp_film_bo, 'P', p*100, 'Air')  # Prandtl number, [unit-less]
    re_ba = dens * ve * 0.1016 / visc                                 # Reynolds number at the base, [unit-less]
    h_ba = k_air_ba * 0.158 * (re_ba ** 0.66) * (pr_ba ** (1 / 3)) / 0.508  # heat transfer coefficient [W/m^2*k]
    l_ba = 0.007974 * (temp_cu - temp_amb) / ((1 / h_ba) + (0.0254 / k_del))  # loss out the base, [W]
    #print(l_ba)
    # loss out the fourier section
    topgrad = (temp_cu_top_ - temp_cu_midd) / 0.018999  # 0.748 in    # map out temperature gradients and average them
    bottomgrad = (temp_cu_midd - temp_cu_bott) / 0.019075  # 0.751 in # to create a reference gradient that will be used
    fullgrad = (temp_cu_top_ - temp_cu_bott) / 0.038075  # 1.499 in   # to find the average temperature on the side
    tgrad = (topgrad + bottomgrad + fullgrad) / 3
    # print(topgrad, bottomgrad, fullgrad, tgrad)
    temp_int = temp_cu_midd - tgrad * 0.037313  # 1.469 in            # top surface temperature estimation
    temp_top = tgrad * 0.064567 + temp_int  # 2.542 in                # bottom of gradient temperature estimation
    l_s = 0.01312 * (((temp_int + temp_top) / 2) - temp_amb) / ((1 / h_ba) + (0.0254 / k_del))  # loss out the side, [W], h = 1.545 in
    #print(l_s)
    loss = l_bo + l_ba + l_s                                     # total loss
    leftover = q - loss                                          # leftover after loss
    qpa = leftover / 25.8064                                     # leftover flux
    print('losses:', q, l_bo, l_ba, l_s, leftover, qpa)

    # Thermal Analysis
    if hsn == 1:                                                 # flat plate correlations
        re = ind.reynolds(ve, dens, g[1], visc)                  # Reynolds number for flat plate
        nu = 0.037 * (re ** (4 / 5)) * (pr ** (1 / 3))           # Nusselt number for turbulent flow over a flat plate
        h = nu * k_air / g[1]                                    # heat transfer coefficient, [W/m^2*k]
        qdp = h * (temp_base - temp_amb)                         # heat flux out the plate, [W/m^2]
        q = qdp * g[2]                                           # heat rate out the plate, [W]
        print('Flat Plate Performance (HSN 1): h=', h, ' qdp=', qdp, 'q=', q)
        perf = (nu, h, qdp, q)                                   # performance tuple
        return perf
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
    eta_o = ind.eta_o(qtot, coeff_hx, g[9], dt)                  # overall efficiency, [unit-less]
    res_o = ind.roverall(dt, qtot)                               # overall resistance, [k/W]

    # Efficiency plot
    plt.plot(v, eta, 'r', v, eta_o, 'b')  # plot the data
    #plt.errorbar()  # https://matplotlib.org/api/_as_gen/matplotlib.pyplot.errorbar.html
    plt.grid(True, alpha=0.25)  # turn on the grid with light lines
    plt.legend(('Fin Efficiency', 'Overall Efficiency'), loc='best')
    plt.xlabel('Velocity (m/s)')  # label the plot areas
    plt.ylabel('Efficiency')
    plt.title('Efficiency vs Air Speed')
    plt.show()

    # Resistance plot
    plt.plot(v, res_fin, 'r', v, res_o, 'b')  # plot the data
    plt.grid(True, alpha=0.25)  # turn on the grid with light lines
    plt.legend(('Fin Resistance', 'Overall Resistance'), loc='best')
    plt.xlabel('Velocity (m/s)')  # label the plot areas
    plt.ylabel('Resistance [k/W]')
    plt.title('Thermal Resistance vs Air Speed')
    plt.show()
    return eta, eta_o, res_fin, res_o


if __name__ == '__main__':
    #  (airspeed, pressure, heat sink number, data)
    e = hs_analysis(8.9408, 1016.2, 2, '0_20_02-28-2019_1129_Lundeen.csv')
    print(e)




