#!/usr/bin/env python3
# Sam Lundeen
# Term Project, analysis

from CoolProp.CoolProp import PropsSI
import math
import pandas
import geo
import ind
from statistics import mean
import datetime
import matplotlib.pyplot as plt
import numpy


def hs_analysis(rdata):

    parse = rdata.split('_')
    hsn = int(parse[0])
    if parse[1] == '14k':
        ve = 31.69
    elif parse[1] == '12k':
        ve = 27.17
    elif parse[1] == '10k':
        ve = 22.64
    elif parse[1] == '8k':
        ve = 18.11
    elif parse[1] == '6k':
        ve = 13.58
    elif parse[1] == '4k':
        ve = 9.06
    else:
        print('Reynolds number input not recognized.')
        return ValueError()
    date = parse[2]
    time = parse[3]
    p_parse = parse[4].split('-')
    if len(p_parse[1]) != 1:
        return ValueError('second argument of power supplied must have 3 digits: 1014-3 ')
    p = int(p_parse[0]) + int(p_parse[1]) / 10
    ps_parse = parse[5].split('-')
    if len(ps_parse[1]) != 3:
        return ValueError('second argument of power supplied must have 3 digits: 50-812 ')
    ps = int(ps_parse[0]) + int(ps_parse[1]) / 1000

    # import data
    data = pandas.read_csv(rdata)                                # read the function input
    cu_top_ = data['top'].values                                # store the temperature in a list
    cu_midd = data['middle'].values
    cu_bott = data['bottom'].values
    cu_left = data['left'].values
    cu_back = data['back'].values
    cu_righ = data['right'].values
    cu_fron = data['front'].values
    de_bott = data['de_bottom'].values
    hs_base = data['base'].values
    ambient = data['ambient'].values

    t = []
    for i in range(1, len(cu_top_) - 301):                      # find the steady state values of the data
        if abs(cu_top_[i] - cu_top_[i + 300]) <= 0.1:
            t.append(i)
    #print(t)

    tcl = mean(cu_left[0:-1])                                   # left copper base temp, [k]
    tcb = mean(cu_back[0:-1])                                   # back copper base temp, [k]
    tcr = mean(cu_righ[0:-1])                                   # right copper base temp, [k]
    tcf = mean(cu_fron[0:-1])                                   # front copper base temp, [k]
    temp_cu = (tcl + tcb + tcf) / 3                             # base temperature of copper, [k]
    temp_cu_top_ = mean(cu_top_[0:-1])                          # top temperature in copper, [k]
    temp_cu_midd = mean(cu_midd[0:-1])                          # middle temperature in copper, [k]
    temp_cu_bott = mean(cu_bott[0:-1])                          # bottom temperature in copper, [k]
    temp_delb = mean(de_bott[0:-1])                             # bottom of delrin temperature, [k]
    temp_base = mean(hs_base[0:-1])                             # heat sink base temperature, [k]
    temp_amb = mean(ambient[0:-1])                              # ambient temperature, [k]
    #print(temp_cu, temp_delb, temp_base, temp_amb)

    # Uncertainties
    u_p_v = 0.027627                                              # velocity prec unc, [m/s] [verified]
    u_b_v = 0.22352                                               # velocity bias unc, [m/s]
    u_v = math.sqrt(u_p_v ** 2 + u_b_v ** 2)                      # velocity unc, [m/s]
    u_b_st = 0.1                                                  # standard bias uncertainty, all but base, [k]
    u_b_temp_base = 0.1                                           # base temp bias uncertainty, [k]
    # u_p_tcf =                                                     # front copper tc prec unc, [k]
    # u_p_tcr =                                                     # right copper tc prec unc, [k]
    # u_p_tcb =                                                     # back copper tc prec unc, [k]
    # u_p_tcl =                                                     # left copper tc prec unc, [k]
    # u_p_tcut =                                                    # top copper tc prec unc, [k]
    # u_p_tcum =                                                    # middle copper tc prec unc, [k]
    # u_p_tcub =                                                    # bottom copper tc prec unc, [k]
    # u_p_tdb =                                                     # delrin bottom tc prec unc, [k]
    # u_p_tba =                                                     # hs base tc prec unc, [k]
    # u_p_tam =                                                     # amb tc prec unc, [k]
    # u_cu = 0.25 * math.sqrt(19.36 + u_p_tcf ** 2 + u_p_tcr ** 2 + u_p_tcb ** 2 + u_p_tcl ** 2)  # avg cu tc unc, [k]

    # Temperature plots
    # time = numpy.arange(0, len(cu_top_), 1)  # x range
    # print(len(time))
    # print(type(cu_top_))

    # plt.plot(time, cu_top_, 'r', time, cu_midd, 'b', time, cu_bott, 'g', time, temp_cu, 'k', time, hs_base, 'c')
    # plt.grid(True, alpha=0.25)  # turn on the grid with light lines
    # plt.legend(('Top', 'Middle', 'Bottom', 'Copper Base', 'Heat Sink Base'), loc='best')
    # plt.xlabel('Time (s)')  # label the plot areas
    # plt.ylabel('Temperature (k)')
    # plt.title('Temperature of the Copper Block Over Time')
    # plt.show()

    # base, ambient temp [k], atmospheric pressure [mb] (https://w1.weather.gov/obhistory/KCVO.html), heat sink #
    # fluid & solid thermal properties
    dt = temp_base - temp_amb                                    # temp difference between base and ambient, [k]
    # u_dt = math.sqrt((u_b_temp_base ** 2) + (u_p_tba ** 2) + (u_b_st ** 2) + (u_p_tam ** 2))
    u_dt = math.sqrt((u_b_temp_base ** 2) + (u_b_st ** 2))
    temp_film = (temp_base + temp_amb) / 2                       # film temperature, [k]
    # u_temp_film = 0.5 * math.sqrt((u_b_temp_base ** 2) + (u_p_tba ** 2) + (u_b_st ** 2) + (u_p_tam ** 2))
    u_temp_film = 0.5 * math.sqrt((u_b_temp_base ** 2) + (u_b_st ** 2))
    k_al = ind.kal(temp_base)                                    # thermal conductivity of aluminum @ base_temp, [W/m*k]
    # u_k_al = max(abs(ind.kal(temp_base) - ind.kal(temp_base + math.sqrt((u_b_temp_base ** 2) + (u_p_tba ** 2)))), abs(ind.kal(temp_base) - ind.kal(temp_base - math.sqrt((u_b_temp_base ** 2) + (u_p_tba ** 2)))))
    u_k_al = max(abs(ind.kal(temp_base) - ind.kal(temp_base + math.sqrt((u_b_temp_base ** 2)))), abs(ind.kal(temp_base) - ind.kal(temp_base - math.sqrt((u_b_temp_base ** 2)))))
    k_air = PropsSI('L', 'T', temp_film, 'P', p*100, 'Air')      # thermal conductivity of air @ temp_film, [W/m*k]
    u_k_air = max(abs(PropsSI('L', 'T', temp_film, 'P', p*100, 'Air') - PropsSI('L', 'T', temp_film + u_temp_film, 'P', p*100, 'Air')), abs(PropsSI('L', 'T', temp_film, 'P', p*100, 'Air') - PropsSI('L', 'T', temp_film - u_temp_film, 'P', p*100, 'Air')))
    k_del = 0.4                                                  # thermal conductivity of delrin, [W/m*k]
    dens = PropsSI('D', 'T', temp_amb, 'P', p*100, 'Air')        # air density @ temp_amb, [kg/m^3]
    # u_dens = max(abs(PropsSI('D', 'T', temp_amb, 'P', p*100, 'Air') - PropsSI('D', 'T', temp_amb + math.sqrt((u_b_st ** 2) + (u_p_tam ** 2)), 'P', p*100, 'Air')), abs(PropsSI('D', 'T', temp_amb, 'P', p*100, 'Air') - PropsSI('D', 'T', temp_amb - math.sqrt((u_b_st ** 2) + (u_p_tam ** 2)), 'P', p*100, 'Air')))
    u_dens = max(abs(PropsSI('D', 'T', temp_amb, 'P', p*100, 'Air') - PropsSI('D', 'T', temp_amb + math.sqrt((u_b_st ** 2)), 'P', p*100, 'Air')), abs(PropsSI('D', 'T', temp_amb, 'P', p*100, 'Air') - PropsSI('D', 'T', temp_amb - math.sqrt((u_b_st ** 2)), 'P', p*100, 'Air')))
    visc = PropsSI('V', 'T', temp_amb, 'P', p*100, 'Air')        # air dynamic viscosity mu @ temp_amb, [N*s/m^2]
    # u_visc = max(abs(PropsSI('V', 'T', temp_amb, 'P', p*100, 'Air') - PropsSI('V', 'T', temp_amb + math.sqrt((u_b_st ** 2) + (u_p_tam ** 2)), 'P', p*100, 'Air')), abs(PropsSI('V', 'T', temp_amb, 'P', p*100, 'Air') - PropsSI('V', 'T', temp_amb - math.sqrt((u_b_st ** 2) + (u_p_tam ** 2)), 'P', p*100, 'Air')))
    u_visc = max(abs(PropsSI('V', 'T', temp_amb, 'P', p*100, 'Air') - PropsSI('V', 'T', temp_amb + math.sqrt((u_b_st ** 2)), 'P', p*100, 'Air')), abs(PropsSI('V', 'T', temp_amb, 'P', p*100, 'Air') - PropsSI('V', 'T', temp_amb - math.sqrt((u_b_st ** 2)), 'P', p*100, 'Air')))
    pr = PropsSI('Prandtl', 'T', temp_film, 'P', p*100, 'Air')   # Prandtl Number @ temp_film, [unit-less]
    u_pr = max(abs(PropsSI('Prandtl', 'T', temp_film, 'P', p*100, 'Air') - PropsSI('Prandtl', 'T', temp_film + u_temp_film, 'P', p*100, 'Air')), abs(PropsSI('Prandtl', 'T', temp_film, 'P', p*100, 'Air') - PropsSI('Prandtl', 'T', temp_film - u_temp_film, 'P', p*100, 'Air')))


    # Heat Sink Geometry
    # geo = (d_h, sr, perimeter, l_c, a_xc, h, a_fin, n, a_base, a_tot, vol, mass)
    # u_geo = (u_b_d_h, u_b_perimeter, u_b_l_c, u_b_a_xc, u_b_st, u_b_a_base, u_b_a_tot)
    g, u_g = geo.size(hsn)                                           # geometry and geometric uncertainties

    # Thermal Loss
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

    # loss out the base
    k_air_ba = PropsSI('L', 'T', temp_film_bo, 'P', p*100, 'Air')     # thermal conductivity of air, [W/m*k]
    pr_ba = PropsSI('Prandtl', 'T', temp_film_bo, 'P', p*100, 'Air')  # Prandtl number, [unit-less]
    re_ba = dens * ve * 0.1016 / visc                                 # Reynolds number at the base, [unit-less]
    h_ba = k_air_ba * 0.158 * (re_ba ** 0.66) * (pr_ba ** (1 / 3)) / 0.0508  # heat transfer coefficient [W/m^2*k]
    l_ba = 0.007974 * (temp_cu - temp_amb) / ((1 / h_ba) + (0.0254 / k_del))  # loss out the base, [W]

    # loss out the fourier section
    topgrad = (temp_cu_top_ - temp_cu_midd) / 0.018999  # 0.748 in    # map out temperature gradients and average them
    bottomgrad = (temp_cu_midd - temp_cu_bott) / 0.019075  # 0.751 in # to create a reference gradient that will be used
    fullgrad = (temp_cu_top_ - temp_cu_bott) / 0.038075  # 1.499 in   # to find the average temperature on the side
    tgrad = (topgrad + bottomgrad + fullgrad) / 3
    temp_int = temp_cu_midd - tgrad * 0.037313  # 1.469 in            # top surface temperature estimation
    temp_top = tgrad * 0.064567 + temp_int  # 2.542 in                # bottom of gradient temperature estimation
    l_s = 0.01312 * (((temp_int + temp_top) / 2) - temp_amb) / ((1 / h_ba) + (0.0254 / k_del))  # loss out the side, [W], h = 1.545 in
    loss = l_bo + l_ba + l_s                                     # total loss
    leftover = ps - loss                                         # leftover after loss
    qpa = leftover / 22.5086                                        # leftover flux
    # print('losses:', q, l_bo, l_ba, l_s, leftover, qpa)

    # Thermal Analysis
    if hsn == 1:                                                 # flat plate correlations
        re = ind.reynolds(ve, dens, g[1], visc)                  # Reynolds number for flat plate
        # u_re_v = 2 * dens * g[1] * math.sqrt(u_p_v ** 2 + u_b_v ** 2) / visc  # v partial for re unc, [unit-less]
        # u_re_dens = 2 * v * g[1] * u_dens / visc                              # dens partial for re unc, [unit-less]
        # u_re_s = 2 * dens * v * u_g[4] / visc                                 # s partial for re unc, [unit-less]
        # u_re_visc = 2 * dens * v * g[1] * u_visc / (visc ** 2)                # visc partial for re unc, [unit-less]
        # u_re = math.sqrt(u_re_v ** 2 + u_re_dens ** 2 + u_re_s ** 2 + u_re_visc ** 2)  # re unc, [unit-less]
        nu = 0.037 * (re ** (4 / 5)) * (pr ** (1 / 3))           # Nusselt number for turbulent flow over a flat plate
        # u_nu_re = 4 * (pr ** (1 / 3)) * u_re / (5 * re ** (1 / 5))  # re partial for nu unc, [?]
        # u_nu_pr = (re ** (4 / 5)) * u_pr / (3 * (pr ** (2 / 3)))  # pr partial for nu unc, [?]
        # u_nu = math.sqrt(u_nu_re ** 2 + u_nu_pr ** 2)  # nu unc, [unit-less]
        h = nu * k_air / g[1]                                    # heat transfer coefficient, [W/m^2*k]
        # u_h = math.sqrt((k_air * u_nu / g[1]) ** 2 + (nu * u_k_air / g[1]) ** 2 + (nu * k_air * u_g[4] / (g[1] ** 2)) ** 2)  # coeff_hx unc, [W/m^2*k]
        qdp = h * (temp_base - temp_amb)                         # heat flux out the plate, [W/m^2]
        # u_qdp = math.sqrt((u_h * (temp_base - temp_amb)) ** 2 + (h * math.sqrt((u_b_temp_base ** 2) + (u_p_tba ** 2))) ** 2 + (h * math.sqrt((u_b_st ** 2) + (u_p_tam ** 2))) ** 2)
        q = qdp * g[2]                                           # heat rate out the plate, [W]
        # u_q = math.sqrt((u_qdp * g[2]) ** 2 + (qdp * g[3]) ** 2)  # q unc, [W]
        print('Flat Plate Performance (HSN 1): h=', h, ' qdp=', qdp, 'q=', q)
        perf = (re, nu, h, qdp, q)                                   # performance tuple
        # unc = (u_re, u_nu, u_h, u_qdp, u_q)                          # uncertainty tuple
        return perf                                # velocity range, [m/s] NEEDS UPDATING
    re = ind.reynolds(ve, dens, g[3], visc)                       # Reynolds Number (l_c), [unit-less] @ T_amb
    u_re_v = 2 * dens * g[3] * math.sqrt(u_v ** 2) / visc  # v partial for re unc, [unit-less]
    u_re_dens = 2 * ve * g[3] * u_dens / visc                              # dens partial for re unc, [unit-less]
    u_re_s = 2 * dens * ve * u_g[4] / visc                                 # s partial for re unc, [unit-less]
    u_re_visc = 2 * dens * ve * g[3] * u_visc / (visc ** 2)                # visc partial for re unc, [unit-less]
    u_re = math.sqrt(u_re_v ** 2 + u_re_dens ** 2 + u_re_s ** 2 + u_re_visc ** 2)  # re unc, [unit-less]

    f = ind.darcy(re, g[1], g[0])                                # Darcy Friction Factor, [unit-less]
    u_f_d_h = 0.0538387 * re * g[1] * math.sqrt(f) * u_g[0] / (g[0] * (g[0] + 0.107677 * re * g[1] * math.sqrt(f)) * (math.log10(g[1] / (3.7 * g[0]) + 2.51 / (re * math.sqrt(f))) ** 3))  # d_h partial for f unc [unit-less]
    u_f_re = 2.65095 * g[0] * u_re / (re * (g[0] + 0.107677 * re * g[1] * math.sqrt(f)) * (math.log10(g[1] / (3.7 * g[0]) + 2.51 / (re * math.sqrt(f))) ** 3))  # re partial for f unc, [unit-less]
    u_f = math.sqrt(u_f_d_h ** 2 + u_f_re ** 2)                  # f unc, [unit-less]

    nu = ind.nusselt(f, re, pr)                                  # Nusselt Number, internal pipe, [unit-less] @ T_film
    u_nu_f = 0.0139194 * pr * (re - 1000) * (math.sqrt(f) * ((pr ** (2 / 3)) - 1) + 0.44542) * u_f / ((math.sqrt(f) * ((pr ** (2 / 3)) - 1) + 0.22271) ** 2)  # f partial for nu unc, [unit-less]
    u_nu_re = f * pr * u_re / (8 * (4.49013 * math.sqrt(f) * ((pr ** (2 / 3)) - 1)) + 1)  # re partial for nu unc [unit-less]
    u_nu_pr = 0.00927962 * (re - 1000) * ((f ** (3 / 2)) * ((pr ** (2 / 3)) - 3) + 0.668132 * f) * u_pr / ((math.sqrt(f) * ((pr ** (2 / 3)) - 1) + 0.22271) ** 2)  # pr partial for nu unc, [unit-less]
    u_nu = math.sqrt(u_nu_f ** 2 + u_nu_re ** 2 + u_nu_pr ** 2)  # nu unc, [unit-less]


    coeff_hx = ind.hxc(nu, k_air, g[3])                          # convective heat transfer coefficient, [h, W/m^2*k]
    u_coeff_hx = math.sqrt((k_air * u_nu / g[3]) ** 2 + (nu * u_k_air / g[3]) ** 2 + (nu * k_air * u_g[2] / (g[3] ** 2)) ** 2)  # coeff_hx unc, [W/m^2*k]

    m = ind.little_m(coeff_hx, g[2], k_al, g[4])                 # temperature profile constant, [m^-1]
    u_m = 0.5 * m * math.sqrt((u_coeff_hx / coeff_hx) ** 2 + (u_g[1] / g[2]) ** 2 + (u_k_al / k_al) ** 2 + (u_g[3] / g[4]) ** 2)  # m unc, [m^-1]

    em = ind.big_m(coeff_hx, g[2], k_al, g[4], dt)               # fin heat rate constant, [W]
    u_em_h = g[2] * k_al * g[4] * dt * u_coeff_hx / (2 * math.sqrt(coeff_hx * g[2] * k_al * g[4]))  # h partial for em unc, [?]
    u_em_p = coeff_hx * k_al * g[4] * dt * u_g[1] / (2 * math.sqrt(coeff_hx * g[2] * k_al * g[4]))  # p partial for em unc, [?]
    u_em_k = g[2] * coeff_hx * g[4] * dt * u_k_al / (2 * math.sqrt(coeff_hx * g[2] * k_al * g[4]))  # k partial for em unc, [?]
    u_em_a_xc = g[2] * k_al * coeff_hx * dt * u_g[3] / (2 * math.sqrt(coeff_hx * g[2] * k_al * g[4]))  # a_xc partial for em unc, [?]
    u_em_dt = math.sqrt(coeff_hx * g[2] * k_al * g[4]) * u_dt  # dt partial for em unc, [?]
    u_em = math.sqrt(u_em_h ** 2 + u_em_p ** 2 + u_em_k ** 2 + u_em_a_xc ** 2 + u_em_dt ** 2)

    qf = ind.q_f(em, m, g[5], coeff_hx, k_al)                    # single fin heat transfer rate, [W]
    print(qf)
    u_qf_em = math.sinh(m * g[5]) + (coeff_hx / (m * k_al)) * math.cosh(m * g[5]) * u_em / math.cosh(m * g[5]) + (coeff_hx / (m * k_al)) * math.sinh(m * g[5])  # em partial for qf unc, [?]
    u_qf_l = em * m * ((m ** 2) - (coeff_hx ** 2) * (k_al ** 2)) * (math.cosh(g[5] * m) ** 2 - math.sinh(g[5] * m) ** 2) * u_g[4] / ((coeff_hx * k_al * math.sinh(g[5] * m) + m * math.cosh(g[5] * m)) ** 2)  # l partial for qf unc, [?]
    u_qf_coeff_hx = k_al * em * m * (math.cosh(g[5] * m) ** 2 - math.sinh(g[5] * m) ** 2) * u_coeff_hx / ((coeff_hx * k_al * math.sinh(g[5] * m) + m * math.cosh(g[5] * m)) ** 2)  # coeff_hx partial for qf unc, [?]
    u_qf_k = coeff_hx * em * m * (math.cosh(g[5] * m) ** 2 - math.sinh(g[5] * m) ** 2) * u_k_al / ((coeff_hx * k_al * math.sinh(g[5] * m) + m * math.cosh(g[5] * m)) ** 2)  # k_al partial for qf unc, [?]
    u_qf_m1 = em * ((-coeff_hx * k_al * math.cosh(g[5] * m) / (m ** 2)) + (coeff_hx * k_al * g[5] * math.sinh(g[5] * m) / m) + (g[5] * math.cosh(g[5] * m))) / ((coeff_hx * k_al * math.sinh(g[5] * m) / m) + math.cosh(g[5] * m))  # first part of m partial for qf unc, [?]
    u_qf_m2 = em * ((coeff_hx * k_al * math.cosh(g[5] * m) / m) + math.sinh(g[5] * m)) * ((-coeff_hx * k_al * math.sinh(g[5] * m) / (m ** 2)) + (coeff_hx * k_al * g[5] * math.cosh(g[5] * m) / m) + (g[5] * math.sinh(g[5] * m))) / (((coeff_hx * k_al * math.sinh(g[5] * m) / m) + math.cosh(g[5] * m)) ** 2)  # 2nd part of m partial for qf unc, [?]
    u_qf_m = (u_qf_m1 - u_qf_m2) * u_m                           # m partial for qf unc, [?]
    u_qf = math.sqrt(u_qf_em ** 2 + u_qf_l ** 2 + u_qf_coeff_hx ** 2 + u_qf_k ** 2 + u_qf_m ** 2)
    print(u_qf_m2)
    print(u_qf)

    eff = ind.effectiveness(qf, coeff_hx, g[4], dt)              # fin effectiveness, [unit-less]
    u_eff_qf = u_qf / (coeff_hx * g[4] * dt)                      # qf partial for eff unc, [?]
    u_eff_coeff_hx = qf * u_coeff_hx / (g[4] * dt * coeff_hx ** 2)  # coeff_hx partial for eff unc, [?]
    u_eff_a_xc = qf * u_g[3] / (coeff_hx * dt * g[4] ** 2)        # a_xc partial for eff unc, [?]
    u_eff_dt = qf * u_dt / (g[4] * coeff_hx * dt ** 2)            # dt partial for eff unc, [?]
    u_eff = math.sqrt(u_eff_qf ** 2 + u_eff_coeff_hx ** 2 + u_eff_a_xc ** 2 + u_eff_dt ** 2)  # eff unc, [unit-less]

    eta = ind.efficiency(qf, coeff_hx, g[6], dt)                 # fin efficiency, [unit-less]
    u_eta_qf = u_qf / (coeff_hx * g[6] * dt)                      # qf partial for eta unc, [?]
    u_eta_coeff_hx = qf * u_coeff_hx / (g[6] * dt * coeff_hx ** 2)  # coeff_hx partial for eta unc, [?]
    u_eta_dt = qf * u_dt / (g[6] * coeff_hx * dt ** 2)            # dt partial for eta unc, [?]
    u_eta = math.sqrt(u_eta_qf ** 2 + u_eta_coeff_hx ** 2 + u_eta_dt ** 2)  # eta unc, [unit-less]

    res_fin = ind.rfin(dt, qf)                                   # single fin resistance, [k/W]
    u_res_fin = math.sqrt((u_dt / qf) ** 2 + (dt * u_qf / (qf ** 2)) ** 2)  # res_fin unc, [k/W]

    res_base = ind.rbase(coeff_hx, g[4])                         # base resistance, [k/W], OBSOLETE

    qtot = ind.q_tot(g[7], eta, coeff_hx, g[6], dt, g[8])        # total heat transfer, [W]
    print(qtot)
    u_qtot_eta = g[7] * coeff_hx * g[6] * dt * u_eta  # eta partial for qtot unc, [?]
    u_qtot_coeff_hx = (g[7] * eta * g[6] * dt + g[8] * dt) * u_coeff_hx  # coeff_hx partial for qtot unc, [?]
    u_qtot_dt = (g[7] * eta * coeff_hx * g[6] + coeff_hx * g[8]) * u_dt  # dt partial for qtot unc, [?]
    u_qtot_a_base = coeff_hx * dt * u_g[5]  # a_base partial for qtot unc, [?]
    u_qtot = math.sqrt(u_qtot_eta ** 2 + u_qtot_coeff_hx ** 2 + u_qtot_dt ** 2 + u_qtot_a_base ** 2)  # qtot unc, [W]

    eta_o = ind.eta_o(qtot, coeff_hx, g[9], dt)                  # overall efficiency, [unit-less]
    u_eta_o_qtot = u_qtot / (coeff_hx * g[9] * dt)                      # qf partial for eta_o unc, [?]
    u_eta_o_coeff_hx = qtot * u_coeff_hx / (g[9] * dt * coeff_hx ** 2)  # coeff_hx partial for eta_o unc, [?]
    u_eta_o_a_tot = qtot * u_g[6] / (coeff_hx * dt * g[9] ** 2)
    u_eta_o_dt = qtot * u_dt / (g[9] * coeff_hx * dt ** 2)            # dt partial for eta_o unc, [?]
    u_eta_o = math.sqrt(u_eta_o_qtot ** 2 + u_eta_o_coeff_hx ** 2 + u_eta_o_dt ** 2 + u_eta_o_a_tot ** 2)  # eta_o unc, [unit-less]

    res_o = ind.roverall(dt, qtot)                               # overall resistance, [k/W]
    print(res_o)
    u_res_o = math.sqrt((u_dt / qtot) ** 2 + (dt * u_qtot / (qtot ** 2)) ** 2)  # res_o unc, [k/W]
    print(u_res_o)

    # Efficiency plot
    # plt.plot(v, eta, 'r', v, eta_o, 'b')  # plot the data
    #plt.errorbar()  # https://matplotlib.org/api/_as_gen/matplotlib.pyplot.errorbar.html
    # plt.grid(True, alpha=0.25)  # turn on the grid with light lines
    # plt.legend(('Fin Efficiency', 'Overall Efficiency'), loc='best')
    # plt.xlabel('Velocity (m/s)')  # label the plot areas
    # plt.ylabel('Efficiency')
    # plt.title('Efficiency vs Air Speed')
    # plt.show()

    # Resistance plot
    # plt.plot(v, res_fin, 'r', v, res_o, 'b')  # plot the data
    # plt.grid(True, alpha=0.25)  # turn on the grid with light lines
    # plt.legend(('Fin Resistance', 'Overall Resistance'), loc='best')
    # plt.xlabel('Velocity (m/s)')  # label the plot areas
    # plt.ylabel('Resistance [k/W]')
    # plt.title('Thermal Resistance vs Air Speed')
    # plt.show()

    id = (hsn, ve, re, p, ps, leftover, qpa)
    pe = (nu, coeff_hx, qf, qtot, eff, eta, eta_o, res_fin, res_o)
    un = (u_v, u_re, u_nu, u_coeff_hx, u_qf, u_qtot, u_eff, u_eta, u_eta_o, u_res_fin, u_res_o)

    dfa = pandas.DataFrame({"File Name": [rdata],
                            "Date": [date],
                            "Time": [time],
                            "hsn": [hsn],
                            "A_tot (m^2)": [g[9]],
                            "Volume (m^3)": [g[10]],
                            "Mass (kg)": [g[11]],
                            "Ra": [g[2]],
                            "ve (m/s)": [ve],
                            "Re": [re],
                            "P (mb)": [p],
                            "ps (W)": [ps],
                            "leftover (W)": [leftover],
                            "q'' (W/m^2)": [qpa],
                            "dTb (k)": [dt],
                            "Nu": [nu],
                            "h (W/m^2*k)": [coeff_hx],
                            "q_fin (W)": [qf],
                            "q_tot (W)": [qtot],
                            "eff": [eff],
                            "eta": [eta],
                            "eta_o": [eta_o],
                            "R_fin (k/W)": [res_fin],
                            "R_0 (k/W)": [res_o],
                            "u_ve": [u_v],
                            "u_re": [u_re],
                            "u_nu": [u_nu],
                            "u_h": [u_coeff_hx],
                            "u_q_fin": [u_qf],
                            "u_q_tot": [u_qtot],
                            "u_eff": [u_eff],
                            "u_eta": [u_eta],
                            "u_eta_o": [u_eta_o],
                            "u_R_fin": [u_res_fin],
                            "u_R_o": [u_res_o]})
    with open('Record.csv', 'a') as f:
        dfa.to_csv(f, header=False, index=False)
    return id, pe


if __name__ == '__main__':
    #  hsn_re_date_time_pressure_power_name.csv
    e = hs_analysis('2_14k_04-6-2019_1516_1013-3_50-799_Lundeen.csv')
