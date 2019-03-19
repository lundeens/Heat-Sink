#!/usr/bin/env python3
# Sam Lundeen
# Term Project, old loss model

re_d = dens * ve * 0.1016 / visc                             # Reynolds Number at width scale [unit-less]
    s = 2 * math.pi * 0.0381 / (0.93 * math.log(4 / 2) - 0.05)   # shape factor for conduction model, [m]
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
    q_leftover = q - l_bo - l_s
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
