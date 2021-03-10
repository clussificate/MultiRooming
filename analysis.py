# -*- coding: utf-8 -*-
"""
@Created at 2021/1/25 15:59
@Author: Kurt
@file:analysis.py
@Desc:
"""
import numpy as np
from main import *
from tqdm import tqdm
from tabulate import tabulate


def parse(result_dict):
    c = result_dict['c']
    profit_u = result_dict['profit_u']
    profit_d = result_dict['profit_d']
    price_u = result_dict['price_u']
    price_d = result_dict['price_d']
    behavior_both = result_dict['behavior_both']
    behavior_off = result_dict['behavior_off']
    behavior_on = result_dict['behavior_on']

    result_string = []
    for s1, s2, s3, s4, s5, s6, s7, s8 in zip(c, price_u, price_d, profit_u, profit_d, behavior_both, behavior_off,
                                              behavior_on):
        result_string.append([s1, s2, s3[0], s3[1], s4, s5, s6, s7, s8])

    return result_string


def utility_compare(Hstore, Honline, Lstore, Lonline):
    if max(Lstore, Lonline) >= 0 and max(Hstore, Honline) >= 0:
        if Lstore >= Lonline:
            return "BS"
        else:
            return "BO"

    elif max(Lstore, Lonline) >= 0 and max(Hstore, Honline) < 0:
        if Lstore >= Lonline:
            return "LS"
        else:
            return "LO"

    elif max(Lstore, Lonline) < 0 and max(Hstore, Honline) >= 0:
        if Hstore >= Honline:
            return "HS"
        else:
            return "HO"

    elif max(Lstore, Lonline) < 0 and max(Hstore, Honline) < 0:
        return "No sales"
    else:
        raise Exception("other cases")


def analyze_behavior(solver, c, con, uniform_price, dual_online_price, dual_offline_price):
    H_u_store, H_u_online = solver.get_utility("H", uniform_price, uniform_price, c, con)
    L_u_store, L_u_online = solver.get_utility("L", uniform_price, uniform_price, c, con)

    H_d_store, H_d_online = solver.get_utility("H", dual_online_price, dual_offline_price, c, con)
    L_d_store, L_d_online = solver.get_utility("L", dual_online_price, dual_offline_price, c, con)

    behavior_uniform = utility_compare(H_u_store, H_u_online, L_u_store, L_u_online)
    behavior_dual = utility_compare(H_d_store, H_d_online, L_d_store, L_d_online)
    return behavior_uniform, behavior_dual


def plot_figures(ALPHA, BETA, V, GAMMA, P_HAT, DELTA_h, DELTA_l, CR, con):
    solver = NumSolver(V, GAMMA, P_HAT, DELTA_h, DELTA_l, CR)
    optimal_profits_u = []
    optimal_prices_u = []

    optimal_profits_d = []
    optimal_prices_d = []

    behavior_tuples_both = []
    behavior_tuples_off = []
    behavior_tuples_on = []

    for current_c in tqdm(np.arange(0, 1, 0.01)):
        optimal_profit_u = 0
        optimal_price_u = 0
        optimal_profit_d = 0
        optimal_online_price_d = 0
        optimal_offline_price_d = 0
        # print("-----------------")
        # print("current c: {:.3f}".format(current_c))
        for p_current in np.arange(0, 1, 0.005):
            # print("current p: {:.3f}".format(p_current))
            profit_current_u = (1 - ALPHA - BETA) * solver.get_profit(p=p_current, poff=p_current, c=current_c,
                                                                      con=con) \
                               + ALPHA * solver.get_profit(p=p_current, poff=p_current, c=0, con=1) \
                               + BETA * solver.get_profit(p=p_current, poff=p_current, c=1, con=0)

            if profit_current_u > optimal_profit_u:
                optimal_profit_u = profit_current_u
                optimal_price_u = p_current
                # print("find a optimal uniform price: {:.3f}, {:.3f}".format(
                #     optimal_price_u, optimal_profit_u))

            for poff_current in np.arange(p_current + con, 1, 0.005):
                profit_current_d = (1 - ALPHA - BETA) * solver.get_profit(p=p_current, poff=poff_current,
                                                                          c=current_c, con=con) \
                                   + ALPHA * solver.get_profit(p=p_current, poff=poff_current, c=0, con=1) \
                                   + BETA * solver.get_profit(p=p_current, poff=poff_current, c=1, con=0)

                if profit_current_d > optimal_profit_d:
                    optimal_profit_d = profit_current_d
                    optimal_online_price_d = p_current
                    optimal_offline_price_d = poff_current
                    # print("find a optimal dual price: {:.3f}, {:.3f}, {:.3f}".format(
                    #     optimal_online_price_d, optimal_offline_price_d, optimal_profit_d))

        behavior_tuple_both = analyze_behavior(solver, current_c, con, optimal_price_u, optimal_online_price_d,
                                               optimal_offline_price_d)
        behavior_tuple_off = analyze_behavior(solver, 0, 1, optimal_price_u, optimal_online_price_d,
                                              optimal_offline_price_d)
        behavior_tuple_on = analyze_behavior(solver, 1, 0, optimal_price_u, optimal_online_price_d,
                                             optimal_offline_price_d)

        optimal_profits_u.append(optimal_profit_u)
        optimal_prices_u.append(optimal_price_u)
        optimal_profits_d.append(optimal_profit_d)
        optimal_prices_d.append((optimal_online_price_d, optimal_offline_price_d))
        behavior_tuples_both.append(behavior_tuple_both)
        behavior_tuples_off.append(behavior_tuple_off)
        behavior_tuples_on.append(behavior_tuple_on)

    res_dict = {"profit_u": optimal_profits_u, "profit_d": optimal_profits_d,
                "price_u": optimal_prices_u, "price_d": optimal_prices_d,
                "behavior_both": behavior_tuples_both, "behavior_off": behavior_tuples_off,
                "behavior_on": behavior_tuples_on}
    return res_dict


if __name__ == "__main__":
    alpha = 0.6
    beta = 0.1
    if alpha + beta > 1:
        raise Exception("alpha + beta should less than 1")
    delta_h = 1
    delta_l = 0.2
    phat = 0.0
    V = 0.7
    gamma = 0.2
    # c = float(self.OfflineCost.get())
    con = 0.3
    cr = 0.1
    DataDict = plot_figures(ALPHA=alpha, BETA=beta, V=V, GAMMA=gamma, P_HAT=phat, DELTA_h=delta_h, DELTA_l=delta_l,
                            CR=cr, con=con)
    keys = DataDict.keys()
    cont = dict([(k, []) for k in keys])
    cont['c'] = []

    flag = False
    for c, profit_u, profit_d, price_u, price_d, \
        behavior_both, behavior_off, behavior_on in zip(np.arange(0, 1, 0.01), DataDict['profit_u'],
                                                        DataDict['profit_d'], DataDict['price_u'],
                                                        DataDict['price_d'],
                                                        DataDict['behavior_both'], DataDict['behavior_off'],
                                                        DataDict['behavior_on']):

        if profit_u > profit_d and price_u < price_d[0]:
            flag = True
            cont['c'].append(c)
            cont['profit_u'].append(profit_u)
            cont['profit_d'].append(profit_d)
            cont['price_u'].append(price_u)
            cont['price_d'].append(price_d)
            cont['behavior_both'].append(behavior_both)
            cont['behavior_off'].append(behavior_off)
            cont['behavior_on'].append(behavior_on)

    if not flag:
        print('Is there a desirable uniform price solution? \n' + " No.")
    else:
        print(
            'Is there a desirable uniform price solution?\n'
            'Yes, e.g.: \n \n ')
        print(
            tabulate(parse(cont),
                     headers=["c", "p", "pon", "poff", "pi_u", "pi_d", "behavior_both", "behavior_off", "behavior_on"]))
        print(
            '\n**Note: BO: both segments buy online; '
            'BS: both segments buy offline;\n'
            'HO: only H type segment buy online; HS: only H type segment buy offline. \n'
            'The 1st element in the behavior tuples is the uniform pricing case; \n'
            'the 2nd element in the behavior tuples is the dual pricing case.')
