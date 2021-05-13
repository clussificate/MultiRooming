# -*- coding: utf-8 -*-
"""
@Created at 2021/1/20 21:11
@Author: Kurt
@file:GUI.py
@Desc:
alpha -> offline only consumers
beta -> online only consumers
Considering the offline only consumers, the retailer may charge a high price such that 1-alpha-beta consumers
conduct showrooming.
"""
from tkinter import *

import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import numpy as np
from tqdm import tqdm
from main import utility_compare

matplotlib.use('TkAgg')


def utility_compare_uniform(Hstore, Honline, Lstore, Lonline):
    if max(Lstore, Lonline) >= 0 and max(Hstore, Honline) >= 0:
        if Lstore > Lonline:
            return "BS"
        else:
            return "BO"

    elif max(Lstore, Lonline) >= 0 and max(Hstore, Honline) < 0:
        if Lstore > Lonline:
            return "LS"
        else:
            return "LO"

    elif max(Lstore, Lonline) < 0 and max(Hstore, Honline) >= 0:
        if Hstore > Honline:
            return "HS"
        else:
            return "HO"

    elif max(Lstore, Lonline) < 0 and max(Hstore, Honline) < 0:
        return "No sales"
    else:
        raise Exception("other cases")


def utility_compare_dual(Hstore, Honline, Hshowrooming, Lstore, Lonline, Lshowrooming):
    # due to the the comparison between the three cases is independent of theta, we have:
    if Honline >= max(Hstore, Hshowrooming, 0):  # also implies Honline >= max(Hstore, Hshowrooming)
        if Lonline >= 0:
            return "BO"
        else:
            return "HO"
    elif Hstore >= max(Hshowrooming, 0):
        if Lstore >= 0:
            return "BS"
        else:
            return "HS"

    elif Hshowrooming >= 0:
        if Lshowrooming >= 0:
            return "BSR"
        else:
            return "HSR"
    elif max(Lstore, Lshowrooming, Lonline) < 0:
        return "NoSales"
    else:
        raise Exception("Other cases occur!")


def find_result(DataDict):
    flag = False
    keys = DataDict.keys()
    cont = dict([(k, []) for k in keys])
    cont['c'] = []

    for c, profit_u, profit_d, price_u, price_d, \
        behavior_both, behavior_off, behavior_on in tqdm(zip(np.arange(0, 1, 0.01), DataDict['profit_u'],
                                                             DataDict['profit_d'], DataDict['price_u'],
                                                             DataDict['price_d'],
                                                             DataDict['behavior_both'], DataDict['behavior_off'],
                                                             DataDict['behavior_on'])):

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
    return flag, cont


class NumSolver:
    def __init__(self, V=0.8, GAMMA=0.2, P_HAT=0.1, DELTA_h=0.6, DELTA_l=0.2, CR=0.7):
        self.V = V
        self.GAMMA = GAMMA
        self.P_HAT = P_HAT
        self.DELTA_h = DELTA_h
        self.DELTA_l = DELTA_l
        self.CR = CR

    def get_utility(self, TYPE, p, poff, c, con):
        uty_res = {"store": -1, "online": -1, "showrooming": -1}

        if p >= self.P_HAT:
            delta = self.DELTA_h
        else:
            delta = self.DELTA_l

        if TYPE == "H":
            theta = 1.00
        elif TYPE == "L":
            theta = self.V
        else:
            raise Exception("Consumer Type Error!")

        u_store = (theta - poff) / 2 - c
        u_online = (theta - 2 * p + delta * p) / 2 - con
        u_showrooming = (theta - p - con) / 2 - c

        u_store = int(u_store * 10000 + 0.5) / 10000
        u_online = int(u_online * 10000 + 0.5) / 10000
        u_showrooming = int(u_showrooming * 10000 + 0.5) / 10000

        return u_store, u_online, u_showrooming

    def get_profit(self, p, poff, c, con, logger=False):

        if p >= self.P_HAT:
            delta = self.DELTA_h
        else:
            delta = self.DELTA_l

        H_u_store, H_u_online, H_u_showrooming = self.get_utility("H", p, poff, c, con)
        L_u_store, L_u_online, L_u_showrooming = self.get_utility("L", p, poff, c, con)

        # Low type
        pi_H_showrooming = 1 / 4 * self.GAMMA * p * (
            [1 if (H_u_showrooming >= 0 and H_u_showrooming > max(H_u_store, H_u_online)) else 0][0])
        pi_H_store = 1 / 4 * self.GAMMA * poff * (
            [1 if (H_u_store >= max(H_u_showrooming, 0) and H_u_store > H_u_online) else 0][0])
        pi_H_online = 1 / 4 * self.GAMMA * (p + (1 - delta) * p - delta * self.CR) * (
            [1 if H_u_online >= max(H_u_showrooming, H_u_store, 0) else 0][0])

        # High type
        pi_L_showrooming = 1 / 4 * (1 - self.GAMMA) * p * (
            [1 if (L_u_showrooming >= 0 and L_u_showrooming > max(L_u_store, L_u_online)) else 0][0])
        pi_L_store = 1 / 4 * (1 - self.GAMMA) * poff * (
            [1 if (L_u_store >= max(L_u_showrooming, 0) and L_u_store > L_u_online) else 0][0])
        pi_L_online = 1 / 4 * (1 - self.GAMMA) * (p + (1 - delta) * p - delta * self.CR) * (
            [1 if L_u_online >= max(L_u_showrooming, L_u_store, 0) else 0][0])

        expected_profit = pi_H_online + pi_H_store + pi_H_showrooming + pi_L_online + pi_L_store + pi_L_showrooming

        return expected_profit


def find_result(DataDict):
    flag = False
    keys = DataDict.keys()
    cont = dict([(k, []) for k in keys])
    cont['c'] = []

    for c, profit_u, profit_d, price_u, price_d, \
        behavior_both, behavior_off, behavior_on in tqdm(zip(np.arange(0, 1, 0.01), DataDict['profit_u'],
                                                             DataDict['profit_d'], DataDict['price_u'],
                                                             DataDict['price_d'],
                                                             DataDict['behavior_both'], DataDict['behavior_off'],
                                                             DataDict['behavior_on'])):

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
    return flag, cont


def parse(result_dict):
    c = result_dict['c']
    profit_u = result_dict['profit_u']
    profit_d = result_dict['profit_d']
    price_u = result_dict['price_u']
    price_d = result_dict['price_d']
    behavior_both = result_dict['behavior_both']
    behavior_off = result_dict['behavior_off']
    behavior_on = result_dict['behavior_on']
    result_string = "".join(["{:.3f}, \t {:.3f}, \t {:.3f}, \t {:.3f}, \t {:.3f}, \t {:.3f}, \t {}, \t {}, \t {} \n"
                            .format(s1, s2, s3[0], s3[1], s4, s5, s6, s7, s8) for s1, s2, s3, s4, s5, s6, s7, s8 in
                             zip(c, price_u, price_d, profit_u, profit_d, behavior_both, behavior_off, behavior_on)])

    final_string = "c, \t p, \t pon, \t poff, \t pi_u, \t pi_d, \t behavior_both, \t behavior_off, \t behavior_on \n" \
                   + result_string
    return final_string


class GUI(Frame):
    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.pack()
        self.create_widgets()

    def create_widgets(self):

        self.PropOffOnly = Scale(self, orient=HORIZONTAL, from_=0, to=1,
                                 label="Proportion of the offline-only consumer alpha",
                                 resolution=0.01, length=250)
        self.PropOffOnly.grid(row=0, column=0)
        self.PropOffOnly.set(0)

        self.PropOnOnly = Scale(self, orient=HORIZONTAL, from_=0, to=1,
                                label="Proportion of the online-only consumer beta",
                                resolution=0.01, length=250)
        self.PropOnOnly.grid(row=1, column=0)
        self.PropOnOnly.set(0)

        self.HighReturnProb = Scale(self, orient=HORIZONTAL, from_=0, to=1,
                                    label="High return probability delta_H",
                                    resolution=0.01, length=250)
        self.HighReturnProb.grid(row=2, column=0)
        self.HighReturnProb.set(0.6)

        self.LowReturnProb = Scale(self, orient=HORIZONTAL, from_=0, to=1,
                                   label="Low return probability delta_L",
                                   resolution=0.01, length=250)
        self.LowReturnProb.grid(row=3, column=0)
        self.LowReturnProb.set(0.2)

        self.PriceThreshold = Scale(self, orient=HORIZONTAL, from_=0, to=1,
                                    label="Threshold of price p^hat. ",
                                    resolution=0.01, length=250)
        self.PriceThreshold.grid(row=4, column=0)
        self.PriceThreshold.set(0.1)

        self.ValLow = Scale(self, orient=HORIZONTAL, from_=0, to=1,
                            label="Valuation of low type segment V",
                            resolution=0.01, length=250)
        self.ValLow.grid(row=5, column=0)
        self.ValLow.set(0.8)

        self.PropHigh = Scale(self, orient=HORIZONTAL, from_=0, to=1,
                              label="Proportion of high type segment gamma ",
                              resolution=0.01, length=250)
        self.PropHigh.grid(row=6, column=0)
        self.PropHigh.set(0.2)

        # self.OfflineCost = Scale(self, orient=HORIZONTAL, from_=0, to=1,
        #                label="Offline cost c",
        #                resolution=0.01, length=250)
        # self.OfflineCost.pack()

        self.OnlineCost = Scale(self, orient=HORIZONTAL, from_=0, to=1,
                                label="Online cost con",
                                resolution=0.01, length=250)
        self.OnlineCost.grid(row=7, column=0)
        self.OnlineCost.set(0.16)

        self.ReturnCostRet = Scale(self, orient=HORIZONTAL, from_=0, to=1,
                                   label="The retailer's return cost cr",
                                   resolution=0.01, length=250)
        self.ReturnCostRet.grid(row=8, column=0)
        self.ReturnCostRet.set(0.7)

        self.startButton = Button(self, text='Start', command=self.read_params)
        self.startButton.grid(row=9, column=0)

        self.quitButton = Button(self, text='Quit', command=self.quit)
        self.quitButton.grid(row=10, column=0)

        self.fig = Figure(figsize=(5, 4), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.get_tk_widget().grid(row=0, column=1, columnspan=1, rowspan=6)

        self.ResLabel = Label(self, justify=LEFT)
        self.ResLabel.grid(row=7, column=1, columnspan=1, rowspan=3)

    def plot_figures(self, ALPHA, BETA, V, GAMMA, P_HAT, DELTA_h, DELTA_l, CR, con):
        self.fig.clear()
        self.solver = NumSolver(V, GAMMA, P_HAT, DELTA_h, DELTA_l, CR)

        optimal_profits_u = []
        optimal_prices_u = []

        optimal_profits_d = []
        optimal_prices_d = []

        behavior_tuples_both = []
        behavior_tuples_off = []
        behavior_tuples_on = []

        for current_c in np.arange(0, 0.5, 0.01):
            optimal_profit_u = 0
            optimal_price_u = 0
            optimal_profit_d = 0
            optimal_online_price_d = 0
            optimal_offline_price_d = 0

            for p_current in np.arange(0, 1, 0.005):
                profit_current_u = (1 - ALPHA - BETA) * self.solver.get_profit(p=p_current, poff=p_current, c=current_c,
                                                                               con=con) \
                                   + ALPHA * self.solver.get_profit(p=p_current, poff=p_current, c=0, con=1) \
                                   + BETA * self.solver.get_profit(p=p_current, poff=p_current, c=1, con=0)
                # find optimal uniform prices
                if profit_current_u > optimal_profit_u:
                    optimal_profit_u = profit_current_u
                    optimal_price_u = p_current

                # poff=pon+con, 1-alpha-beta segment does not conduct showrooming poff=V, 1-alpha-beta segment
                # conducts showrooming. The offline optimally serve all offline-only consumers poff=1, 1-alpha-beta
                # segment conducts showrooming. The offline optimally serve high type offline-only consumers
                OfflinePriceScenarios = [p_current + con, V, 1]
                for poff_current in OfflinePriceScenarios:
                    profit_current_d = (1 - ALPHA - BETA) * self.solver.get_profit(p=p_current, poff=poff_current,
                                                                                   c=current_c, con=con) \
                                       + ALPHA * self.solver.get_profit(p=p_current, poff=poff_current, c=0, con=1) \
                                       + BETA * self.solver.get_profit(p=p_current, poff=poff_current, c=1, con=0)

                    if profit_current_d > optimal_profit_d:
                        optimal_profit_d = profit_current_d
                        optimal_online_price_d = p_current
                        optimal_offline_price_d = poff_current

            behavior_tuple_both = self.analyze_behavior(current_c, con, optimal_price_u, optimal_online_price_d,
                                                        optimal_offline_price_d)
            behavior_tuple_off = self.analyze_behavior(0, 1, optimal_price_u, optimal_online_price_d,
                                                       optimal_offline_price_d)
            behavior_tuple_on = self.analyze_behavior(1, 0, optimal_price_u, optimal_online_price_d,
                                                      optimal_offline_price_d)

            optimal_profits_u.append(optimal_profit_u)
            optimal_prices_u.append(optimal_price_u)
            optimal_profits_d.append(optimal_profit_d)
            optimal_prices_d.append((optimal_online_price_d, optimal_offline_price_d))
            behavior_tuples_both.append(behavior_tuple_both)
            behavior_tuples_off.append(behavior_tuple_off)
            behavior_tuples_on.append(behavior_tuple_on)

        for current_c, price_u, price_d, profit_u, profit_d in zip(
                np.arange(0, 1, 0.01), optimal_prices_u, optimal_prices_d, optimal_profits_u, optimal_profits_d):
            behavior_u_both, behavior_d_both = self.analyze_behavior(current_c, con, price_u, price_d[0], price_d[1])
            behavior_u_Off, behavior_d_Off = self.analyze_behavior(0, 1, price_u, price_d[0], price_d[1])
            behavior_u_On, behavior_d_On = self.analyze_behavior(1, 0, price_u, price_d[0], price_d[1])
            print("current c: {}, price_u: {}, price_d: {}".format(current_c, price_u, price_d))
            print("behavior of both-channel u:{}, d:{}".format(behavior_u_both, behavior_d_both))
            print("behavior of offline-only u:{}, d:{}".format(behavior_u_Off, behavior_d_Off))
            print("behavior of online-only u:{}, d:{}".format(behavior_u_On, behavior_d_On))
            print("profit_u: {:.3f}, profit_d: {:.3f}".format(profit_u, profit_d))

        ax1 = self.fig.add_subplot(2, 1, 1)
        ax1.plot(np.arange(0, 0.5, 0.01), optimal_profits_u, 'g', label="Uniform")
        ax1.plot(np.arange(0, 0.5, 0.01), optimal_profits_d, 'k--', label="Dual")
        ax1.set_title("profits")
        ax1.legend(prop=dict(size=6), frameon=False)

        ax2 = self.fig.add_subplot(2, 1, 2)
        ax2.plot(np.arange(0, 0.5, 0.01), optimal_prices_u, 'g', label="Uniform")
        ax2.plot(np.arange(0, 0.5, 0.01), [x[0] for x in optimal_prices_d], 'k--', label="Online of Dual")
        ax2.plot(np.arange(0, 0.5, 0.01), [x[1] for x in optimal_prices_d], linestyle='-.', label="Offline of Dual")
        ax2.set_title("prices")
        ax2.legend(prop=dict(size=6), frameon=False)

        self.fig.tight_layout()
        self.canvas.draw()

        res_dict = {"profit_u": optimal_profits_u, "profit_d": optimal_profits_d,
                    "price_u": optimal_prices_u, "price_d": optimal_prices_d,
                    "behavior_both": behavior_tuples_both, "behavior_off": behavior_tuples_off,
                    "behavior_on": behavior_tuples_on}
        return res_dict

    def read_params(self):
        self.fig.clear()
        alpha = float(self.PropOffOnly.get())
        beta = float(self.PropOnOnly.get())
        if alpha + beta > 1:
            raise Exception("alpha + beta should be less than 1")

        delta_h = float(self.HighReturnProb.get())
        delta_l = float(self.LowReturnProb.get())
        phat = float(self.PriceThreshold.get())
        V = float(self.ValLow.get())
        gamma = float(self.PropHigh.get())
        # c = float(self.OfflineCost.get())
        con = float(self.OnlineCost.get())
        cr = float(self.ReturnCostRet.get())

        data_dict = self.plot_figures(ALPHA=alpha, BETA=beta, V=V, GAMMA=gamma, P_HAT=phat, DELTA_h=delta_h,
                                      DELTA_l=delta_l, CR=cr, con=con)

        result = find_result(data_dict)

        if not result[0]:
            self.ResLabel.config(text='Is there a desirable uniform price solution? \n' + " No.")
        else:
            self.ResLabel.config(
                text='Is there a desirable uniform price solution?\n'
                     'Yes, e.g.: \n \n {} \n'
                     '**Note: BO: both segments buy online; '
                     'BS: both segments buy offline;\n'
                     'HO: only H type segment buy online; HS: only H type segment buy offline. \n'
                     'The 1st element in the behavior tuple is the uniform pricing case; \n'
                     'the 2nd element in the behavior tuple is the dual pricing case'.format(parse(result[1])))

    def analyze_behavior(self, c, con, uniform_price, dual_online_price, dual_offline_price):
        H_u_store, H_u_online, _ = self.solver.get_utility("H", uniform_price, uniform_price, c, con)
        L_u_store, L_u_online, _ = self.solver.get_utility("L", uniform_price, uniform_price, c, con)
        behavior_uniform = utility_compare_uniform(H_u_store, H_u_online, L_u_store, L_u_online)

        H_d_store, H_d_online, H_d_showrooming = self.solver.get_utility("H", dual_online_price, dual_offline_price, c,
                                                                         con)
        L_d_store, L_d_online, L_d_showrooming = self.solver.get_utility("L", dual_online_price, dual_offline_price, c,
                                                                         con)

        behavior_dual = utility_compare_dual(H_d_store, H_d_online, H_d_showrooming, L_d_store, L_d_online,
                                             L_d_showrooming)

        return behavior_uniform, behavior_dual


app = GUI()
app.master.title("Numerical Experiment")
# app.master.geometry('1020x800')
mainloop()
