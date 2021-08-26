# -*- coding: utf-8 -*-
"""
@Created at 2021/1/15 21:11
@Author: Kurt
@file:GUI.py
@Desc:
Basic discrete model without consideration of price-related return and showrooming-free consumers.
"""
from tkinter import *
from main import utility_compare, NumSolver
# from main_asymmetric import utility_compare, NumSolver  # main_asymmetric is the forgetting return consumers case.
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import numpy as np

matplotlib.use('TkAgg')


def find_result(DataTuple):
    optimal_profits_u, optimal_prices_u, optimal_profits_d, optimal_prices_d, behavior_tuples = DataTuple
    flag = False
    Res = {"c": [], "online_price_uniform": [], "online_price_dual": [], "profit_uniform": [],
           "profit_dual": [], "behavior_u": [], "behavior_d": []}
    for c, profit_u, profit_d, price_u, price_d, behavior_tuple \
            in zip(np.arange(0, 1, 0.01), optimal_profits_u, optimal_profits_d,
                   optimal_prices_u, optimal_prices_d, behavior_tuples):
        if profit_u > profit_d and price_u < price_d:
            flag = True
            Res["c"].append(c)
            Res["online_price_uniform"].append(price_u)
            Res["online_price_dual"].append(price_d)
            Res["profit_uniform"].append(profit_u)
            Res["profit_dual"].append(profit_d)
            Res['behavior_u'].append(behavior_tuple[0])
            Res['behavior_d'].append(behavior_tuple[1])

    return flag, Res


def parse(result_dict):
    cs = result_dict['c']
    ps = result_dict["online_price_uniform"]
    pons = result_dict["online_price_dual"]
    profit_us = result_dict["profit_uniform"]
    profit_ds = result_dict["profit_dual"]
    behavior_us = result_dict['behavior_u']
    behavior_ds = result_dict['behavior_d']
    result_string = "".join(["{:.3f}, \t {:.3f}, \t {:.3f}, \t {:.3f}, \t {:.3f}, \t {}, \t    {}, \n".format(
        s1, s2, s3, s4, s5, s6, s7) for s1, s2, s3, s4, s5, s6, s7 in
        zip(cs, ps, pons, profit_us, profit_ds, behavior_us, behavior_ds)])

    final_string = "c, \t p, \t pon, \t profit_u, \t profit_d, behavior_u, behavior_d \n" + result_string
    return final_string


class GUI(Frame):
    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.pack()
        self.create_widgets()

    def create_widgets(self):

        self.HighReturnProb = Scale(self, orient=HORIZONTAL, from_=0, to=1,
                                    label="High return probability delta_H",
                                    resolution=0.01, length=250)
        self.HighReturnProb.grid(row=0, column=0)
        self.HighReturnProb.set(0.89)

        self.LowReturnProb = Scale(self, orient=HORIZONTAL, from_=0, to=1,
                                   label="Low return probability delta_L",
                                   resolution=0.01, length=250)
        self.LowReturnProb.grid(row=1, column=0)
        self.LowReturnProb.set(0.89)

        self.PriceThreshold = Scale(self, orient=HORIZONTAL, from_=0, to=1,
                                    label="Threshold of price p^hat. ",
                                    resolution=0.01, length=250)
        self.PriceThreshold.grid(row=2, column=0)
        self.PriceThreshold.set(0.1)

        self.ValLow = Scale(self, orient=HORIZONTAL, from_=0, to=1,
                            label="Valuation of low type segment V",
                            resolution=0.01, length=250)
        self.ValLow.grid(row=3, column=0)
        self.ValLow.set(0.52)

        self.PropHigh = Scale(self, orient=HORIZONTAL, from_=0, to=1,
                              label="Proportion of high type segment gamma ",
                              resolution=0.01, length=250)
        self.PropHigh.grid(row=4, column=0)
        self.PropHigh.set(0.25)

        # self.OfflineCost = Scale(self, orient=HORIZONTAL, from_=0, to=1,
        #                label="Offline cost c",
        #                resolution=0.01, length=250)
        # self.OfflineCost.pack()

        self.OnlineCost = Scale(self, orient=HORIZONTAL, from_=0, to=1,
                                label="Online cost con",
                                resolution=0.01, length=250)
        self.OnlineCost.grid(row=5, column=0)
        self.OnlineCost.set(0.14)

        self.ReturnCostRet = Scale(self, orient=HORIZONTAL, from_=0, to=1,
                                   label="The retailer's return cost cr",
                                   resolution=0.01, length=250)
        self.ReturnCostRet.grid(row=6, column=0)
        self.ReturnCostRet.set(0.11)

        self.startButton = Button(self, text='Start', command=self.read_params)
        self.startButton.grid(row=7, column=0)

        self.quitButton = Button(self, text='Quit', command=self.quit)
        self.quitButton.grid(row=8, column=0)

        self.fig = Figure(figsize=(5, 4), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.get_tk_widget().grid(row=0, column=1, columnspan=1, rowspan=4)

        self.ResLabel = Label(self, justify=LEFT)
        self.ResLabel.grid(row=5, column=1, columnspan=1, rowspan=3)

    def plot_figures(self, V, gamma, P_HAT, DELTA_h, DELTA_l, CR, con):
        self.fig.clear()
        self.solver = NumSolver(V, gamma, P_HAT, DELTA_h, DELTA_l, CR)

        optimal_profits_u = []
        optimal_prices_u = []

        optimal_profits_d = []
        optimal_prices_d = []

        behavior_tuples = []

        for current_c in np.arange(0, 1, 0.01):
            optimal_profit_u = 0
            optimal_price_u = 0
            optimal_profit_d = 0
            optimal_price_d = 0

            for p_current in np.arange(0, 1, 0.005):
                profit_current_u = self.solver.get_profit(p=p_current, poff=p_current, c=current_c, con=con)
                if profit_current_u > optimal_profit_u:
                    optimal_profit_u = profit_current_u
                    optimal_price_u = p_current

                profit_current_d = self.solver.get_profit(p=p_current, poff=p_current + con, c=current_c, con=con)
                if profit_current_d > optimal_profit_d:
                    optimal_profit_d = profit_current_d
                    optimal_price_d = p_current

            behavior_tuple = self.analyze_behavior(current_c, con, optimal_price_u,
                                                   optimal_price_d, optimal_price_d + con)

            optimal_profits_u.append(optimal_profit_u)
            optimal_prices_u.append(optimal_price_u)
            optimal_profits_d.append(optimal_profit_d)
            optimal_prices_d.append(optimal_price_d)
            behavior_tuples.append(behavior_tuple)

        for current_c, price_u, price_d, profit_u, profit_d in zip(
                np.arange(0, 1, 0.01), optimal_prices_u, optimal_prices_d, optimal_profits_u, optimal_profits_d):
            behavior_uniform, behavior_dual = self.analyze_behavior(current_c, con, price_u,
                                                                    price_d, price_d + con)
            print("current c: {:.5f}, price_u: {:.5f}, price_d: {:.5f}".format(current_c, price_u, price_d))
            print("profit_u: {:.5f}, profit_d: {:.5f}".format(profit_u, profit_d))
            print("behavior u:{}, d:{}".format(behavior_uniform, behavior_dual))

        ax1 = self.fig.add_subplot(2, 1, 1)
        ax1.plot(np.arange(0, 1, 0.01), optimal_profits_u, 'g', label="Uniform")
        ax1.plot(np.arange(0, 1, 0.01), optimal_profits_d, 'k--', label="Dual")

        ax2 = self.fig.add_subplot(2, 1, 2)
        ax2.plot(np.arange(0, 1, 0.01), optimal_prices_u, 'g', label="Uniform")
        ax2.plot(np.arange(0, 1, 0.01), optimal_prices_d, 'k--', label="Online of Dual")
        ax2.plot(np.arange(0, 1, 0.01), [x + con for x in optimal_prices_d], linestyle='-.', label="Offline of Dual")

        ax1.set_title("profits")
        ax2.set_title("prices")

        ax1.legend(prop=dict(size=6), frameon=False)
        ax2.legend(prop=dict(size=6), frameon=False)
        self.fig.tight_layout()
        self.canvas.draw()

        return optimal_profits_u, optimal_prices_u, optimal_profits_d, optimal_prices_d, behavior_tuples

    def read_params(self):
        self.fig.clear()

        delta_h = float(self.HighReturnProb.get())
        delta_l = float(self.LowReturnProb.get())
        phat = float(self.PriceThreshold.get())
        V = float(self.ValLow.get())
        gamma = float(self.PropHigh.get())
        # c = float(self.OfflineCost.get())
        con = float(self.OnlineCost.get())
        cr = float(self.ReturnCostRet.get())

        data = self.plot_figures(V=V, gamma=gamma, P_HAT=phat, DELTA_h=delta_h, DELTA_l=delta_l, CR=cr, con=con)

        result = find_result(data)

        if not result[0]:
            self.ResLabel.config(text='Is there a desirable uniform price solution? \n' + " No.")
        else:
            self.ResLabel.config(
                text='Is there a desirable uniform price solution?\n'
                     'Yes, e.g.: \n \n {} \n'
                     '**Note: BO: both segments buy online; '
                     'BS: both segments buy offline;\n'
                     'HO: only H type segment buy online; HS: only H type segment buy offline.'.format(
                    parse(result[1])))

    def analyze_behavior(self, c, con, uniform_price, dual_online_price, dual_offline_price):
        H_u_store, H_u_online = self.solver.get_utility("H", uniform_price, uniform_price, c, con=con)
        L_u_store, L_u_online = self.solver.get_utility("L", uniform_price, uniform_price, c, con=con)

        H_d_store, H_d_online = self.solver.get_utility("H", dual_online_price, dual_offline_price, c, con=con)
        L_d_store, L_d_online = self.solver.get_utility("L", dual_online_price, dual_offline_price, c, con=con)

        behavior_uniform = utility_compare(H_u_store, H_u_online, L_u_store, L_u_online)
        behavior_dual = utility_compare(H_d_store, H_d_online, L_d_store, L_d_online)
        return behavior_uniform, behavior_dual


app = GUI()
app.master.title("Numerical Experiments")
# app.master.geometry('1020x800')
mainloop()
