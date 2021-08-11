# -*- coding: utf-8 -*-
"""
@Created at 2021/4/27 15:00
@Author: Kurt
@file:main.py
@Desc:
"""


def utility_compare(Hstore, Honline, Lstore, Lonline):
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


class NumSolver:
    def __init__(self, V=0.8, GAMMA=0.2, P_HAT=0.1, DELTA_h=0.6, DELTA_l=0.2, CR=0.7):
        self.V = V
        self.GAMMA = GAMMA
        self.P_HAT = P_HAT
        self.DELTA_h = DELTA_h
        self.DELTA_l = DELTA_l
        self.CR = CR

    def get_utility(self, TYPE, p, poff, c, con):

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
        # return np.round(u_store, 3), np.round(u_online, 3)

        # quicker round
        u_store = int(u_store * 1000000 + 0.5) / 1000000
        u_online = int(u_online * 1000000 + 0.5) / 1000000
        return u_store, u_online

    def get_profit(self, p, poff, c, con, logger=False):

        if p >= self.P_HAT:
            delta = self.DELTA_h
        else:
            delta = self.DELTA_l

        H_u_store, H_u_online = self.get_utility("H", p, poff, c, con)
        L_u_store, L_u_online = self.get_utility("L", p, poff, c, con)

        if logger:
            if max(L_u_store, L_u_online) >= 0 and max(H_u_store, H_u_online) >= 0:
                if L_u_store > L_u_online:
                    print("sell to both H and L segments offline")
                else:
                    print("sell to both H and L segments online")

            elif max(L_u_store, L_u_online) >= 0 and max(H_u_store, H_u_online) < 0:
                if L_u_store > L_u_online:
                    print("sell to L segment offline")
                else:
                    print("sell to L segment online")

            elif max(L_u_store, L_u_online) < 0 and max(H_u_store, H_u_online) >= 0:
                if H_u_store > H_u_online:
                    print("sell to H segment offline")
                else:
                    print("sell to H segment online")

            elif max(L_u_store, L_u_online) < 0 and max(H_u_store, H_u_online) < 0:
                print("No sales")
            else:
                raise Exception("other cases")

        pi_H_online = 1 / 4 * self.GAMMA * (2 * p - delta * (p + self.CR)) * (
            [1 if (H_u_online >= 0 and H_u_online >= H_u_store) else 0][0])
        pi_H_store = 1 / 4 * self.GAMMA * poff * ([1 if (H_u_store >= 0 and H_u_online < H_u_store) else 0][0])

        pi_L_online = 1 / 4 * (1 - self.GAMMA) * (2 * p - delta * (p + self.CR)) * (
            [1 if (L_u_online >= 0 and L_u_online >= L_u_store) else 0][0])
        pi_L_store = 1 / 4 * (1 - self.GAMMA) * poff * ([1 if (L_u_store >= 0 and L_u_online < L_u_store) else 0][0])

        expected_profit = pi_H_online + pi_H_store + pi_L_online + pi_L_store

        return expected_profit
