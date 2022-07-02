import numpy as np
"""
Module for computing the inverse kinematics of the delta robot
"""


class InverseKinematics:
    def __init__(self):
        self.L = 66.8
        self.l = 96.55
        self.d12, self.d13, self.d23 = 60.62, 60.62, 60.62
        self.d78, self.d79, self.d89 = 52.68, 52.68, 52.68
        self.alp1 = 0
        sign = 1
        self.alp2 = self.alp1 + sign * 2.0 * np.pi / 3.0
        self.alp3 = self.alp1 - sign * 2.0 * np.pi / 3.0

        self.wb = np.sqrt(self.d12 ** 2.0 / (2.0 - 2.0 * np.cos(2.0 * np.pi / 3.0)))
        self.up = np.sqrt(self.d78 ** 2.0 / (2.0 - 2.0 * np.cos(2.0 * np.pi / 3.0)))

    def compute_ik(self, ee_pos):

        x = ee_pos[0]
        y = ee_pos[1]
        z = ee_pos[2]

        a1 = x - self.wb * np.cos(self.alp1) + self.up * np.cos(self.alp1)
        b1 = y - self.wb * np.sin(self.alp1) + self.up * np.sin(self.alp1)
        a2 = x - self.wb * np.cos(self.alp2) + self.up * np.cos(self.alp2)
        b2 = y - self.wb * np.sin(self.alp2) + self.up * np.sin(self.alp2)
        a3 = x - self.wb * np.cos(self.alp3) + self.up * np.cos(self.alp3)
        b3 = y - self.wb * np.sin(self.alp3) + self.up * np.sin(self.alp3)

        E1 = -2.0 * a1 * self.L * np.cos(self.alp1) - 2.0 * b1 * self.L * np.sin(self.alp1)
        F1 = 2.0 * z * self.L
        G1 = a1 ** 2 + b1 ** 2.0 + z ** 2.0 + self.L ** 2.0 - self.l ** 2.0
        E2 = -2.0 * a2 * self.L * np.cos(self.alp2) - 2.0 * b2 * self.L * np.sin(self.alp2)
        F2 = 2.0 * z * self.L
        G2 = a2 ** 2.0 + b2 ** 2.0 + z ** 2.0 + self.L ** 2.0 - self.l ** 2.0
        E3 = -2.0 * a3 * self.L * np.cos(self.alp3) - 2.0 * b3 * self.L * np.sin(self.alp3)
        F3 = 2.0 * z * self.L
        G3 = a3 ** 2.0 + b3 ** 2.0 + z ** 2.0 + self.L ** 2.0 - self.l ** 2.0

        t1p = (-F1 + np.sqrt(E1 ** 2.0 + F1 ** 2.0 - G1 ** 2.0)) / (G1 - E1)
        t1m = (-F1 - np.sqrt(E1 ** 2.0 + F1 ** 2.0 - G1 ** 2.0)) / (G1 - E1)
        t2p = (-F2 + np.sqrt(E2 ** 2.0 + F2 ** 2.0 - G2 ** 2.0)) / (G2 - E2)
        t2m = (-F2 - np.sqrt(E2 ** 2.0 + F2 ** 2.0 - G2 ** 2.0)) / (G2 - E2)
        t3p = (-F3 + np.sqrt(E3 ** 2.0 + F3 ** 2.0 - G3 ** 2.0)) / (G3 - E3)
        t3m = (-F3 - np.sqrt(E3 ** 2.0 + F3 ** 2.0 - G3 ** 2.0)) / (G3 - E3)

        angles = [2*np.rad2deg(np.arctan(t1m)), 2*np.rad2deg(np.arctan(t2m)), 2*np.rad2deg(np.arctan(t3m))]

        return angles
