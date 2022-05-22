import os

from matplotlib import animation
import matplotlib.pyplot as plt

from absl import app
from absl import flags

import math

import numpy as np
import mpl_toolkits.mplot3d as p3d


def kochanek_bartels_surface(u, v, inputs, grid):
    # TODO: calculate i, j from grid
    i = 1
    j = 1

    # get world positions
    P = np.zeros((4, 4))
    for indexU in range(4):
        for indexV in range(4):
            P[indexU, indexV] = inputs['world_pos'][point_accessor(
                i + indexU - 1, j + indexV - 1, 10, 10)]

    # P00 = inputs['world_pos'][point_accessor(i, j, 10, 10)]

    # TODO: parse t, b, c
    tU = np.ones((2, 2))
    bU = np.ones((2, 2))
    cU = np.ones((2, 2))

    tV = np.ones((2, 2))
    bV = np.ones((2, 2))
    cV = np.ones((2, 2))

    U = np.array([u**3, u**2, u, 1])
    V = np.array([v**3, v**2, v, 1])
    M_H = np.matrix([2, -2, 1, 1], [-3, 3, -2, -1], [0, 0, 1, 0], [1, 0, 0, 0])

    #
    DV00 = (1 - tV[0, 0]) * (1 + cV[0, 0]) * (1 + bV[0, 0]) / 2.0 * (P(1, 1) - P(1, 0)) + \
        (1 - tV[0, 0]) * (1 - cV[0, 0]) * \
        (1 - bV[0, 0]) / 2.0 * (P(1, 2) - P(1, 1))
    # 
    SV01 = (1 - tV[0, 0]) * (1 + cV[0, 0]) * (1 + bV[0, 0]) / 2.0 * (P(1, 1) - P(1, 0)) + \
        (1 - tV[0, 0]) * (1 - cV[0, 0]) * \
        (1 - bV[0, 0]) / 2.0 * (P(1, 2) - P(1, 1))

    G = np.matrix([])


def point_accessor(i, j, m, n):
    return n * i + j


def main(unused_argv):
    print()


if __name__ == '__main__':
    app.run(main)
