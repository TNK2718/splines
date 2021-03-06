import os
from unittest import result

import bisect

from matplotlib import animation
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from absl import app
from absl import flags

import math

import numpy as np
import mpl_toolkits.mplot3d as p3d

# Implementation of Peter Comninos, "An interpolatingpiecewise bicubic surface with shapeparameters", 2001
class kochanek_bartels_surface():
    def __init__(self, grid_size_i, grid_size_j, grid_i, grid_j, grid_flag):
        self.gridsize_i = grid_size_i
        self.gridsize_j = grid_size_j
        self.grid_i = grid_i
        self.grid_j = grid_j
        self.grid_flag = grid_flag
        self.DU = []
        self.SU = []
        self.DV = []
        self.SV = []
        self.T = []

    # Compute and store vectors(U, V, T) for computing
    def update_geometry(self, world_pos, tU, cU, bU, tV, cV, bV):
        self.calc_UV(world_pos, tU, cU, bU, tV, cV, bV)
        self.calc_T(world_pos)

    # Compute position at (u, v) of the surface 
    def get_value(self, u, v, world_pos):
        # TODO: indetify (i, j) with binary search
        i = bisect.bisect_right(self.grid_i, u) - 1
        j = bisect.bisect_right(self.grid_j, v) - 1
        # print((i, j))

        # For normalize
        hU = self.grid_i[i + 1] - self.grid_i[i]
        hV = self.grid_j[j + 1] - self.grid_j[j]

        #
        if i < 0 or i >= self.gridsize_i or j < 0 or j >= self.gridsize_j:
            print("(u, v) out of range!!")
            return None

        #
        G = [[world_pos[self.get_index(i, j)], world_pos[self.get_index(i, j + 1)], self.DV[self.get_index(i, j)], self.SV[self.get_index(i, j + 1)]],
             [world_pos[self.get_index(i + 1, j)], world_pos[self.get_index(i + 1, j + 1)],
              self.DV[self.get_index(i + 1, j)], self.SV[self.get_index(i + 1, j + 1)]],
             [self.DU[self.get_index(i, j)], self.DU[self.get_index(
                 i, j + 1)], self.T[self.get_index(i, j)], self.T[self.get_index(i, j + 1)]],
             [self.SU[self.get_index(i + 1, j)], self.SU[self.get_index(i + 1, j + 1)],
              self.T[self.get_index(i + 1, j)], self.T[self.get_index(i + 1, j + 1)]]
             ]

        result = np.zeros((3))
        for k in range(3):
            for l in range(3):
                # result += self.cubic_hermite(k, u - self.grid_i[i]) * self.cubic_hermite(
                #     l, v - self.grid_j[j]) * G[k][l]
                result += self.cubic_hermite(k, (u - self.grid_i[i]) / hU) * self.cubic_hermite(
                    l, (v - self.grid_j[j]) / hV) * G[k][l]

        return result

    def calc_T(self, world_pos):
        #
        self.T = []
        #
        for j in range(self.gridsize_j):
            for i in range(self.gridsize_i):
                if self.grid_flag[self.get_index(i, j)]:
                    self.T.append(self.getT(i, j, world_pos))
                else:
                    # TODO: null
                    self.T.append(np.zeros((3)))

    def calc_UV(self, world_pos, tU, cU, bU, tV, cV, bV):
        #
        self.DU = []
        self.SU = []
        self.DV = []
        self.SV = []
        #
        for j in range(self.gridsize_j):
            for i in range(self.gridsize_i):
                if self.grid_flag[self.get_index(i, j)]:
                    self.DU.append(self.getDU(i, j, world_pos, tU, cU, bU))
                    self.SU.append(self.getSU(i, j, world_pos, tU, cU, bU))
                    self.DV.append(self.getDV(i, j, world_pos, tV, cV, bV))
                    self.SV.append(self.getSV(i, j, world_pos, tV, cV, bV))
                else:
                    # TODO: null
                    self.DU.append(np.zeros((3)))
                    self.SU.append(np.zeros((3)))
                    self.DV.append(np.zeros((3)))
                    self.SV.append(np.zeros((3)))

    def getDU(self, i, j, world_pos, tU, cU, bU):
        h = 2.0
        # if i - 1 < 0:
        #     h = self.grid_i[i + 1] - self.grid_i[i]
        # elif i + 1 >= self.gridsize_i:
        #     h = self.grid_i[i] - self.grid_i[i - 1]
        # else:
        #     h = self.grid_i[i + 1] - self.grid_i[i - 1]

        result = (1 - tU[self.get_index(i, j)]) * (1 + cU[self.get_index(i, j)]) * (1 + bU[self.get_index(i, j)]) / h * (world_pos[self.get_index(i, j)] - world_pos[self.get_index(i - 1, j)]) + \
            (1 - tU[self.get_index(i, j)]) * (1 - cU[self.get_index(i, j)]) * \
            (1 - bU[self.get_index(i, j)]) / h * \
            (world_pos[self.get_index(i + 1, j)] -
             world_pos[self.get_index(i, j)])
        return result

    def getSU(self, i, j, world_pos, tU, cU, bU):
        h = 2.0
        # if i - 1 < 0:
        #     h = self.grid_i[i + 1] - self.grid_i[i]
        # elif i + 1 >= self.gridsize_i:
        #     h = self.grid_i[i] - self.grid_i[i - 1]
        # else:
        #     h = self.grid_i[i + 1] - self.grid_i[i - 1]

        result = (1 - tU[self.get_index(i, j)]) * (1 - cU[self.get_index(i, j)]) * (1 + bU[self.get_index(i, j)]) / h * (world_pos[self.get_index(i, j)] - world_pos[self.get_index(i - 1, j)]) + \
            (1 - tU[self.get_index(i, j)]) * (1 + cU[self.get_index(i, j)]) * \
            (1 - bU[self.get_index(i, j)]) / h * \
            (world_pos[self.get_index(i + 1, j)] -
             world_pos[self.get_index(i, j)])
        return result

    def getDV(self, i, j, world_pos, tV, cV, bV):
        h = 2.0
        # if j - 1 < 0:
        #     h = self.grid_j[j + 1] - self.grid_j[j]
        # elif j + 1 >= self.gridsize_j:
        #     h = self.grid_j[j] - self.grid_j[j - 1]
        # else:
        #     h = self.grid_j[j + 1] - self.grid_j[j - 1]

        result = (1 - tV[self.get_index(i, j)]) * (1 + cV[self.get_index(i, j)]) * (1 + bV[self.get_index(i, j)]) / h * (world_pos[self.get_index(i, j)] - world_pos[self.get_index(i, j - 1)]) + \
            (1 - tV[self.get_index(i, j)]) * (1 - cV[self.get_index(i, j)]) * \
            (1 - bV[self.get_index(i, j)]) / h * \
            (world_pos[self.get_index(i, j + 1)] -
             world_pos[self.get_index(i, j)])
        return result

    def getSV(self, i, j, world_pos, tV, cV, bV):
        h = 2.0
        # if j - 1 < 0:
        #     h = self.grid_j[j + 1] - self.grid_j[j]
        # elif j + 1 >= self.gridsize_j:
        #     h = self.grid_j[j] - self.grid_j[j - 1]
        # else:
        #     h = self.grid_j[j + 1] - self.grid_j[j - 1]

        result = (1 - tV[self.get_index(i, j)]) * (1 - cV[self.get_index(i, j)]) * (1 + bV[self.get_index(i, j)]) / h * (world_pos[self.get_index(i, j)] - world_pos[self.get_index(i, j - 1)]) + \
            (1 - tV[self.get_index(i, j)]) * (1 + cV[self.get_index(i, j)]) * \
            (1 - bV[self.get_index(i, j)]) / h * \
            (world_pos[self.get_index(i, j + 1)] -
             world_pos[self.get_index(i, j)])
        return result

    def getT(self, i, j, world_pos):
        # TODO: hU
        hU = 0.0
        if i - 1 < 0:
            hU = 2.0 * (self.grid_i[i + 1] - self.grid_i[i])
        elif i + 1 >= self.gridsize_i:
            hU = 2.0 * (self.grid_i[i] - self.grid_i[i - 1])
        else:
            hU = self.grid_i[i + 1] - self.grid_i[i - 1]
        # TODO: hV
        hV = 0.0
        if j - 1 < 0:
            hV = 2.0 * (self.grid_j[j + 1] - self.grid_j[j])
        elif j + 1 >= self.gridsize_j:
            hV = 2.0 * (self.grid_j[j] - self.grid_j[j - 1])
        else:
            hV = self.grid_j[j + 1] - self.grid_j[j - 1]

        #
        # Peter Comninos, "An interpolatingpiecewise bicubic surface with shapeparameters", 2001 -> eq(26) CONTAINS MISPRINT
        first = (self.SV[self.get_index(i + 1, j)] -
                 self.DV[self.get_index(i - 1, j)]) / hU
        second = (self.SU[self.get_index(i, j + 1)] -
                  self.DU[self.get_index(i, j - 1)]) / hV
        third = None
        # TODO: boundary check and interpolation
        if j - 1 < 0:
            third = 2.0 * (world_pos[self.get_index(i - 1, 1)] -
                           world_pos[self.get_index(i - 1, 0)]) / hU * hV
        elif j + 1 >= self.gridsize_j:
            third = 2.0 * (world_pos[self.get_index(i - 1, self.gridsize_j-1)] -
                           world_pos[self.get_index(i - 1, self.gridsize_j-1 - 1)]) / hU * hV
        else:
            third = (world_pos[self.get_index(i - 1, j + 1)] -
                     world_pos[self.get_index(i - 1, j - 1)]) / hU * hV

        fourth = None
        # TODO: boundary check and interpolation
        if j - 1 < 0:
            fourth = -2.0 * (world_pos[self.get_index(i + 1, 1)] -
                             world_pos[self.get_index(i + 1, 0)]) / hU * hV
        elif j + 1 >= self.gridsize_j:
            fourth = -2.0 * (world_pos[self.get_index(i + 1, self.gridsize_j-1)] -
                             world_pos[self.get_index(i + 1, self.gridsize_j-1 - 1)]) / hU * hV
        else:
            fourth = -(world_pos[self.get_index(i + 1, j + 1)] -
                       world_pos[self.get_index(i + 1, j - 1)]) / hU * hV

        return first + second + third + fourth

    def get_index(self, i, j):
        # boundary check
        if i < 0 or j < 0:
            if i < 0 and j < 0:
                return 0
            elif i < 0:
                return j * self.gridsize_j
            elif j < 0:
                return i
        if i >= self.gridsize_i or j >= self.gridsize_j:
            if i >= self.gridsize_i and j >= self.gridsize_j:
                return (self.gridsize_j - 1) * self.gridsize_i + (self.gridsize_i - 1)
            elif i >= self.gridsize_i:
                return j * self.gridsize_i + (self.gridsize_i - 1)
            elif j >= self.gridsize_j:
                return (self.gridsize_j - 1) * self.gridsize_i + i
        #
        return j * self.gridsize_i + i

    def cubic_hermite(self, k, u):
        if k == 0:
            return 2 * u**3 - 3 * u**2 + 1
        elif k == 1:
            return -2 * u**3 + 3 * u**2
        elif k == 2:
            return u**3 - 2 * u**2 + u
        elif k == 1:
            return u**3 - u**2


def main(unused_argv):
    gridsize_i = 5
    gridsize_j = 5
    # Add epsilon to max value to draw the boundary
    grid_i = np.linspace(0.0, 10.0 + 0.000001, gridsize_i)
    grid_j = np.linspace(0.0, 10.0 + 0.000001, gridsize_j)

    grid_flag = []
    for i in range(gridsize_i * gridsize_j):
        grid_flag.append(True)

    world_pos = []
    tU, cU, bU, tV, cV, bV = [], [], [], [], [], []
    for j in range(gridsize_j):
        for i in range(gridsize_i):
            # world_pos.append(np.array(
            #     [10.0 * i / (gridsize_i - 1), 10.0 * j / (gridsize_j - 1), np.sin(2 * pi * i / (gridsize_i - 1)) * np.cos(2 * pi * j / (gridsize_j - 1))]))
            world_pos.append(np.array(
                [10.0 * i / (gridsize_i - 1), 10.0 * j / (gridsize_j - 1), math.sqrt((5.0 * (i / (gridsize_i - 1) - 0.5))**2 + (5.0 * (j / (gridsize_j - 1) - 0.5))**2)]))            # world_pos.append(np.array(
            #     [10.0 * i / (gridsize_i - 1), 5.0 * j / (gridsize_j - 1), np.sin(2 * pi * i / gridsize_i)]))
            # world_pos.append(np.array(
            #     [10.0 * i / (gridsize_i - 1), 10.0 * j / (gridsize_j - 1), 1.0 * i**2 / gridsize_i]))

            tU.append(0.0)
            cU.append(0.0)
            bU.append(0.0)
            tV.append(0.0)
            cV.append(0.0)
            bV.append(0.0)
    
    surface_gen = kochanek_bartels_surface(
        gridsize_i, gridsize_j, grid_i, grid_j, grid_flag)
    surface_gen.update_geometry(world_pos, tU, cU, bU, tV, cV, bV)

    size = 50
    pu = np.linspace(0.0, 10.0, size)
    pv = np.linspace(0.0, 10.0, size)
    x = np.zeros((size * size, 3))
    for j in range(len(pv)):
        for i in range(len(pu)):
            x[j * len(pu) + i] = surface_gen.get_value(pu[i], pv[j], world_pos)

    fig = plt.figure("Kochanek_bartels", figsize=(6, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter3D(np.array(world_pos)[:, 0], np.array(world_pos)[
                 :, 1], np.array(world_pos)[:, 2], c="red", s=70)
    ax.scatter3D(x[:, 0], x[:, 1], x[:, 2])
    ax.set_xlabel("u")
    ax.set_ylabel("v")
    ax.set_zlabel("z")
    ax.set_title("Kochanek_bartels")
    plt.show()


if __name__ == '__main__':
    app.run(main)
