#!/usr/bin/python
import numpy as np
import operator
import math
import sys
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


def LenVector(v):
    return sum([i ** 2 for i in v])


def Periodic_Boundary(a, b, Boundary):
    v = []
    for i in range(3):
        tmp = a[i] - b[i] + math.floor(b[i] / Boundary[i]) * Boundary[i]
        v.append(tmp - round(tmp / Boundary[i]) * Boundary[i])
    return v


def Probability(items):
    pre, value, count, prob = items[0], [], 0, []
    value.append(pre)
    for i in items:
        if pre == i:
            count += 1
        else:
            prob.append([pre, count])
            pre = i
            count = 1
    prob = np.array(prob)
    return np.array([[i, j / sum(prob[:, 1])] for i, j in prob])


def Compute_R2(CH3, Boundary):
    R2 = []
    for i in range(0, len(CH3), 2):
        a = CH3[i][1][1:4]
        b = CH3[i+1][1][1:4]
        R2.append(LenVector(Periodic_Boundary(a, b, Boundary)))
    if len(R2): return sum(R2)/len(R2)
    else: return 0


def ReadLmpTrj(f):
    lines = np.array([l for l in f])
    lines = [l.split() for l in lines]
    allframe_R2, Boundary = [], []
    CH3 = {}
    for i in range(len(lines)):
        if len(lines[i]) == 6 and lines[i][1] == "BOX":
            Boundary=[float(lines[i + k + 1][1]) - float(lines[i + k + 1][0]) for k in range(3)]
        if len(lines[i]) == 9 and int(lines[i][2]) == 4:  # format: ATOMS_id mol type x y z xu yu zu
            CH3[lines[i][0]] = [int(lines[i][1]), float(lines[i][3]), float(lines[i][4]), float(lines[i][5])]
        if len(lines[i]) > 1 and lines[i][1] == "TIMESTEP":
            CH3 = sorted(CH3.items(), key=operator.itemgetter(1))
            R2 = Compute_R2(CH3, Boundary)
            if R2: allframe_R2.append(R2)
            CH3, Boundary = {}, []
    return allframe_R2


if __name__ == "__main__":
    file=open(sys.argv[1], "r")
    R2=ReadLmpTrj(file)
    print round(sum(R2)/len(R2)/100,4),"nm"
