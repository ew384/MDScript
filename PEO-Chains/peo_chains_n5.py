#!/usr/bin/python
import sys
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt


def Ini(chains):
    global CC
    CC = 1.40
    if chains % 2:
        return CC * (chains + 1) * 2 / 3
    else:
        return CC / 3 * (4 * (chains + 1) ** 2 + 2) ** 0.5


def BuildChain(begin, end, chains):
    n = chains + 2 
    end = np.array(end)
    begin = np.array(begin)
    normz = np.array([0, 0, 1])
    end2 = end + CC / 2 * normz
    line = []
    for i in range(n):
        if chains % 2:
            bead = begin + (end - begin) * i * 1.0 / (n - 1)
        else:
            bead = begin + (end2 - begin) * i * 1.0 / (n - 1)
        if i % 2: bead = bead + CC / 2 * normz
        line.append(bead)
    return line


def PEO(chains):
    oxygen = [1 + i * 3 for i in range(chains + 1)]
    chains = chains * 3 + 1
    Scaling = Ini(chains)
    begin = [0, 0, 0]
    end = [0, Scaling, 0]
    peo = BuildChain(begin, end, chains)
    bond = [[i, i + 1] for i in range(1, len(peo))]
    return np.array(peo), oxygen, bond


def PrintLT(xyz, oxygen, bond):
    f = open('PEO.lt', 'wr')
    f.write("import \"trappe_modified.lt\"\n\nPEO{\n    write('Data Atoms') {\n")
    for i in range(len(xyz)):
        if i == 0 or i == len(xyz) - 1:
            f.write('{}{}  {}  {}  {}  {}\n'.format("    $atom:a", i + 1, " $mol:. @atom:TraPPE/CH3   0.25", xyz[i][0],
                                                    xyz[i][1], xyz[i][2]))
        elif i in oxygen:
            f.write('{}{}  {}  {}  {}  {}\n'.format("    $atom:a", i + 1, " $mol:. @atom:TraPPE/O   -0.5", xyz[i][0],
                                                    xyz[i][1], xyz[i][2]))
        else:
            f.write('{}{}  {}  {}  {}  {}\n'.format("    $atom:a", i + 1, " $mol:. @atom:TraPPE/CH2   0.25", xyz[i][0],
                                                    xyz[i][1], xyz[i][2]))
    f.write("  }\n\n  write('Data Bonds') {\n")
    for i in range(len(bond)):
        if bond[i][0] == 1 or bond[i][1] == len(xyz):
            f.write('{}{}{}{}{}{}\n'.format("    $bond:b", i + 1, " @bond:TraPPE/OCH3  $atom:a", bond[i][0], " $atom:a",
                                            bond[i][1]))
        elif (bond[i][0]-1) in oxygen or (bond[i][1]-1) in oxygen:
            f.write('{}{}{}{}{}{}\n'.format("    $bond:b", i + 1, " @bond:TraPPE/OCH2  $atom:a", bond[i][0], " $atom:a",
                                            bond[i][1]))
        else:
            f.write('{}{}{}{}{}{}\n'.format("    $bond:b", i + 1, " @bond:TraPPE/CH2CH2  $atom:a", bond[i][0],
                                            " $atom:a", bond[i][1]))

    f.write("  }\n}")
    f.close()
    return


def PrintSystemLT(xhi, yhi, zhi, xlo, ylo, zlo,length):
    f = open('system.lt', 'wr')
    f.write("import \"trappe_modified.lt\"\nimport \"PEO.lt\"\nwrite_once(\"Data Boundary\") {\n ")
    f.write('{} {} {} {}\n'.format(xlo, xhi, "xlo", "xhi"))
    f.write('{} {} {} {}\n'.format(ylo, yhi, "ylo", "yhi"))
    f.write('{} {} {} {}\n'.format(zlo, zhi, "zlo", "zhi"))
    f.write("  }\n")
    f.write('{}{}{}\n'.format("peo = new PEO [", 15, "].move(10,0,0)"))
    f.write('{}{}{}{}\n'.format("		   [2", "].move(0,",length+4,",0)"))
    f.write('{}{}{}\n'.format("		   [", 5, "].move(0,0,10)"))
    f.close()
    return


if __name__ == "__main__":
    chains = int(sys.argv[1])
    Coordination, oxygen, bond = PEO(chains)
    xhi = max(Coordination[:, 0])*6+60 + 0.1
    yhi = max(Coordination[:, 1])*5+5*40 + 0.1
    zhi = max(Coordination[:, 2])*5+50 + 0.1
    xlo = min(Coordination[:, 0]) - 0.1
    ylo = min(Coordination[:, 1]) - 0.1
    zlo = min(Coordination[:, 2]) - 0.1
    length = max(Coordination[:, 1])-min(Coordination[:, 1])
    #PrintSystemLT(xhi, yhi, zhi, xlo, ylo, zlo,length)
    PrintLT(Coordination, oxygen, bond)
