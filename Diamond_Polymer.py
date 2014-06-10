#!/usr/bin/python
import numpy as np
import sys
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

#################################################################
# Endian: wangn.0919@gmail.com                                  #
#                                                               #
# This file can create a polymer of diamond-like structure with #
# periodic boundary. This polymer will have Nx, Ny, Nz units of #
# cell along the x, y and z coordinates. Each chain has M       #
# number of CH2 United Atoms.                                   #
#                                                               #
# Excuting this file: ./Diamond_Plymer.py Nx Ny Nz  M           #
#################################################################

Cell=[
        [ 0,  0,  0],
        [ 1,        0,        -1/2**0.5],
        [-1,        0,        -1/2**0.5],
        [ 0,        1,         1/2**0.5],
        [ 0,       -1,         1/2**0.5],
        [-1,        1,        -2/2**0.5],
        [-1,       -1,        -2/2**0.5],
        [ 1,       -1,         2/2**0.5],
        [-1,       -1,         2/2**0.5]
]

Tetrahedron=[[0,1,-1/2**0.5],[0,-1,-1/2**0.5],[1,0,1/2**0.5],[-1,0,1/2**0.5]]

def Ini(chains):
        global CC
        CC=1.4
        if chains%2:
                return CC/2*3**0.5*(chains+1)/1.5**0.5
        else:
                return ((CC/2*3**0.5*(chains+1))**2+(CC/2)**2)**0.5/1.5**0.5

def Network(Cell,Nx, Ny, Nz):
        global Q
        Q, x, xy, xyz = [], [], [], []
        for i in Cell:
                for j in range(N):
                        x.append(i+np.array([j,0,0])*2*Scaling)
        x=np.array(x)
        for i in x:
                for j in range(N):
                        xy.append(i+np.array([0,j,0])*2*Scaling)
        xy=np.array(xy)
        for i in xy:
                for j in range(N):
                        xyz.append(i+np.array([0,0,j])*4/2**0.5*Scaling)
        xyz=np.array(xyz)
        for i in xyz:
                if Search(i,Q):continue
                else:Q.append(i)
        return

def Search(node,T):
        for i in range(len(T)):
                k=0
                for j in range(3):
                        if node[j]-T[i][j] < 0.001 and node[j]-T[i][j]> -0.001:k=k+1
                if k==3: return i+1
        return 0

def Connection(T):
        global connected, LenT
        LenT=len(T)
        T1=np.array(Cell[1:5])
        T2=np.array(Tetrahedron)*Scaling
        connected=np.zeros(shape=(LenT+1,4))
        for i in range(LenT):
                ID=[Search(j,T) for j in T1+T[i]]
                if sum(ID)==0:
                        ID=[Search(j,T) for j in T2+T[i]]
                connected[i+1]=ID
        connected=np.array(connected)
        return

def Reset(i,k):
        for j in range(4):
                if connected[k][j] == i:
                        connected[k][j]=0
        return

def InArray(x,Array):
        for i in Array:
                if x == i:return 1
        return 0

def IndexofArray(x,Array):
        for i in range(len(Array)):
                if x == Array[i]:return i+1
        return 0

def BuildChain(begin,end,chains):
        n=chains+2
        end=np.array(end)
        begin=np.array(begin)
        normz=np.array([0,0,1])
        end2=end+CC/2*normz
        line=[]
        for i in range(n):
                if chains%2:
                        bead=begin+(end-begin)*i*1.0/(n-1)
                else:
                        bead=begin+(end2-begin)*i*1.0/(n-1)
                if i%2: bead=bead+CC/2*normz
                line.append(bead)
        return line

def MirrorEdge(edgeR,edgeL,m,T):
        for i in edgeL:
                for j in edgeR:
                        p=0
                        for k in range(3):
                                if k==m:continue
                                if T[j-1][k]-T[i-1][k] < 0.001 and T[j-1][k]-T[i-1][k]> -0.001:
                                        p=p+1
                        if p==2: mirror[m][j]=i
        return

def Periodize(T):
        T=np.array(T)
        XR=max(T[:,0])
        YR=max(T[:,1])
        ZR=max(T[:,2])
        XL=min(T[:,0])
        YL=min(T[:,1])
        ZL=min(T[:,2])
        global edge, edgeXR, edgeYR, edgeXL, edgeYL, edgeZR, edgeZL
        edgeXR, edgeYR, edgeXL, edgeYL, edgeZR, edgeZL = [], [], [], [], [], []
        for i in range(len(T)):
                if T[i][0] == XR: edgeXR.append(i+1)
                if T[i][1] == YR: edgeYR.append(i+1)
                if T[i][2] == ZR: edgeZR.append(i+1)
                if T[i][0] == XL: edgeXL.append(i+1)
                if T[i][1] == YL: edgeYL.append(i+1)
                if T[i][2] == ZL: edgeZL.append(i+1)
        edge=edgeXR+edgeYR+edgeZR
        edgeL=edgeXL+edgeYL+edgeZL
        global mirror
        mirror=np.zeros(shape=(3,max(max(edge),max(edgeL))+1))
        MirrorEdge(edgeXR,edgeXL,0,T)
        MirrorEdge(edgeYR,edgeYL,1,T)
        MirrorEdge(edgeZR,edgeZL,2,T)
        return

def CheckBoundary(a,b):
        if InArray(a,edge) and InArray(b,edge):
                return 0
        return 1

def ConnectBeads(T,chains):
        T=np.array(T)
        nbond, n = 0, LenT+1
        for i in range(1,LenT+1):
                for j in range(4):
                        if connected[i][j]:
                                Reset(i,connected[i][j])
                                nbond+=1
        bond=np.zeros(nbond*(chains+2))
        reverse=[]
        core=np.zeros(shape=(LenT+1,4))
        for i in range(1,LenT+1):
                p=0
                for j in range(4):
                        if connected[i][j] and CheckBoundary(i,connected[i][j]):
                                begin=i
                                end=connected[i][j]
                                line=BuildChain(T[begin-1],T[end-1],chains)
                                core[begin][p]=n
                                p=p+1
                                for k in range(1,len(line)-2):
                                        bond[n]=n+1
                                        n=n+1
                                        Q.append(line[k])
                                Q.append(line[len(line)-2])
                                bond[n]=end
                                reverse.append([end,n])
                                n=n+1
        return bond, core, reverse

def UpdateTopology(Q,bond,core,reverse):
        RL=[]
        for i in range(3):
                for j in range(len(mirror[0,:])):
                        if mirror[i][j]:
                                RL.append([j, mirror[i][j]])
        RL=np.array(sorted(RL))
        for i in range(len(RL)):
                for j in range(len(RL)):
                        if RL[i][1]==RL[j][0]:
                                RL[i][1]=RL[j][1]
        Coordination, updateIndex, updateRL = [], [], []
        for i in range(len(Q)):
                if InArray(i+1,RL[:,0]):continue
                else:
                        Coordination.append(Q[i])
                        updateIndex.append(i+1)
        ReIndex=np.zeros(max(updateIndex)+1)
        for i in range(len(updateIndex)):
                ReIndex[updateIndex[i]]=i+1
        for i, j in RL:
                for k in range(len(reverse)):
                        if i == reverse[k][0]:
                                updateRL.append([reverse[k][1],j])
        RL=[[ReIndex[int(i)],ReIndex[int(j)]] for i, j in updateRL]
        bond=[[ReIndex[i], ReIndex[int(bond[i])]] for i in range(len(bond)) if ReIndex[bond[i]]]
        tmp=[]
        for i in range(len(core)):
                for j in range(4):
                        if core[i][j]:
                                tmp.append([ReIndex[i], ReIndex[core[i][j]]])
        AllBond=np.array([[int(i), int(j)] for i, j in sorted(tmp+bond+RL)])
        Coordination=np.array([[round(x,6), round(y,6), round(z,6)] for x, y, z in Coordination])
        return Coordination, AllBond

def Plot(T,bond,core):
        fig=plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        T=np.array(T)
        ax.scatter(T[:,0],T[:,1],T[:,2])
        for i in range(1,LenT+1):
                for j in range(4):
                        if core[i][j]:
                                b=np.array([T[i-1],T[core[i][j]-1]])
                                ax.plot(b[:,0],b[:,1],b[:,2],c='k')
        for i in range(len(bond)):
                if bond[i]:
                        b=np.array([T[i-1],T[bond[i]-1]])
                        ax.plot(b[:,0],b[:,1],b[:,2],c='k')
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        fig.show()
        return

def Plot2(T,bond):
        fig=plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        T=np.array(T)
        ax.scatter(T[:,0],T[:,1],T[:,2])
        for i in range(len(bond)):
                b=np.array([T[bond[i][0]-1],T[bond[i][1]-1]])
                ax.plot(b[:,0],b[:,1],b[:,2],c='k')
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.set_title("Network of Periodic Boundary")
        fig.show()
        return

def PrintData(Coordination, Bond):
        f1=open('Coordination.data', 'wr')
        for i in Coordination: print >>f1, i
        f1.close()
        f2=open('Bond.data', 'wr')
        for i in Bond: print >>f2, i
        f2.close()
        return

if __name__ == "__main__":
        Cell=np.array(Cell)
        Nx=int(sys.argv[1])
        Ny=int(sys.argv[2])
        Nz=int(sys.argv[3])
        chains=int(sys.argv[4])
        Scaling=Ini(chains)
        Cell=Scaling*Cell
        Network(Cell,Nx,Ny,Nz)
        Connection(Q)
        Periodize(Q)
        bond, core, reverse = ConnectBeads(Q,chains)
        Plot(Q, bond, core)
        Coordination, AllBond = UpdateTopology(Q,bond,core,reverse)
        Plot2(Coordination, AllBond)
        PrintData(Coordination, AllBond)
        raw_input()
