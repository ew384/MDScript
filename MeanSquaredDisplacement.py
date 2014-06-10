#!/usr/bin/python
import numpy as np
import math
import sys
############################################
# United Atoms Model:                      #    
# C, CH, CH2, CH3                          #
mass=[0,12.01,13.018,14.17,15.25]
# LAMMPS dump file format:                 #
#                                          #
# ITEM: TIMESTEP                           #
# "TimeBegin"                              #
# ITEM: NUMBER OF ATOMS                    #
# "NumberOfAtoms"                          #
# ITEM: BOX BOUNDS pp pp pp                #
# "xlo xhi"                                #
# "ylo yhi"                                #
# "zlo zhi"                                #
# ITEM: ATOMS id mol type x y z xu yu zu   #
#                                          #
############################################

def COM(index,P):
        x, y, z, masstotal = 0, 0, 0, 0
        index=index*Natom
        for i in range(Natom):
                masstotal+=mass[P[index+i][0]]
                x=mass[P[index+i][0]]*P[index+i][1]
                y=mass[P[index+i][0]]*P[index+i][2]
                z=mass[P[index+i][0]]*P[index+i][3]
        x/=masstotal
        y/=masstotal
        z/=masstotal
        com=[x,y,z]
        return np.array(com)

def Distance(P):
        d=np.zeros(shape=(Ntime+1,Ntime+1))
        for i in range(Ntime):
                for j in range(Ntime):
                                if i<j:
                                        tmp=COM(j,P)-COM(i,P)
                                        d[i][j]=sum([k**2 for k in tmp])**0.5
                                        print d[i][j]
        return d

def MSD(d):
        msd=[]
        for dt in range(Ntime):
                tmp, N = 0, 0
                for i in range(Ntime):
                        if i+dt<=Ntime:
                                tmp+=d[i][i+dt]
                                N+=1
                        msd.append(tmp/N)
        return msd

def Read(f):
        global Natom, Ntime
        Natom, Ntime = 0, 0
        lines=np.array([l for l in f])
        lines=[l.split() for l in lines]
        Natom=int(lines[3][0])
        Positions=[]
        for i in range(1,len(lines)):
                if len(lines[i])>1 and lines[i][1]=="TIMESTEP":
                        timestep=int(lines[i+1][0])-int(lines[1][0])
                        break
        for i in range(len(lines)):
                if len(lines[i])==9:
                        Positions.append([int(lines[i][2]),float(lines[i][6]),float(lines[i][7]),float(lines[i][8])])
                if len(lines[i])>1 and lines[i][1]=="TIMESTEP":
                        Ntime+=1
        time=[i*timestep for i in range(Ntime)]
        msd=MSD(Distance(Positions))
        return time,msd

def Plot(time,msd):
        fig=plt.figure()
        ax=fig.add_subplot(111)
        ax.set_xlabel("Time (fs)")
        ax.set_ylabel("MSD (Angstrom^2)")
        ax.set_title("Mean Squared Displacement")
        ax.plot(time,msd)
        fig.show()
        raw_input()

if __name__ == "__main__":
        f=open(sys.argv[1],'r')
        f2=open(sys.argv[2],'wr')
        time,msd=Read(f)
        f.close()
        for i in range(len(time)):
                print >>f2,time[i],msd[i]
        f2.close()
