#!/usr/bin/python
import numpy as np
import math
import sys
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

#############################################
#TraPPE-UA Force Field
#Temperature T=330K
global T, k, beta, kcal
T=330
k=1.3806e-23
beta=-1.0/(k*T)
kcal=6.59e-21
#############################################

def frange(x,y,j):
        while x<=y:
                yield x
                x+=j

def Energy_Angle (x):
        K=62.1
        theta0=114*math.pi/180.0
        E=K*(x-theta0)**2
        return E

def Energy_Bond (r):
        K=120.0
        r0=1.40
        E=K*(r-r0)**2
        return E

def Energy_Dihedral (x):
        K1= 1.34012
        K2=-0.271016
        K3= 3.1653
        K4= 0.0
        E=0.5*(K1*(1+math.cos(x))+K2*(1-math.cos(2*x))+K3*(1+math.cos(3*x))+K4*(1-math.cos(4*x)))
        return E

def Probability_Bond (x,y,j):
        P=[]
        for i in frange(x,y,j):
                P.append([i,math.exp(beta*Energy_Bond(i)*kcal)])
	P=np.array(P)
        Sum=sum(P[:,1])*1.0
        P=[[i, j/Sum] for i, j in P]
        return np.array(P)

def Probability_Angle(x,y,j):
        P=[]
        Sum=0
        for i in frange(x,y,j):
                P.append([i, math.exp(beta*Energy_Angle(i*math.pi/180.0)*kcal)])
	P=np.array(P)
        Sum=sum(P[:,1])*1.0
        P=[[i, j/Sum] for i, j in P]
        return np.array(P)

def Probability_Dihedral(x,y,j):
        P=[]
        for i in frange(x,y,j):
                P.append([i,math.exp(beta*Energy_Dihedral(i*math.pi/180.0)*kcal)])
	P=np.array(P)
        Sum=sum(P[:,1])*1.0
        P=[[i, j/Sum] for i, j in P]
        return np.array(P)

def LenVector(v):
	return sum([i**2 for i in v])**0.5

def normalVector(v1,v2):
        n, tmp = np.cross(v1,v2), 1
       	for i in n: tmp*=i 
        if tmp > 0:return 1
        return 0
def Periodic_Boundary(a,b,Boundary):
	v=[]
	for i in range(3):
		tmp=a[i]-b[i]+math.floor(b[i]/Boundary[i])*Boundary[i]
		v.append(tmp-round(tmp/Boundary[i])*Boundary[i])
	return v
		
def Torsion(v1,v2,v):
	v1=np.cross(v,v1)
	v2=np.cross(v,v2)	
	tmp=math.acos(np.dot(v1,v2)/(LenVector(v1)*LenVector(v2)))*180/math.pi
	if normalVector(v1,v2):	return tmp
	return 360-tmp

def Dihedral_Distribution(bond,coordination,Boundary):
	dihedral=[]
	for index in bond:
		v1=Periodic_Boundary(coordination[index[1]],coordination[index[0]],Boundary)
		v=Periodic_Boundary(coordination[index[2]],coordination[index[1]],Boundary)
		v2=Periodic_Boundary(coordination[index[2]],coordination[index[3]],Boundary)
		dihedral.append(Torsion(v1,v2,v))
	return dihedral

def Bond_Distribution(bond,coordination,Boundary):
	Lenbond=[]
	for index in bond:
		Lenbond.append(LenVector(Periodic_Boundary(coordination[index[1]],coordination[index[0]],Boundary)))
	return Lenbond

def Angle_Distribution(bond,coordination,Boundary):
	angle=[]
	for index in bond:	
		v1=Periodic_Boundary(coordination[index[0]],coordination[index[1]],Boundary)
		v2=Periodic_Boundary(coordination[index[2]],coordination[index[1]],Boundary)
		angle.append((math.acos(np.dot(v1,v2)/(LenVector(v1)*LenVector(v2)))*180/math.pi))
	return angle

def Probability(items):
	pre, value, count, prob = items[0], [], 0, []
	value.append(pre)
	for i in items:
		if pre == i: count+=1
		else:
			prob.append([pre,count])
			pre=i
			count=1
	prob=np.array(prob)
	return np.array([[i,j/sum(prob[:,1])] for i, j in prob])

def ReadBond(f):
	lines=np.array([l for l in f])
        lines=[l.split() for l in lines]
	bond, angle, dihedral = [], [], []
	global Natom
	for i in range(len(lines)):
		if len(lines[i])==2 and lines[i][1] == "atoms":
			Natom=int(lines[i][0])
			break
	for i in range(len(lines)):
		if len(lines[i])==1 and lines[i][0] == "Bonds":
			index=i+2
			break
	for i in range(index,len(lines)):
		if len(lines[i])==0:break
		bond.append([int(lines[i][2]),int(lines[i][3])])

	for i in range(len(lines)):
		if len(lines[i])==1 and lines[i][0] == "Angles":
			index=i+2
			break
	for i in range(index,len(lines)):
		if len(lines[i])==0:break
		angle.append([int(lines[i][2]),int(lines[i][3]), int(lines[i][4])])
		
	for i in range(len(lines)):
		if len(lines[i])==1 and lines[i][0] == "Dihedrals":
			index=i+2
			break
	for i in range(index,len(lines)):
		if len(lines[i])==0:break
		dihedral.append([int(lines[i][2]),int(lines[i][3]), int(lines[i][4]),int(lines[i][5])])
	return bond, angle, dihedral

def ReadLmpTrj(f):
	lines=np.array([l for l in f])
        lines=[l.split() for l in lines]
	coordination, allframe, Boundary =np.zeros(shape=(Natom+1,3)), [], []
	for i in range(len(lines)):
		if len(lines[i])==6 and lines[i][1]=="BOX":
			Boundary.append([float(lines[i+k+1][1])-float(lines[i+k+1][0]) for k in range(3)])
		if len(lines[i])==9: #format: ATOMS_id mol type x y z xu yu zu
			coordination[int(lines[i][0])][0]=float(lines[i][3])
			coordination[int(lines[i][0])][1]=float(lines[i][4])
			coordination[int(lines[i][0])][2]=float(lines[i][5])
                if len(lines[i])>1 and lines[i][1]=="TIMESTEP":
                        allframe.append(coordination)
	return allframe, Boundary

def Compute():
	f1=open(sys.argv[1],"r") # topology file generated from write_data command
	f2=open(sys.argv[2],"r") # lammps trjectory from dump command
	bond, angle, dihedral = ReadBond(f1)
	allframe, Boundary = ReadLmpTrj(f2)
	Bond, Angle, Dihedral = [], [], []
	for i in range(len(allframe)):
		Bond.append(Bond_Distribution(bond, allframe[i], Boundary[i]))
		Angle.append(Angle_Distribution(angle, allframe[i], Boundary[i]))
		Dihedral.append(Dihedral_Distribution(dihedral, allframe[i], Boundary[i]))
	Bond, Angle, Dihedral = np.array(Bond), np.array(Angle), np.array(Dihedral)
	pb=Probability(sorted([round(i,2) for i in Bond.flatten()]))
	pa=Probability(sorted([round(i) for i in Angle.flatten()]))
	pd=Probability(sorted([round(i) for i in Dihedral.flatten()]))
	return pb, pa, pd

def Plot(simu,theo,type):
	x_labels=[
		"Bond Length (Angstrom)",
		"CH2-C-CH2 Angle (Degree)",
		"Torsional Angle (Degree)"]
	title_labels=[
		"TraPPE_UA Bond Length Probability Distribution",
		"TraPPE_UA Angle Probability Distribution",
		"TraPPE_UA Torsional Angle Probability Distribution"]
	scatter_labels=[
		"Bond Length Data Points",
		"Bending Angle Data Points",
		"Torsional Angle Data Points"]
	color=['g','r','y','k']
	fig=plt.figure()
	roof=max(max(simu[:,1]),max(theo[:,1]))
	ax=fig.add_subplot(111,ylim=(0,roof+roof/3))
	ax.plot(theo[:,0],theo[:,1],'k',label="Boltzmann Distribution")
	ax.scatter(simu[:,0],simu[:,1],c='k',marker=".",label=scatter_labels[type])
	ax.set_xlabel(x_labels[type])
	ax.set_ylabel("Probability")
	ax.set_title(title_labels[type])
	ax.legend()
	fig.show()
	return

if __name__=="__main__":
	pb, pa, pd = Compute()
	PB=Probability_Bond(min(pb[:,0]),max(pb[:,0]),0.01)
	PA=Probability_Angle(20.0,180.0,1)	
	PD=Probability_Dihedral(0,360,1)
	Plot(pb,PB,0)
	Plot(pa,PA,1)
	Plot(pd,PD,2)
	raw_input()
