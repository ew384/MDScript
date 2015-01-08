#!/usr/bin/python
import numpy as np
import math
import sys
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
def is_number(s):
	try:
		float(s)
		return True
	except ValueError:
		pass
	return False

def Compute_Density(f):
	mass_CH3O_CH2CH2O_9_CH3 = 15.2507*2+16*10+2*9*14.1707  
	mass_CH3O_CH2CH2O_5_CH3 = 15.2507*2+16*10+2*5*14.1707  
	mass9 =mass_CH3O_CH2CH2O_9_CH3 * 150
	mass5 =mass_CH3O_CH2CH2O_5_CH3 * 50
	mass = mass5
	myline=[]
	time=0
	for line in f:
		lines=np.array(line.split())
                if len(lines)==0:continue
                if lines[0]=='thermo' and len(lines)==2:
                        thermo=int(lines[1])
		if is_number(lines[0]) and len(lines)==1:
			print lines[0]
                        myline.append([time*thermo,mass*10/6.02/float(lines[0])])
			time+=1
        return np.array(myline)

def loop(n):
        for i in range(n):
                f=open(sys.argv[i+2],'r')
                line=Compute_Density(f)
		last_10_items=line[-10:,1]
		average_density=sum(last_10_items)/len(last_10_items)
		max_density=max(last_10_items)
		print round(average_density,4),"g/cm3"
                ax.plot(line[:,0],line[:,1],color=c[i%6])
        return

def MultiPlot():
        global ax, m, c
        fig=plt.figure()
        ax=fig.add_subplot(111,ylim=(0.3,1.25))
        ax.set_xlabel("Time (fs)")
        ax.set_ylabel("Density (g/cm3)")
        ax.set_title("Density vs. Time")
        m=['^','.','o','*','+','-']
        c=['r','y','b','g','m','k']
        n=int(sys.argv[1])
        loop(n)
        fig.show()
	raw_input()

if __name__=="__main__":
        MultiPlot()
