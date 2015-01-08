#!/usr/bin/python
import numpy as np
import math
import sys
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def ReadFile(f):
        lines=np.array([l for l in f])
        lines=[l.split() for l in lines]
        for i in range(len(lines)):
                if len(lines[i])==0:continue
		if lines[i][0]=='thermo':
			thermo = float(lines[i][1])
		if lines[i][0]=='Volume':
			index=i+1
			break
	vol=[]
        for i in range(index,len(lines)):
                if lines[i][0]=='Loop':
                        break
                else:
                    time = i*thermo
                    vol.append([time, float(lines[i][0])])
        return np.array(vol)

def Plot(n):
        global ax, m, c
        fig=plt.figure()
        m=['^','.','o','*','+','-']
        c=['r','y','b','g','m','k']
        for i in range(n):
                f=open(sys.argv[i+2],'r')
                line=ReadFile(f)
		ymin = min(line[:,1])
		ymin = ymin-0.5*ymin
		ymax = ymin*6
		print min(line[:,0])
        	ax=fig.add_subplot(111,ylim=(ymin,ymax))
                ax.plot(line[:,0],line[:,1],color=c[i%len(c)])
		ax.set_xlabel("Time (fs)")
		ax.set_ylabel("Volume (Angstrom3)")
		ax.set_title("Volume vs. Time")
        fig.show()

if __name__=="__main__":
	Plot(int(sys.argv[1]))
	raw_input()
