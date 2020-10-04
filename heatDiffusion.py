#!/usr/bin/env python
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import evolve
from Pile import Pile
import math

def interiorAverage(grid, topZ):
	sample = grid[1:-1, 1:-1, 1:topZ]
	return np.average(sample)
	
def main():

	for h in range(0,5):
		runSimu(0.5+0.2*h, 0.4)#H: 0.5-1.5
		#pile = Pile(1, 0.8+h*0.6/5, 0.05, 0.05, 273+20)
		#print("height,mass,rho")
		#for x in range(0, 40):	
			#print("%.3f,%.4f,%.3f" % (pile.height, pile.mass, pile.bulkRho))
			#pile.eatMass(pile.mass*0.04)

def plotTemp(pile, step_n, index, filename, dt, time_stamps, Temp, ATemp, start=0):
	fig = plt.figure(figsize=(16, 12))
	ax1 = fig.add_subplot(111)
	im1 = ax1.imshow(np.transpose(np.fliplr(pile.meshTemp[index])))
	fig.colorbar(im1)
	ax1.set_title("Temperature %f Kelvin." % (interiorAverage(pile.meshTemp, pile.topEdgeOfZ)))
	ax2 = fig.add_subplot(122)
	ax2.scatter(time_stamps, Temp)
	ax2.scatter(time_stamps, ATemp)
	#plt.axhline(y=273+initialT, color='r', linestyle='-')
	ax2.set_xlabel("time (min)")
	plt.savefig('%s-Temperature-%d-%d-%d-ave.png'%(filename, start, start+step_n, dt), transparent=True)
			
def plotOxygen(pile, step_n, index, filename, dt, time_stamps, Oxygen, ref, start=0):
	fig = plt.figure(figsize=(16, 12))
	ax1 = fig.add_subplot(121)
	#plt.subplot(1, 2, 1)
	im1 =ax1.imshow(np.transpose(np.fliplr(pile.meshO2[index])))
	fig.colorbar(im1)
	ax1.set_title("Oxygen %f kg/m3." % (interiorAverage(pile.meshO2, pile.topEdgeOfZ)))
	ax2 = fig.add_subplot(122)
	ax2.scatter(time_stamps, Oxygen)
	ax2.set_xlabel("time (min)")
	plt.axhline(y=ref, color='r', linestyle='-')
	plt.savefig('%s-Oxygen-%d-%d-%d-ave.png'%(filename, start, start+step_n, dt), transparent=True)

def plotBacteria(pile, step_n, index, filename, dt, time_stamps, Bacteria, start=0):
	fig = plt.figure(figsize=(16, 12))
	ax1 = fig.add_subplot(121)
	im1 = ax1.imshow(np.transpose(np.fliplr(pile.meshX[index])))
	fig.colorbar(im1)
	ax1.set_title("Bacteria %f mol/m3." % (interiorAverage(pile.meshX, pile.topEdgeOfZ)))
	ax2 = fig.add_subplot(122)
	ax2.scatter(time_stamps, Bacteria)
	ax2.set_xlabel("time (min)")
	plt.savefig('%s-Bacteria-%d-%d-%d-ave.png'%(filename, start, start+step_n, dt), transparent=True)


def runSimu(Lz, ds):
	Lx = 1 # length in x
	Ly = 1 # length in y
	dx = 0.07 # grid spacing m
	dt = 30 # seconds
	#initialT = 20
	filename = "3d-dis/-H-%.2fm-ds-%.2f%%-2mol-20-collapse" %(Lz, ds)
	print(filename)

	#meshTemp = np.full((round(Lx/dx), round(Ly/dx)), 20, dtype='float64') # initial temperature in C
	#setBoundaryCondition(meshTemp, 10, 15, 10, 10)
	pile = Pile(Lx, Ly, Lz, dx, dx, dx, 273+20, ds)
	initialMass = pile.mass
	#pile.loadFields()
	steps = 800
	start = 0
	log_step = 25
	time_stamps=[]
	ATemp = []
	Temp = []
	Oxygen = []
	Bacteria = []
	Mass = []
	print("percentage mass,TIME(h),HEIGHT,TEMP0,TEMP,OXY,BAC,MASS")
	timeNow = 0
	for i in range(start, start+steps):
		ambientT = 10*np.sin(2*np.pi/(24*3600)*timeNow) + 15 + 273
		evolve.timeEvolve(pile, dt, ambientT)
		timeNow += dt
		#output = ''
		#if(i>40 and i<75 and i%2==0):
			

		if(i%log_step==0):
			time_stamps.append(i*dt/3600)
			t=interiorAverage(pile.meshTemp, pile.topEdgeOfZ)
			o=interiorAverage(pile.meshO2, pile.topEdgeOfZ)
			b=interiorAverage(pile.meshX, pile.topEdgeOfZ)
			#t = pile.meshTemp[pile.meshTemp.shape[0]//2, pile.topEdgeOfY//2]
			#o = pile.meshO2[pile.meshO2.shape[0]//2, pile.topEdgeOfY//2]
			#b = pile.meshX[pile.meshX.shape[0]//2, pile.topEdgeOfY//2]
			ATemp.append(ambientT)
			Temp.append(t)
			
			Oxygen.append(o)
			Bacteria.append(b)
			Mass.append(pile.mass)

			plotTemp(pile, i, int(Lx/dx)//2, filename, dt, time_stamps, Temp, ATemp, start=0)
			plotOxygen(pile, i, int(Lx/dx)//2, filename, dt, time_stamps, Oxygen, 0.272, start=0)
			plotBacteria(pile, i, int(Lx/dx)//2, filename, dt, time_stamps, Bacteria, start=0)
			
		if(i%(log_step*4)==0):
			print("%.4f%%,%f,%f,%f, %f, %f, %f, %f" % (pile.mass/initialMass*100, i*dt/3600, pile.height,ambientT,t,o,b, pile.mass))
			plt.close('all')
			if(pile.mass*pile.height<0):
				break

	#print (pile.meshTemp)
	#pile.saveFields()

	#https://stackoverflow.com/questions/23876588/matplotlib-colorbar-in-each-subplot
	plotTemp(pile, steps, int(Lx/dx)//2, filename, dt, time_stamps, Temp, ATemp, start=0)
	plotOxygen(pile, steps, int(Lx/dx)//2, filename, dt, time_stamps, Oxygen, 0.272, start=0)
	plotBacteria(pile, steps, int(Lx/dx)//2, filename, dt, time_stamps, Bacteria, start=0)
	
	fig = plt.figure(figsize=(16, 12))
	plt.scatter(time_stamps, Mass)
	plt.title("Remaining Mass: %f%%"%(pile.mass/initialMass*100))
	plt.xlabel("time (min)")
	fig.savefig('%s-Mass-%d-%d-%d-ave.png'%(filename, start, start+steps, dt), transparent=True)

if __name__ == '__main__':
	sys.exit(main())
