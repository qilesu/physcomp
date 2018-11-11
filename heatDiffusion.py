#!/usr/bin/env python
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import evolve
from Pile import Pile
import math

def interiorAverage(grid):
	sample = grid[1:len(grid), 1:len(grid[0])-1]
	return np.average(sample)
	
def main():
	for h in range(1):
		runSimu(0.8+h*0.05)


def runSimu(Ly):
	Lx = 1 # length in x
	#Ly = 1 # length in y
	dx = 0.05 # grid spacing m
	dt = 10 # seconds
	initialT = 35

	#meshTemp = np.full((round(Lx/dx), round(Ly/dx)), 20, dtype='float64') # initial temperature in C
	#setBoundaryCondition(meshTemp, 10, 15, 10, 10)
	pile = Pile(Lx, Ly, dx, dx, 273+initialT)
	initialMass = pile.mass
	#pile.loadFields()
	steps = 8000
	start = 0
	log_step = 25
	time_stamps=[]
	Temp = []
	Oxygen = []
	Bacteria = []
	Mass = []
	for i in range(start, start+steps):
		evolve.timeEvolve(pile, dt)
		if(i%log_step==0):
			time_stamps.append(i*dt/60)
			Temp.append(interiorAverage(pile.meshTemp))
			Oxygen.append(interiorAverage(pile.meshO2))
			Bacteria.append(interiorAverage(pile.meshX))
			Mass.append(pile.mass)
		if(i%(log_step*4)==0):
			print("%f%%" % (pile.mass/initialMass*100))

	#print (pile.meshTemp)
	#pile.saveFields()

	filename = "H-%.2fm-2mol-%d" %(Ly, initialT)
	print(filename)

	#https://stackoverflow.com/questions/23876588/matplotlib-colorbar-in-each-subplot
	fig = plt.figure(figsize=(16, 12))
	ax1 = fig.add_subplot(121)
	im1 = ax1.imshow(np.transpose(np.fliplr(pile.meshTemp)))
	fig.colorbar(im1)
	ax1.set_title("Temperature %f Kelvin." % (interiorAverage(pile.meshTemp)))
	ax2 = fig.add_subplot(122)
	ax2.scatter(time_stamps, Temp)
	plt.axhline(y=273+initialT, color='r', linestyle='-')
	ax2.set_xlabel("time (min)")
	#fig = plt.figure(figsize=(244.0/72, 140.0/72))
	plt.savefig('%s-Temperature-%d-%d-%d.png'%(filename, start, start+steps, dt), transparent=True)
	
	fig = plt.figure(figsize=(16, 12))
	ax1 = fig.add_subplot(121)
	#plt.subplot(1, 2, 1)
	im1 =ax1.imshow(np.transpose(np.fliplr(pile.meshO2)))
	fig.colorbar(im1)
	ax1.set_title("Oxygen %f kg/m3." % (interiorAverage(pile.meshO2)))
	ax2 = fig.add_subplot(122)
	ax2.scatter(time_stamps, Oxygen)
	ax2.set_xlabel("time (min)")
	plt.axhline(y=0.272, color='r', linestyle='-')
	plt.savefig('%s-Oxygen-%d-%d-%d.png'%(filename, start, start+steps, dt), transparent=True)

	fig = plt.figure(figsize=(16, 12))
	ax1 = fig.add_subplot(121)
	im1 = ax1.imshow(np.transpose(np.fliplr(pile.meshX)))
	fig.colorbar(im1)
	ax1.set_title("Bacteria %f mol/m3." % (interiorAverage(pile.meshX)))
	ax2 = fig.add_subplot(122)
	ax2.scatter(time_stamps, Bacteria)
	ax2.set_xlabel("time (min)")
	plt.savefig('%s-Bacteria-%d-%d-%d.png'%(filename, start, start+steps, dt), transparent=True)

	fig = plt.figure(figsize=(16, 12))
	plt.scatter(time_stamps, Mass)
	plt.title("Remaining Mass: %f%%"%(pile.mass/initialMass*100))
	plt.xlabel("time (min)")
	fig.savefig('%s-Mass-%d-%d-%d.png'%(filename, start, start+steps, dt), transparent=True)

if __name__ == '__main__':
	sys.exit(main())
