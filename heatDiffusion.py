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
	Lx = 1 # length in x
	Ly = 1 # length in y
	dx = 0.05 # grid spacing m
	dt = 5 # seconds
	initialT = 20

	#meshTemp = np.full((round(Lx/dx), round(Ly/dx)), 20, dtype='float64') # initial temperature in C
	#setBoundaryCondition(meshTemp, 10, 15, 10, 10)
	pile = Pile(Lx, Ly, dx, dx, 273+initialT)
	steps = 1000
	Temp = []
	Oxygen = []
	Bacteria = []
	Mass = []
	for i in range(steps):
		evolve.timeEvolve(pile, dt)
		if(i%20==0):
			Temp.append(interiorAverage(pile.meshTemp))
			Oxygen.append(interiorAverage(pile.meshO2))
			Bacteria.append(interiorAverage(pile.meshX))
			Mass.append(pile.mass)

	print (pile.meshTemp)

	#https://stackoverflow.com/questions/23876588/matplotlib-colorbar-in-each-subplot
	fig = plt.figure(figsize=(16, 12))
	ax1 = fig.add_subplot(121)
	im1 = ax1.imshow(np.transpose(np.fliplr(pile.meshTemp)))
	fig.colorbar(im1)
	ax1.set_title("Temperature %f Kelvin." % (interiorAverage(pile.meshTemp)))
	ax2 = fig.add_subplot(122)
	ax2.scatter(range(0,math.floor(steps/20)), Temp)
	ax2.set_xlabel("time (s)")
	#fig = plt.figure(figsize=(244.0/72, 140.0/72))
	plt.savefig('Diffusion-Run-Temperature-%d-%d.png'%(steps, dt), transparent=True)
	
	fig = plt.figure(figsize=(16, 12))
	ax1 = fig.add_subplot(121)
	#plt.subplot(1, 2, 1)
	im1 =ax1.imshow(np.transpose(np.fliplr(pile.meshO2)))
	fig.colorbar(im1)
	ax1.set_title("Oxygen %f kg/m3." % (interiorAverage(pile.meshO2)))
	ax2 = fig.add_subplot(122)
	ax2.scatter(range(0,math.floor(steps/20)), Oxygen)
	ax2.set_xlabel("time (s)")
	plt.savefig('Diffusion-Run-Oxygen-%d-%d.png'%(steps, dt), transparent=True)

	fig = plt.figure(figsize=(16, 12))
	ax1 = fig.add_subplot(121)
	im1 = ax1.imshow(np.transpose(np.fliplr(pile.meshX)))
	fig.colorbar(im1)
	ax1.set_title("Bacteria %f mol/m3." % (interiorAverage(pile.meshX)))
	ax2 = fig.add_subplot(122)
	ax2.scatter(range(0,math.floor(steps/20)), Bacteria)
	ax2.set_xlabel("time (s)")
	plt.savefig('Diffusion-Run-Bacteria-%d-%d.png'%(steps, dt), transparent=True)

	fig = plt.figure(figsize=(16, 12))
	plt.scatter(range(0, math.floor(steps/20)), Mass)
	plt.title("Mass versus Time")
	plt.xlabel("time (s)")
	fig.savefig('Diffusion-Run-Mass-%d-%d.png'%(steps, dt), transparent=True)

	#fig = plt.figure(figsize=(244.0/72, 140.0/72))
	#plt.savefig('Diffusion-Run-Oxygen-Bacteria-%d-%d.png'%(steps, dt), transparent=True)

if __name__ == '__main__':
	sys.exit(main())
