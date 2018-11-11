#!/usr/bin/env python
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import evolve
from Pile import Pile
import math

def interiorAverage(grid, topY):
	sample = grid[1:len(grid), 1:topY]
	return np.average(sample)
	
def main():

	for h in range(0,1):
		runSimu(0.82, 0.1)#H: 0.5-1.5
		#pile = Pile(1, 0.8+h*0.6/5, 0.05, 0.05, 273+20)
		#print("height,mass,rho")
		#for x in range(0, 40):	
			#print("%.3f,%.4f,%.3f" % (pile.height, pile.mass, pile.bulkRho))
			#pile.eatMass(pile.mass*0.04)

def runSimu(Ly, ds):
	Lx = 2 # length in x
	#Ly = 1 # length in y
	dx = 0.05 # grid spacing m
	dt = 30 # seconds
	#initialT = 20
	filename = "grid/H-%.2fm-ds-%.2f%%-2mol-15-collapse" %(Ly, ds)
	print(filename)

	#meshTemp = np.full((round(Lx/dx), round(Ly/dx)), 20, dtype='float64') # initial temperature in C
	#setBoundaryCondition(meshTemp, 10, 15, 10, 10)
	pile = Pile(Lx, Ly, dx, dx, 273+15, ds)
	initialMass = pile.mass
	#pile.loadFields()
	steps = 2000
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
		output = ''
		if(i%log_step==0):
			time_stamps.append(i*dt/3600)
			t=interiorAverage(pile.meshTemp, pile.topEdgeOfY)
			o=interiorAverage(pile.meshO2, pile.topEdgeOfY)
			b=interiorAverage(pile.meshX, pile.topEdgeOfY)
			ATemp.append(ambientT)
			Temp.append(t)
			
			Oxygen.append(o)
			Bacteria.append(b)
			Mass.append(pile.mass)
		if(i%(log_step*4)==0):
			print("%.4f%%,%f,%f,%f, %f, %f, %f, %f" % (pile.mass/initialMass*100, i*dt/3600, pile.height,ambientT,t,o,b, pile.mass))
			if(pile.mass*pile.height<0):
				break

	#print (pile.meshTemp)
	#pile.saveFields()

	#https://stackoverflow.com/questions/23876588/matplotlib-colorbar-in-each-subplot
	fig = plt.figure(figsize=(16, 12))
	ax1 = fig.add_subplot(121)
	im1 = ax1.imshow(np.transpose(np.fliplr(pile.meshTemp)))
	fig.colorbar(im1)
	ax1.set_title("Temperature %f Kelvin." % (interiorAverage(pile.meshTemp, pile.topEdgeOfY)))
	ax2 = fig.add_subplot(122)
	ax2.scatter(time_stamps, Temp)
	ax2.scatter(time_stamps, ATemp)
	#plt.axhline(y=273+initialT, color='r', linestyle='-')
	ax2.set_xlabel("time (min)")
	#fig = plt.figure(figsize=(244.0/72, 140.0/72))
	plt.savefig('%s-Temperature-%d-%d-%d.png'%(filename, start, start+steps, dt), transparent=True)
	
	fig = plt.figure(figsize=(16, 12))
	ax1 = fig.add_subplot(121)
	#plt.subplot(1, 2, 1)
	im1 =ax1.imshow(np.transpose(np.fliplr(pile.meshO2)))
	fig.colorbar(im1)
	ax1.set_title("Oxygen %f kg/m3." % (interiorAverage(pile.meshO2, pile.topEdgeOfY)))
	ax2 = fig.add_subplot(122)
	ax2.scatter(time_stamps, Oxygen)
	ax2.set_xlabel("time (min)")
	plt.axhline(y=0.272, color='r', linestyle='-')
	plt.savefig('%s-Oxygen-%d-%d-%d.png'%(filename, start, start+steps, dt), transparent=True)

	fig = plt.figure(figsize=(16, 12))
	ax1 = fig.add_subplot(121)
	im1 = ax1.imshow(np.transpose(np.fliplr(pile.meshX)))
	fig.colorbar(im1)
	ax1.set_title("Bacteria %f mol/m3." % (interiorAverage(pile.meshX, pile.topEdgeOfY)))
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
