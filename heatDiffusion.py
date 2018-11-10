#!/usr/bin/env python
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import evolve
from Pile import Pile


	
def main():
	Lx = 1 # length in x
	Ly = 1 # length in y
	dx = 0.01 # grid spacing m
	dt = 1 # seconds

	#meshTemp = np.full((round(Lx/dx), round(Ly/dx)), 20, dtype='float64') # initial temperature in C
	#setBoundaryCondition(meshTemp, 10, 15, 10, 10)
	pile = Pile(Lx, Ly, dx, dx, 20)
	steps = 1000
	for i in range(steps):
		evolve.timeEvolve(pile, dt)

	print (pile.meshTemp)
	plt.subplot(1, 2, 1)
	plt.imshow(np.transpose(np.fliplr(pile.meshTemp)))
	plt.subplot(1, 2, 2)
	plt.imshow(np.transpose(np.fliplr(pile.meshO2)))
	#fig = plt.figure(figsize=(244.0/72, 140.0/72))
	plt.savefig('Diffusion-square-1000.png', transparent=True)

if __name__ == '__main__':
	sys.exit(main())
