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
	dt = 0.01 # seconds
	#alpha = 0.00145 # diffusivity

	#meshTemp = np.full((round(Lx/dx), round(Ly/dx)), 20, dtype='float64') # initial temperature in C

	#setBoundaryCondition(meshTemp, 10, 15, 10, 10)
	pile = Pile(Lx, Ly, dx, dx, 20, 1)
	steps = 100
	for i in range(steps):
		evolve.timeEvolve(pile.meshTemp, evolve.heatDiffusionFunc, pile.dx, pile.dy, dt)

	print (pile.meshTemp)
	plt.imshow(pile.meshTemp)
	#fig = plt.figure(figsize=(244.0/72, 140.0/72))
	plt.savefig('heatDiffusion.png', transparent=True)

if __name__ == '__main__':
	sys.exit(main())
