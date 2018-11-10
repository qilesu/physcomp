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
        dx = 0.05 # grid spacing m
        dt = 5 # seconds

        #meshTemp = np.full((round(Lx/dx), round(Ly/dx)), 20, dtype='float64') # initial temperature in C
        #setBoundaryCondition(meshTemp, 10, 15, 10, 10)
        pile = Pile(Lx, Ly, dx, dx, 273+35)
        pile.loadFields()
        steps = 100
        for i in range(steps):
                evolve.timeEvolve(pile, dt)

        print (pile.meshTemp)
        plt.subplot(1, 3, 1)
        plt.imshow(np.transpose(np.fliplr(pile.meshTemp)))
        plt.subplot(1, 3, 2)
        plt.imshow(np.transpose(np.fliplr(pile.meshO2)))
        plt.subplot(1, 3, 3)
        plt.imshow(np.transpose(np.fliplr(pile.meshX)))
        plt.colorbar()
        #fig = plt.figure(figsize=(244.0/72, 140.0/72))
        plt.savefig('Diffusion-run.png', transparent=True)
        pile.saveFields()

if __name__ == '__main__':
        sys.exit(main())
