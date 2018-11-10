#!/usr/bin/env python
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

Lx = 1 # length in x
Ly = 1 # length in y
dx = 0.2 # grid spacing m
dt = 1 # seconds
alpha = 0.00145 # diffusivity

meshTemp= np.full((round(Lx/dx), round(Ly/dx)), 20, dtype='float64') # initial temperature in C

def timeEvolve(grid, diffusivity, dx, dy, dt, steps):
    laplacian = np.zeros_like(grid) # stores lapacian value
    for i in range(steps):
        # calculate laplacian of temperature at non-edge cell
        # with finite difference
        for i in range(1, len(grid)-1):
            for j in range(1, len(grid[0])-1):
                laplacian[i][j] = (grid[i-1][j]+grid[i+1][j]-2*grid[i][j])/(dx*dx)
                laplacian[i][j] += (grid[i][j-1]+grid[i][j+1]-2*grid[i][j])/(dy*dy)

        # evolve non-edge cells
        for i in range(1, len(grid)-1):
            for j in range(1, len(grid[0])-1):
                grid[i][j] += diffusivity*dt*laplacian[i][j]

def setBoundaryCondition(grid, top, bottom, left, right):
    for j in range(1, len(grid[0])-1):
        grid[0][j] = top # top edge
        grid[len(meshTemp)-1][j] = bottom # bottom edge
    for i in range(1, len(meshTemp)-1):
        grid[i][0] = left # left edge
        grid[i][len(meshTemp[0])-1] = right # right edge
    
def main():
    fig = plt.figure(figsize=(244.0/72, 140.0/72))
    plt.savefig('heatDiffusion.pdf', transparent=True)

if __name__ == '__main__':
    sys.exit(main())
