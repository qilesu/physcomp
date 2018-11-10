import numpy as np
#
#def heatDiffusionFunc(laplacian):
	#alpha = 0.00145 # diffusivity
	#return alpha*laplacian; 
def laplacianFunc(grid, dx, dy):
	laplacian = np.zeros_like(grid) # stores lapacian value
	#for i in range(steps):
	# calculate laplacian of temperature at non-edge cell
	# with finite difference
	for i in range(1, len(grid)-1):
		for j in range(1, len(grid[0])-1):
			laplacian[i][j] = (grid[i-1][j]+grid[i+1][j]-2*grid[i][j])/(dx*dx)
			laplacian[i][j] += (grid[i][j-1]+grid[i][j+1]-2*grid[i][j])/(dy*dy)

	return laplacian


def timeEvolve(pile, dt):
	laplacianT = laplacianFunc(pile.meshTemp, pile.dx, pile.dy)
	alpha = 0.00145 # diffusivity
	# evolve non-edge cells
	for i in range(1, len(pile.meshTemp)-1):
		for j in range(1, len(pile.meshTemp[0])-1):
			pile.meshTemp[i][j] += alpha*laplacianT[i][j]*dt #diffusivity*dt*laplacian[i][j] #specific function of (laplacian, alpha, grid)