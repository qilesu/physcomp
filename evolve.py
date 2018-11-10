import numpy as np

def heatDiffusionFunc(laplacian):
	alpha = 0.00145 # diffusivity
	return alpha*laplacian; 

def timeEvolve(grid, diffusionFunc, dx, dy, dt):

	laplacian = np.zeros_like(grid) # stores lapacian value
	#for i in range(steps):
	# calculate laplacian of temperature at non-edge cell
	# with finite difference
	for i in range(1, len(grid)-1):
		for j in range(1, len(grid[0])-1):
			laplacian[i][j] = (grid[i-1][j]+grid[i+1][j]-2*grid[i][j])/(dx*dx)
			laplacian[i][j] += (grid[i][j-1]+grid[i][j+1]-2*grid[i][j])/(dy*dy)

	# evolve non-edge cells
	for i in range(1, len(grid)-1):
		for j in range(1, len(grid[0])-1):
			result = diffusionFunc(laplacian[i][j])
			grid[i][j] += result*dt #diffusivity*dt*laplacian[i][j] #specific function of (laplacian, alpha, grid)