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
	laplacianO2 = laplacianFunc(pile.meshO2, pile.dx, pile.dy)
	ew = pile.voidFraction
	alpha = (ew*(0.026)+(1-ew)*(0.3))/(ew*1.17*1005+(1-ew)*ew*1150*3320)#0.00145 # diffusivity
	dO2 = (0.176/10000)/np.power(273+25, 3/2)*np.power(273+45, 3/2)*ew #https://en.wikipedia.org/wiki/Mass_diffusivity
	# evolve non-edge cells
	for i in range(1, len(pile.meshTemp)-1):
		for j in range(1, len(pile.meshTemp[0])-1):
			pile.meshTemp[i][j] += alpha*laplacianT[i][j]*dt #diffusivity*dt*laplacian[i][j] #specific function of (laplacian, alpha, grid)
			pile.meshO2[i][j] += dO2*laplacianO2[i][j]*dt

	for x in range(0, len(pile.meshTemp)-1):
		pile.meshO2[x][0] = pile.meshO2[x][1]