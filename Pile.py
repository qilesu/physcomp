import numpy as np

class Pile: 
	def __init__(self, Lx, Ly, dx, dy, initialT, initialM):
		self.meshTemp = np.full((round(Lx/dx), round(Ly/dy)), initialT, dtype='float64') # initial temperature in C
		self.setBoundaryCondition(self.meshTemp, 10, 15, 10, 10)
		self.dx = dx
		self.dy = dy


	def setBoundaryCondition(self, grid, top, bottom, left, right):
		for j in range(1, len(grid[0])-1):
			grid[0][j] = top # top edge
			grid[len(grid)-1][j] = bottom # bottom edge
		for i in range(1, len(grid)-1):
			grid[i][0] = left # left edge
			grid[i][len(grid[0])-1] = right # right edge
