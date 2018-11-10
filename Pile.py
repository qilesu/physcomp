import numpy as np

class Pile: 
	def __init__(self, Lx, Ly, dx, dy, initialT, initialH):
		self.meshTemp = np.full((round(Lx/dx), round(Ly/dy)), initialT, dtype='float64') # initial temperature in C
		self.setBoundaryCondition(self.meshTemp, 10, 15, 10, 10)
		self.dx = dx
		self.dy = dy
		self.rhoMin = 100 #kg m-3
		self.resistance = 73/9.8 # E/g
		self.dryFraction = 0.3 #mass fraction
		self.area = 1
		self.height = initialH
		self.initMass()
		self.averageRho()
		self.voidFraction() #volume fraction

	def initMass(self):
		self.mass = self.area*self.rhoMin*(np.exp(self.height/(self.resistance*self.dryFraction))-1)*self.resistance*self.dryFraction

	def computeHeight(self):
		self.height = self.dryFraction*self.resistance*np.log(1+self.mass/(self.area*self.rhoMin*self.resistance*self.dryFraction))

	def averageRho(self):
		self.bulkRho = self.mass/(self.area*self.height)

	def eatMass(self, decrement):
		self.mass -= decrement
		self.computeHeight()
		self.averageRho()

	def voidFraction(self):
		return 1-(self.dryFraction/1150+1/1000-self.dryFraction/1000)*self.bulkRho

	def setBoundaryCondition(self, grid, top, bottom, left, right):
		for j in range(1, len(grid[0])-1):
			grid[0][j] = top # top edge
			grid[len(grid)-1][j] = bottom # bottom edge
		for i in range(1, len(grid)-1):
			grid[i][0] = left # left edge
			grid[i][len(grid[0])-1] = right # right edge
