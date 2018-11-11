import numpy as np

class Pile: 
	def __init__(self, Lx, initialH, dx, dy, initialT):
		initialO2 = 0.272 #kg m-3 simple paper
		initialX = 2e-2 #mol/m3
		self.meshTemp = np.full((round(Lx/dx), round(initialH/dy)), initialT, dtype='float64') # initial temperature in C
		self.meshO2 = np.full((round(Lx/dx), round(initialH/dy)), initialO2+0.01, dtype='float64') 
		self.meshX = np.full((round(Lx/dx), round(initialH/dy)), initialX, dtype='float64')
		self.setBoundaryCondition(self.meshTemp, 273+10, 273+15, 273+10, 273+10)
		self.setBoundaryCondition(self.meshO2, initialO2, initialO2, initialO2, initialO2)
		self.dx = dx
		self.dy = dy
		self.rhoMin = 100 #kg m-3
		self.resistance = 73/9.8 # E/g
		self.dryFraction = 0.3 #mass fraction
		self.area = 1
		self.height = initialH
		self.initMass()
		self.averageRho()
		self.computeVoidFraction() #volume fraction

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

	def computeVoidFraction(self):
		self.voidFraction = 1-(self.dryFraction/1150+1/1000-self.dryFraction/1000)*self.bulkRho

	def setBoundaryCondition(self, grid, top, bottom, left, right):
		for j in range(1, len(grid[0])-1):
			grid[0][j] = left # top edge
			grid[len(grid)-1][j] = right # bottom edge
		for i in range(1, len(grid)-1):
			grid[i][0] = bottom # left edge
			grid[i][len(grid[0])-1] = top # right edge
