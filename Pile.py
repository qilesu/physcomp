import numpy as np

class Pile: 
	def __init__(self, Lx, initialH, dx, dy, initialT):
		initialO2 = 0.272 #kg m-3 simple paper
		initialX = 2e-2 #mol/m3
		self.meshTemp = np.full((round(Lx/dx)+2, round(initialH/dy))+2, initialT, dtype='float64') # initial temperature in C
		self.meshO2 = np.full((round(Lx/dx)+2, round(initialH/dy))+2, initialO2+0.01, dtype='float64') 
		self.meshX = np.full((round(Lx/dx)+2, round(initialH/dy))+2, initialX, dtype='float64')
		self.setBoundaryCondition(self.meshTemp, 273+10, 273+15, 273+10, 273+10)
		self.setBoundaryCondition(self.meshO2, initialO2, initialO2, initialO2, initialO2)
		self.dx = dx
		self.dy = dy
		self.rhoMin = 100 #kg m-3
		self.resistance = 73/9.8 # E/g
		self.dryFraction = 0.3 #mass fraction
		self.area = 1
		self.height = initialH
                self.topEdgeOfY = round(initialH/dy)
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
                # calculate number of grid in y direction that matters
                numGridInY = round(pile.height/pile.dy)
		for j in range(1, len(grid[0])-1):
			grid[0][j] = left # left edge
			grid[len(grid)-1][j] = right # right edge
		for i in range(1, len(grid)-1):
			grid[i][0] = bottom # bottom edge
			grid[i][numGridInY] = top # top edge

	def saveFields(self):
		np.savetxt("meshTemp.dat", self.meshTemp)
		np.savetxt("meshO2.dat", self.meshO2)
		np.savetxt("meshX.dat", self.meshX)
		# store loose variables
		empty = np.array([])
		headerString = str(self.dx) + "\n"
		headerString += str(self.dy) + "\n"
		headerString += str(self.rhoMin) + "\n"
		headerString += str(self.resistance) + "\n"
		headerString += str(self.dryFraction) + "\n"
		headerString += str(self.area) + "\n"
		headerString += str(self.height) + "\n"
                headerString += str(self.topEdgeOfY) + "\n"
		headerString += str(self.mass) + "\n"
		headerString += str(self.bulkRho) + "\n"
		headerString += str(self.voidFraction)
		np.savetxt("otherVariables.dat", empty, header=headerString, comments='')
	def loadFields(self):
		try:
				self.meshTemp = np.loadtxt("meshTemp.dat")
				self.meshO2 = np.loadtxt("meshO2.dat")
				self.meshX = np.loadtxt("meshX.dat")
				A = np.loadtxt("otherVariables.dat")
				self.dx = A[0]
				self.dy = A[1]
				self.rhoMin = A[2]
				self.resistance = A[3]
				self.dryFraction = A[4]
				self.area = A[5]
				self.height = A[6]
                                self.topEdgeOfY = A[7]
				self.mass = A[8]
				self.bulkRho = A[9]
				self.voidFraction = A[10]
		except:
				print("ERROR OCCURED WHEN LOADING SAVE FILE.")
