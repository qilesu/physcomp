import numpy as np

class Pile: 
        def __init__(self, Lx, initialH, dx, dy, initialT, ds):
                self.iniO2 = 0.272 #kg m-3 simple paper
                self.iniT = initialT # K
                self.iniX = 2#2 #mol/m3
                self.meshTemp = np.full((round(Lx/dx)+2, round(initialH/dy)+2), self.iniT, dtype='float64') # initial temperature in C
                self.meshO2 = np.full((round(Lx/dx)+2, round(initialH/dy)+2), self.iniO2, dtype='float64') 
                self.meshX = np.full((round(Lx/dx)+2, round(initialH/dy)+2), self.iniX, dtype='float64')
                self.dx = dx
                self.dy = dy
                self.rhoMin = 100 #kg m-3
                self.resistance = 73/9.8 # E/g
                self.dryFraction = ds#0.3 #mass fraction
                self.area = 1
                self.height = initialH
                self.topEdgeOfY = int(round(initialH/dy))
                self.initMass()
                self.averageRho()
                self.computeVoidFraction() #volume fraction
                self.setBoundaryCondition(self.meshTemp, self.iniT, self.iniT, self.iniT, self.iniT)
                self.setBoundaryCondition(self.meshO2, self.iniO2, self.iniO2, self.iniO2, self.iniO2)
                
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
                self.topEdgeOfY = int(round(self.height/self.dy))

        def computeVoidFraction(self):
                self.voidFraction = 1-(self.dryFraction/1150+1/1000-self.dryFraction/1000)*self.bulkRho

        def setBoundaryCondition(self, grid, top, bottom, left, right):
                # j in range(1, len(grid[0])-1):
                #        grid[0][j] = left # left edge
                #        grid[len(grid)-1][j] = right # right edge
                grid[0, 1:-1] = left
                grid[-1, 1:-1] = right
                #for i in range(1, len(grid)-1):
                #        grid[i][0] = bottom # bottom edge
                #        grid[i][self.topEdgeOfY] = top # top edge
                grid[1:-1, 0] = bottom
                grid[1:-1, self.topEdgeOfY:] = top
                #grid[:, self.topEdgeOfY+1:] = top

        def saveFields(self):
                np.savetxt("meshTemp.dat", self.meshTemp)
                np.savetxt("meshO2.dat", self.meshO2)
                np.savetxt("meshX.dat", self.meshX)
                # store loose variables
                empty = np.array([])
                headerString = str(self.iniO2) + "\n"
                headerString = str(self.iniT) + "\n"
                headerString = str(self.iniX) + "\n"
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
                        self.iniO2 = A[0]
                        self.iniT = A[1]
                        self.iniX = A[2]
                        self.dx = A[3]
                        self.dy = A[4]
                        self.rhoMin = A[5]
                        self.resistance = A[6]
                        self.dryFraction = A[7]
                        self.area = A[8]
                        self.height = A[9]
                        self.topEdgeOfY = A[10]
                        self.mass = A[11]
                        self.bulkRho = A[12]
                        self.voidFraction = A[13]
                except:
                        print("ERROR OCCURED WHEN LOADING SAVED FILE.")
