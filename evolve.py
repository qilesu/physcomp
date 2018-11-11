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

def f1(temp): #temp in Kelvin
        t = temp-273
        return -3.11e-4*t*t+3.48e-2*t + 0.0265

def f2(temp):
        t = temp-273
        return 2.142e-4*t*t-2.356e-2*t+1.348

def timeEvolve(pile, dt):   
        laplacianT = laplacianFunc(pile.meshTemp, pile.dx, pile.dy)
        laplacianO2 = laplacianFunc(pile.meshO2, pile.dx, pile.dy)
        ew = pile.voidFraction
        pCeff = (ew*1.17*1005+(1-ew)*ew*1150*3320)
        alpha = (ew*(0.026)+(1-ew)*(0.3))/pCeff#0.00145 # diffusivity
        dO2 = (0.176/10000)/np.power(273+25, 3/2)*np.power(273+45, 3/2)*ew #https://en.wikipedia.org/wiki/Mass_diffusivity
        Kp = 0.056 #kg/m3
        Ko = 10e-2 #mg/L
        Yo = 1.12 #mol/mol yield rate 
        Yt = 14e6/pCeff #K/kg proportional to O2 consumption
        # evolve non-edge cells
        eatenM = 0
        for i in range(1, len(pile.meshTemp)-1):
                for j in range(1, pile.topEdgeOfY):
                        dissolvedO2 = pile.meshO2[i][j]/(0.272)/ew*9.3 # mg/L
                        Xp = (1.0e-4)*(pile.bulkRho*pile.dryFraction)/(pile.bulkRho*pile.dryFraction+Kp)*dissolvedO2/(dissolvedO2+Ko)*pile.meshX[i][j]*f1(pile.meshTemp[i][j])*dt
                        Xm = (7.6e-5) * pile.meshX[i][j] * f2(pile.meshTemp[i][j]) * dt
                        pile.meshO2[i][j] += dO2*laplacianO2[i][j]*dt - Xp*(0.032)/Yo
                        pile.meshTemp[i][j] += alpha*laplacianT[i][j]*dt + Xp * (0.032)/Yo * Yt
                        dM = Xp/Yo*0.180/6
                        eatenM += dM
                        pile.meshX[i][j] += Xp-Xm
                        #if(j==5 and i==5): 
                                #print ("f2(T)", f2(pile.meshTemp[i][j]))
                                #print ("dissolved O2: ", dissolvedO2)
                                #print ("Kp factor: ", (pile.bulkRho*pile.dryFraction)/(pile.bulkRho*pile.dryFraction+Kp))
                                #print ("Xp:", Xp)
                                #print ("Xm:", Xm)
        #print("eatenM:", eatenM)

        for x in range(1, len(pile.meshTemp)-1):
                pile.meshO2[x][0] = pile.meshO2[x][1]

        pile.eatMass(eatenM)
        setBoundaryConditions(pile)

def setBoundaryConditions(pile):
        pile.setBoundaryCondition(pile.meshTemp, pile.iniT, pile.iniT, pile.iniT, pile.iniT)
        pile.setBoundaryCondition(pile.meshTemp, pile.iniO2, pile.iniO2, pile.iniO2, pile.iniO2)
