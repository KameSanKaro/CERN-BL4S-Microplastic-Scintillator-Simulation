import math
import random
from matplotlib import pyplot as plt

class screp:
    def __init__(self, mainpart, magnitude):
        self.mp = mainpart
        self.mag = magnitude
    def tofloat(self):
        return self.mp*10**(self.mag)
    def srtag(self):
        return "(" + str(self.mp) + "," + str(self.mag) + ")"
    
def toscrep(conv):
    convstr = str(conv)
    newscr =  screp(conv/(10**(len(convstr)-1)),(len(convstr)-1))
    return newscr
    
def mulsr(sr1,sr2):
    if sr1.mp*sr2.mp < 10:
        return screp(sr1.mp*sr2.mp,sr1.mag+sr2.mag)
    else:
        return screp(sr1.mp*sr2.mp/10,sr1.mag+sr2.mag+1)

c = 299792458
h = screp(6.62607015,-34)
n = 1.49

class Scintillator:
    def __init__(self, ppmconcentration, x, y, z, resolution):
        self.ppmconcentration = ppmconcentration
        self.x = x
        self.y = y
        self.z = z
        self.resolution = resolution
        """Resolution should be formatted in the form of qty. of grid elements per cm. For instance, a 1mm resolution would formatted as 10:1, or 10.0."""
        def pollute(ppmconc,xval,yval,zval,res):

            scdict = {}

            for i in range(xval*res):
                for j in range(yval*res):
                    for k in range(zval*res):
                        if random.random() < ppmconc/10**6:
                            scdict[(i/res,j/res,k/res)] = 1
                        else:
                            scdict[(i/res,j/res,k/res)] = 0
            
            return scdict
        
        self.grid = pollute(self.ppmconcentration,self.x,self.y,self.z,self.resolution)

    def getConcentration(self):
        return self.ppmconcentration

    def getDimension(self):
        return [self.x, self.y, self.z]
    
    def getResolution(self):
        return self.resolution

    def scintillatorTag(self):
        return "A IM Liquid - PMMA scintillator with a PMMA-MP concentration of " + str(self.ppmconcentration) +" ppm and dimensions " + str((self.x,self.y,self.z)) + " cm. (Resolution: " + str(self.resolution) + " lattice elements per cm)"
    
    def printGrid(self):
        grid = self.returnGrid()

        for k in range(self.z*self.resolution):
            print("Z = " + str(k))
            for j in range(self.y*self.resolution):
                linestr = ""
                for i in range(self.x*self.resolution):
                    linestr = linestr + str(grid[(i/self.resolution,self.y-(j+1)/self.resolution,k/self.resolution)])
                print(linestr)
    
    def returnGrid(self):
        return self.grid

    def returnMPcount(self):

        count = 0
        grid = self.returnGrid()
        
        for i in range(self.x*self.resolution):
            for j in range(self.y*self.resolution):
                for k in range(self.z*self.resolution):
                    if grid[(i/self.resolution,j/self.resolution,k/self.resolution)] == 1:
                        count += 1
        
        return count
    
    def getGridZoneCoords(self):

        gridzonecoordlist = []
        gridx = []
        gridy = []
        gridz = []

        for i in range(self.x*self.resolution):
            gridx.append(i/self.resolution)
        for j in range(self.y*self.resolution):
            gridy.append(j/self.resolution)
        for k in range(self.z*self.resolution):
            gridz.append(k/self.resolution)

        gridzonecoordlist = [gridx,gridy,gridz]

        return gridzonecoordlist
    
def MPcountchecker(numTrials,conc,x,y,z,res):

    totalcount = 0

    for trial in range(numTrials):
        gridcode = "sc" + str(trial)

        gridcode = Scintillator(conc,x,y,z,res)

        totalcount += gridcode.returnMPcount()

    return totalcount/numTrials

class Proton:
    def __init__(self,velocityx,velocityy,velocityz,xcoord,ycoord,zcoord):
        self.velocityx = velocityx #all in m/s, make sure to convert to cm/s before updating position
        self.velocityy = velocityy
        self.velocityz = velocityz
        self.xco = xcoord
        self.yco = ycoord
        self.zco = zcoord
        self.isPhoton = False
        self.isRemoved = False
        self.mass = screp(1.67262192,-27)

    def protonTag(self):
        if self.isPhoton == False:
            return "A proton travelling at " + str(math.sqrt((self.velocityx*100)**2+(self.velocityy*100)**2+(self.velocityz*100)**2)) + " cm/s. Position: " + str((self.xco,self.yco,self.zco))
        else:
            return "A photon travelling at " + str(math.sqrt((self.velocityx*100)**2+(self.velocityy*100)**2+(self.velocityz*100)**2)) + " cm/s. Position: " + str((self.xco,self.yco,self.zco))
    
    def getVelocity(self):
        return self.velocityx,self.velocityy,self.velocityz, "meters per second."

    def getPosition(self):
        return (self.xco,self.yco,self.zco)

    def updatePos(self,dt):
        "dt in seconds"
        self.xco = self.xco + (self.velocityx*100)*dt #self.velocities in m/s, position in cm on grid
        self.yco = self.yco + (self.velocityy*100)*dt
        self.zco = self.zco + (self.velocityz*100)*dt

    def convertToPhoton(self):
        self.isPhoton = True
        self.mass = 0

    def scintillate(self,dissipation,wavelength,quantum_efficiency,listtoadd):
        photonls = listtoadd

        m = self.mass.tofloat()
        old_sp = math.sqrt(self.velocityx**2+self.velocityy**2+self.velocityz**2) #speed in m/s
        #print("old speed: " + str(old_sp))
        lorentz = 1/math.sqrt(1-(old_sp/c)**2)
        #print("lorentz factor: " + str(lorentz))
        resting_energy = m*c**2
        #print("resting energy: " + str(resting_energy))
        old_energy = math.sqrt((lorentz*m*old_sp*c)**2+(m*c**2)**2)
        #print("old energy: " + str(old_energy))
        new_energy = old_energy - dissipation
        #print("new energy: " + str(new_energy))
        if new_energy < resting_energy:
            new_energy = resting_energy
        new_sp = c * math.sqrt(1-(resting_energy/new_energy)**2)
        #print("new speed: " + str(new_sp))
        
        if random.random() < quantum_efficiency:
            change_coef = new_sp/old_sp
            self.velocityx = self.velocityx*change_coef
            self.velocityy = self.velocityy*change_coef
            self.velocityz = self.velocityz*change_coef
            
            if new_energy > resting_energy:

                emissionval = mulsr(h,wavelength)
                dp = emissionval.tofloat()

                old_sp = math.sqrt(self.velocityx**2+self.velocityy**2+self.velocityz**2)
                lorentz = 1/math.sqrt(1-(old_sp/c)**2)
                new_momentum = lorentz*m*old_sp - dp
                if new_momentum < 0:
                    new_momentum = 0
                new_energy = math.sqrt((new_momentum*c)**2+(m*c**2)**2)
                new_sp = c * math.sqrt(1-(resting_energy/new_energy)**2)

                change_coef = new_sp/old_sp
                self.velocityx = self.velocityx*change_coef
                self.velocityy = self.velocityy*change_coef
                self.velocityz = self.velocityz*change_coef
                
                if new_energy > resting_energy:

                    pr_sp = math.sqrt(self.velocityx**2+self.velocityy**2+self.velocityz**2)
                    photon_speedcoef = c/(pr_sp*n)
                    ph_tag = self.protonTag() + str(self)
                    phvx = self.velocityx*photon_speedcoef
                    phvy = self.velocityy*photon_speedcoef
                    phvz = self.velocityz*photon_speedcoef
                    phx = self.xco
                    phy = self.yco
                    phz = self.zco
                    ph_tag = Proton(phvx,phvy,phvz,phx,phy,phz)
                    ph_tag.convertToPhoton()

                    photonls.append(ph_tag)
        else:
            change_coef = new_sp/old_sp
            self.velocityx = self.velocityx*change_coef
            self.velocityy = self.velocityy*change_coef
            self.velocityz = self.velocityz*change_coef

        return photonls

    def isitaPhoton(self):
        return self.isPhoton

    def isitRemoved(self):
        return self.isRemoved
    
    def removeproton(self):
        self.isRemoved = True

def protonBeam(numOfProtons, average_velocityx, average_velocityy, average_velocityz, velocity_variationx, velocity_variationy, velocity_variationz, position_margins):
    protonlist = []

    for pr_num in range(numOfProtons):
        pr_str = "pr" + str(pr_num)
        pr_str = Proton(random.normalvariate(average_velocityx,velocity_variationx),random.normalvariate(average_velocityy,velocity_variationy), random.normalvariate(average_velocityz,velocity_variationz), position_margins[0]+(position_margins[1]-position_margins[0])*random.random(), position_margins[2]+(position_margins[3]-position_margins[2])*random.random(), position_margins[4]+(position_margins[5]-position_margins[4])*random.random())
        protonlist.append(pr_str)

    return protonlist

def simulateTrial(numOfProtons, average_velocityx, average_velocityy, average_velocityz, velocity_variationx, velocity_variationy, velocity_variationz, beamMargins, pmma_conc, sc_x, sc_y, sc_z, sc_resolution, dtsr, dissipationsr, wavelength, quantum_efficiency, durationiniter):
    scTr = Scintillator(pmma_conc,sc_x,sc_y,sc_z,sc_resolution)
    protonTr = protonBeam(numOfProtons,average_velocityx,average_velocityy,average_velocityz,velocity_variationx,velocity_variationy,velocity_variationz,beamMargins)
    photonTr = []

    photondetectioncount = 0
    passcount = 0
    detectedphotons = []
    dt = dtsr.tofloat()
    dissipation = dissipationsr.tofloat()

    sc_dim = scTr.getDimension()
    sc_coordlist = scTr.getGridZoneCoords()
    sc_grid = scTr.returnGrid()

    for dur in range(durationiniter):
        for proton in protonTr:
            removal = proton.isitRemoved()
            if removal == False:

                proton.updatePos(dt)

                proton_coord = proton.getPosition()

                if proton_coord[0] < beamMargins[0]:
                    proton.removeproton()
                elif proton_coord[1] > sc_dim[1]:
                    proton.removeproton()
                elif proton_coord[1] < 0:
                    proton.removeproton()
                elif proton_coord[2] > sc_dim[2]:
                    proton.removeproton()
                elif proton_coord[2] < 0:
                    proton.removeproton()
                else:
                    passcount += 1

            if removal == False:
                smallerx = 0
                smallery = 0
                smallerz = 0

                for xcoord in sc_coordlist[0]:
                    if xcoord < proton_coord[0]:
                        smallerx = xcoord
                for ycoord in sc_coordlist[1]:
                    if ycoord < proton_coord[1]:
                        smallerx = ycoord
                for zcoord in sc_coordlist[2]:
                    if zcoord < proton_coord[2]:
                        smallerx = zcoord

                if sc_grid[(smallerx,smallery,smallerz)] == 1:
                    photonTr = proton.scintillate(dissipation,wavelength,quantum_efficiency,photonTr)
                    if math.sqrt(proton.velocityx**2+proton.velocityy**2+proton.velocityz**2) == 0:
                        proton.removeproton()
        
        if len(photonTr) > 0:
            for photon in photonTr:
                removal = photon.isitRemoved()
                if removal == False:

                    photon.updatePos(dt)

                    photon_coord = photon.getPosition()

                    if photon_coord[0] < beamMargins[0]:
                        photon.removeproton()
                    elif photon_coord[1] > sc_dim[1]:
                        photon.removeproton()
                    elif photon_coord[1] < 0:
                        photon.removeproton()
                    elif photon_coord[2] > sc_dim[2]:
                        photon.removeproton()
                    elif photon_coord[2] < 0:
                        photon.removeproton()
                    else:
                        passcount += 1

                if removal == False:
                    if photon_coord[0] > sc_dim[0]:
                        photondetectioncount += 1
                        detectedphotons.append(photon)
                        photon.removeproton()

    return photondetectioncount, detectedphotons

def grapher(trial_repetition, graph_nodesnum, mean_calcnum, numOfProtons, average_velocityx, average_velocityy, average_velocityz, velocity_variationx, velocity_variationy, velocity_variationz, beamMargins, sc_x, sc_y, sc_z, sc_resolution, dt, dissipation, wavelength, quantum_efficiency, durationiniter):
    yvalaggregated = []

    for graphing in range(trial_repetition):
        xvallist = []
        yvallist = []

        for concval in range(graph_nodesnum):
            concx = (concval*1000000)/graph_nodesnum
            resultlist = []
            counter = 0

            for trial in range(mean_calcnum):
                result = simulateTrial(numOfProtons, average_velocityx, average_velocityy, average_velocityz, velocity_variationx, velocity_variationy, velocity_variationz, beamMargins, concx, sc_x, sc_y, sc_z, sc_resolution, dt, dissipation, wavelength, quantum_efficiency, durationiniter)
                resultlist.append(result[0])

            for value in resultlist:
                counter += value

            yvallist.append((counter/len(resultlist)))
            xvallist.append(concx/10000)

        yvalaggregated.append(yvallist)

    plt.title("Concentrations Compared - Monte Carlo Stochastic Simulation (Quantum Efficiency: 20%)")
    plt.xlabel("Concentration (%)")
    plt.ylabel("Number of Detected Photons")
    for graphing in range(trial_repetition):
        plt.plot(xvallist,yvalaggregated[graphing])
    plt.show()

#use the grapher function to initiate trial
#the simulation requires IMMENSE computational power for realistic inputs, ONLY RUN ON HIGH-END DEVICES

#trial_repetition: the number of graphs in the final result
#graph_nodesnum: the number of nodes in each graph (For instance, graph_nodesnum = 10 would generate a graph where calculations were done with PMMA concentrations of 0%, 10% ... 90%)
#mean_calcnum: the number of trials conducted for each data point. results are averaged to obtain realistic results. ANY VALUE GREATER THAN 10 IS NOT RECOMMENDED
#numOfProtons: number of protons
#average_velocityx: average speed of protons in the x direction
#average_velocityy: average speed of protons in the y direction
#average_velocityz: average speed of protons in the z direction
#velocity_variationx: the standard deviation of the velocities in the x direction
#velocity_variationy: the standard deviation of the velocities in the y direction
#velocity_variationz: the standard deviation of the velocities in the z direction
#beamMargins: a tuple of form (x1,x2,y1,y2,z1,z2) dictating where the proton beam is generated. For instance, (0,10,0,10,0,10) would generate a proton beam where
#the protons are uniformly distributed over x = 0-10, y = 0-10, and z = 0-10.
#sc_x = dimension of the scintillator in the x direction (in cm)
#sc_y = dimension of the scintillator in the y direction (in cm)
#sc_z = dimension of the scintillator in the z direction (in cm)
#sc_resolution: the resolution of the scintillator, check the Scintillator class for more information.
#dt: the differential timeframe; dictates the exhaustiveness of the calculation. MUST BE GIVEN AS A screp() (SCIENTIFIC REPRESENTATION) object!
#dissipation: the approximate energy loss of a proton in PMMA per dt. MUST BE GIVEN AS A screp() (SCIENTIFIC REPRESENTATION) object!
#wavelength: the wavelength of the emitted photons. MUST BE GIVEN AS A screp() (SCIENTIFIC REPRESENTATION) object!
#quantum_efficiency: the quantum efficiency of the PMMA microparticles
#durationiniter: number of iterations per trial, serves as a control mechanism to prevent infinite loops and dictates the duration of the simulation
#durationiniter*dt = simulation duration

#grapher() ENTER PARAMETERS TO START SIMULATION