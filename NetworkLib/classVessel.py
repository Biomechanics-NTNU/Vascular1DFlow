import pprint
from copy import deepcopy

from classGrids import *
from classCompliance import *

class Vessel(object):
    '''
    Class representing a vessel in a vascular Network
    '''
    number = 0 
    def __init__(self, Id= None, name = None ,start = None, end = None, leftDaughter = None, rightDaughter = None, angleToMother = None, 
                 geom = None, length = None, radiusA = None, radiusAInterval = None, radiusB=None, radiusBInterval = None, N = None,NInterval = None,
                 comp = None, Pfunc = None, Ps = None, As = None, wallThickness = None, wallThicknessInterval = None, youngModulus = None,
                 youngModulusInterval = None, beta = None, betaInterval = None, my = None, rho = None, gamma = None,dlt = None, pref = None,
                 lengthInterval = None, myInterval= None, rhoInterval = None):
        '''
        Constructor
        '''
        ### Input properties
        
        # Topology properties
        self.Id = Id                                            # id of the vessel
        self.name = name                                        # name of the vessel
        self.start = start                                      # node id of starting node  
        self.end = end                                          # node id of ending node            
        self.leftDaughter = leftDaughter                        # id of left daugther vessel
        self.rightDaughter = rightDaughter                      # id of right daugther vessel
        self.angleToMother = angleToMother                      # angle to the x-axis of his mother vessel
                
        # GRID properties
        self.geom = geom                                        # geommetry type: 'Uni', 'Cone', 'Cons'
        self.length = length                                    # length of the vessel
        self.lengthInterval = lengthInterval
        self.radiusA = radiusA                                  # radius at vessel begin // in hole Vessel if radiusB = 0
        self.radiusAInterval = radiusAInterval
        self.radiusB = radiusB                                  # radius at vessel end
        self.radiusBInterval = radiusBInterval
        self.N = N                                              # number of gridpoints
        self.NInterval = NInterval                              
        
        # SOLID properties
        self.comp = comp                                        # complianceType: 'Exponential', 'Laplace', 'Reymond'
        self.Pfunc = Pfunc                                      # pressure function: 'True', 'False
        self.Ps = Ps                                            # pressure at state S
        self.As = As                                            # Area at state S
        self.wallThickness = wallThickness                      # wallthickness
        self.wallThicknessInterval = wallThicknessInterval
        self.youngModulus = youngModulus                        # YoungModulus
        self.youngModulusInterval= youngModulusInterval
        self.beta = beta                                        # beta value"
        self.betaInterval = betaInterval
        
        # FLUID properties 
        self.my = my                                            # blood viscosity
        self.myInterval = myInterval
        self.rho = rho                                          # blood density
        self.rhoInterval = rhoInterval 
        self.gamma = gamma                                      # velocity profile gammaection: 'True' 'False'
        self.dlt = dlt                                          # 
        self.pref = pref                                        # reference pressure
        
        ### Calculated properties
        # GRID properties
        self.z  = None                                          # axial coordinates of the vessel
        self.dz = None                                          # axial grid spacing
        self.A0 = None                                          # initial Area
        
        # SOLID properties   
        self.P0 = None    
        self.C  = None                                          # Compliance function C(P) for all nodes
        self.C_nID = None                                          # Compliance function C(P,NodeId) for specific Node with NodeID
        self.A  = None                                          # Area function A(P) for all nodes
        self.A_nID = None                                          # Area function A(P,NodeId) for specific Node with NodeID
        self.c  = None                                          # wave speed function c(rho,A,P) for all nodes
        self.c_nID = None                                          # wave speed function c(rho,A,P,NodeId) for specific Node with NodeID
        
        self.R = None                                           # resistance in the vessel
        self.alpha = None                                       # Womersley number
        
        self.quiet = False
        Vessel.number += 1
        
        
        
    def __del__(self):  
        '''
        Class Destructor
        '''
        if self.quiet == False: print "vessel ",self.Id," removed."
        Vessel.number -= 1 
        
        
    def initialize(self, globalFluid):
        '''
        Initialisation of the vessel.
        This method calculates and set up, the compliance, gird, fluid and resistance of the Vessel.
        '''
        # initialze fluid
        for key,value in globalFluid.iteritems():
            try:
                testArg = self.__getattribute__(key)
                if testArg == None:                
                    self.__setattr__(key,value)
            except:
                raise "ValueError: Fluid initialisation: %s, could not update varibale" %key
        # initialize gird
        
        try:
            grid = eval(self.geom)({'geom':self.geom,'l':self.length,'R':self.radiusA, 'r': self.radiusB, 'N': self.N})
            
            self.z,self.dz,self.A0 =  grid()
        except:
                raise ValueError("Error in Grid initialisation")
        # initialize compliance
        if True:
            
            if self.comp == 'Laplace':
                self.As  = np.ones(len(self.A0))*self.As
                self.wallThickness = np.ones(len(self.A0))*self.wallThickness
                self.youngModulus  = np.ones(len(self.A0))*self.youngModulus

            if self.comp == 'Laplace2':
                self.As = self.A0
                self.beta = np.ones(len(self.A0))*self.beta

            if self.comp == 'Exponential':
                self.As   = self.A0*self.Solid['As']
                self.beta = np.ones(len(self.A0))*self.Solid['beta']
                self.P0 = eval(self.comp)({'comp': self.comp,'Pfunc': self.Pfunc,'Ps': self.Ps,'As': self.As,'beta': self.beta}).P00(self.A0)

            solid = {'comp': self.comp,'Pfunc': self.Pfunc,'Ps': self.Ps,'As': self.As,'beta': self.beta,'H': self.wallThickness,'E': self.youngModulus}
            self.compliance = eval(self.comp)(solid,self.rho)
            if self.Pfunc:
                self.C  = self.compliance.C
                self.A  = self.compliance.A
                self.c  = self.compliance.c
                self.C_nID = self.compliance.C_
                self.A_nID = self.compliance.A_
                self.c_nID = self.compliance.c_
                
            else:
                self.C  = self.compliance.C0
                self.A  = self.compliance.A
                self.c  = self.compliance.c0
                self.C_nID = self.compliance.C0_
                self.A_nID = self.compliance.A_
                self.c_nID = self.compliance.c0_
        #except:
        #        raise ValueError("Error in Compliance set up")
            
        # set up Resistance
        
        # calculate Womersley number

    def getIntervalDictionary(self):
        '''
        retruns a dictionary containing all variable intervals
        the keys are the gammaesponding variable names
        '''
        intervalDictionary = {}
        intervalls = {"radiusA": self.radiusAInterval,"radiusB": self.radiusBInterval,"length": self.lengthInterval,
                               "wallThickness": self.wallThicknessInterval,"youngModulus": self.youngModulusInterval, 
                               "beta": self.betaInterval, "N":self.NInterval}
        for name, interval in intervalls.iteritems():
            if interval != None:
                intervalDictionary[name] = interval
        return deepcopy(intervalDictionary)
    
    
    def updateDataDict(self, Dict):
        '''
        updates the vessel data using a dictionary in from of 
	    dataDict = {'variableName': value}
        '''
        for key,value in Dict.iteritems():
            try:
                self.__getattribute__(key)
                self.__setattr__(key,value)
            except:
                ValueError("Wrong key: %s, could not update varibale" %key)
        
    def evaluateDaughters(self,vessels):
        '''
        find the daugthers of the vessel through the connectivity start and end
        '''
        if self.start != None and self.end != None:
            daughters = []
            # find the ids of the daughter vessel
            for vessel_i in vessels.itervalues():
                if vessel_i.start == self.end and self.end != None:
                    daughters.append(vessel_i.Id)
                    
            # place the left an the right
            if len(daughters) == 0:
                self.leftDaughter = None
                self.rightDaughter = None
            elif len(daughters) == 1:
                self.leftDaughter = daughters[0]
                self.rightDaughter = None
            elif len(daughters) == 2:
                if self.leftDaughter == daughters[0]:
                    self.rightDaughter = daughters[1]
                elif self.leftDaughter == daughters [1]:
                    self.rightDaughter = daughters[0]
                else:
                    self.leftDaughter = daughters[0]
                    self.rightDaughter = daughters[1]
            return True
        
        else: return False
        
        
    def printToConsole(self):
        '''
        writes all variables and their values to the console
        '''
        print "----------------"
        print "    Vessel %d"%self.Id,"\n"
        print " name:           ",self.name
        print " start:          ",self.start
        print " end:            ",self.end
        print " Ldaugther       ",self.leftDaughter
        print " Rdaugther       ",self.rightDaughter
        print " AngleToMother   ",self.angleToMother,"\n"
        print " geom            ",self.geom                                 
        print " length          ",self.length                                    
        print " lengthInterval  ",self.lengthInterval 
        print " radiusA         ",self.radiusA                                  
        print " radiusAInterval ",self.radiusAInterval 
        print " radiusB         ",self.radiusB                        
        print " radiusBInterval ",self.radiusBInterval  
        print " N               ",self.N
        print " N interval      ",self.NInterval,"\n"
        print " comp            ",self.comp                                       
        print " Pfunc           ",self.Pfunc                                  
        print " Ps              ",self.Ps                                           
        print " As              ",self.As                                           
        print " wallThickness   ",self.wallThickness                      
        print " wallThick.Int.  ",self.wallThicknessInterval 
        print " youngModulus    ",self.youngModulus          
        print " youngMod. Int.  ",self.youngModulusInterval
        print " beta            ",self.beta                         
        print " betaInterval    ",self.betaInterval,"\n"
        print " my              ",self.my                                        
        print " rho             ",self.rho                                           
        print " gamma           ",self.gamma
        print " dlt             ",self.dlt
        print " pref            ",self.pref,"\n"
        print " z               ",self.z                                       
        print " dz              ",self.dz                                     
        print " A0              ",self.A0,"\n"
              
        print " R               ",self.R ,"\n"
