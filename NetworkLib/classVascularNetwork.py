
import sys,os
# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath( __file__ ) )
sys.path.append(cur+'/../NetworkLib')

from classVessel import Vessel
from copy import deepcopy

import pprint


class VascularNetwork(object):
    '''
    Class representing a vascular Network
    The vascular network consits out of vessels defined in classVessel::Vessel()
    Additional Topologie, BoundaryConditions and the SimulationContext are saved.
    '''
    def __init__(self, quiet=False):
        '''
        Constructor
        '''
        # name of the network
        self.name = 'vascularNetwork'
        
        # dictionary with containing all vessel data
        # key = vessel id; value = vessel::Vessel()
        self.vessels = {}
        self.root = []
        self.ends = []
        
        ## Dictionary with all cumultative resistances
        self.Rcum = {}
        
        # Boundary Conditions in form of {vessel.Id : {type: values=[] }}
        self.boundaryConditions = {}
        self.boundaryConditionIntervals = {}
        
        # dictionary containing the global fluid data if defined
        self.globalFluid = {'my': None,'rho': None,'gamma': None,'dlt':None ,'pref':None}
        
        # dictionary containg the simulation Context
        self.simulationContext = {'totalTime': 1.0, 'CFL': 0.95}
        
        # dirctonary containing the solver calibration 
        self.networkSolver = {'NumScheme':'MacCormack','EqSystem':'02','CharSystem':0}
     
        self.quiet = quiet
        
        # tree traverse list 
        self.calcList = []
     
    def __del__(self):  
        '''
        Class Destructor
        '''
        for key in self.vessels.keys():
            ves = self.vessels.pop(key)
            ves.quiet = self.quiet
            del ves
        if self.quiet == False: print "VascularNetwork removed."
        
    
    # all classes concerning vessel
    def addVessel(self, Id = None, dataDict = False):
        '''
        adds vessel to the Network
	if no id, a random id is choosen
	if no DataDict, no values are assigned
	'''
        # set id to 1 + highest id of existing vessels
        if Id == None: 
            try: Id = max(self.vessels.keys())+1
            except: Id = 0
            
        # check Id
        if Id not in self.vessels:
                # create vessel with given variables
                vessel = Vessel(Id = Id , name = ('vessel_'+str(Id)))
                # set vesselData if available
                if dataDict:
                    vessel.updateDataDict(dataDict)
                # add vessel to network
                self.vessels[vessel.Id] = vessel
        # raise error if Id is set doubled
        else:
            raise ValueError("vessel Id =%d exists already! Could not add vessel" %vessel.Id)
 
    def deleteVessel(self, inputID):
        '''
        Remove vessel from network and delete it
        '''
        if inputID in self.vessels:
            del self.vessels[inputID]
        else:
            raise ValueError("vessel with Id =%d does not exist! Could not remove vessel" %inputID)
     
    def updateNetwork(self,updateDict):
        '''
        Update vascular Network with an Dictionary: updateDict
        
        updateDict = {'networkData': networkDataDict, vesselData':vesselDataDict}

			networkDataDict = {' boundaryConditions': netBC,'globalFluid':netFluid,
                                            'simulationContext': netSimCon }
			vesselDataDict = { vessel.id : DataDict}
        '''
        if 'networkData' in updateDict: 
            self.boundaryConditions.update((updateDict['networkData'])['boundaryConditions'])
            self.globalFluid.update((updateDict['networkData'])['globalFluid']) 
            self.simulationContext.update((updateDict['networkData'])['simulationContext'])
            self.networkSolver.update((updateDict['networkData'])['networkSolver'])
        
        if 'vesselData' in updateDict:
            for Id, vesselData in (updateDict['vesselData']).iteritems():
                if Id in self.vessels:
                    self.vessels[Id].updateDataDict(vesselData)
                else:
                    self.addVessel(Id, vesselData)    
                
    def showVessels(self):
        '''
        writes Network properties (without vesselData) to console
        '''
        print " Vessels in Network:"
        for vessel_i in self.vessels.itervalues():
            vessel_i.printToConsole()
    
    def showNetwork(self):
        '''
        writes the Vesseldata for each vessel to console (calls printToConsole() from each vessel)
	    '''
        print "-------------------"
        print " vascularNetwork ",self.name,"\n"
        print " boundaryConditions:"
        pprint.pprint(self.boundaryConditions)
        print " gloabl Fluid properties:"
        pprint.pprint(self.globalFluid)
        print " simulation Context"
        pprint.pprint(self.simulationContext)
        print " network Solver"
        pprint.pprint(self.networkSolver)
        
        
    def initialize(self):
        for vessel_i in self.vessels.itervalues():
            vessel_i.initialize(self.globalFluid)
              
        for id,bcs in self.boundaryConditions.iteritems():
            for bc in bcs:
                if bc.name in ['_Ue','Ue']: bc.area = self.vessels[id].A0
                
    def evaluateConnections(self):
        '''
        This functions evaluates the connection between vessels: 
        if nodes are defined the daughters are set if not defined;
        if daughters are defined the nodes are set;
        In addition the root vessel and the open ends are found and saved in
        the self.root and self.ends.
        Futher are the position of the boundaryConditions determined and set.
        '''   
        ids = []
        daughters = []
        noDaughters = []
        startEnds = []
        for vessel_i in self.vessels.itervalues():
            startEnds.append(vessel_i.evaluateDaughters(self.vessels))  
            if vessel_i.leftDaughter != None: 
                daughters.append(vessel_i.leftDaughter)
                if vessel_i.rightDaughter != None: daughters.append(vessel_i.rightDaughter)
            else:
                noDaughters.append(vessel_i.Id)
            ids.append(vessel_i.Id)
        if daughters != []:
            for daugther in daughters:
                if daugther in ids:
                    ids.remove(daugther)
        self.root = ids
        self.ends = noDaughters
        
                
        # set start and end nodes for the network,parsing through the network as a binary tree
        
        if False in startEnds:
            self.calcList = []
            # set node count
            count = 0
            # vessels with non visualized daughters:
            viz = []
            #find root
            root = self.root[0]
            # set root start node
            self.vessels[root].start = count
            #increase count
            count = count + 1
            #set root end node
            self.vessels[root].end = count
            #increase count
            count = count + 1
            # add root to the viz vessels if root has daughters:
            if self.vessels[root].leftDaughter != None:
                viz.append(root)
            else:
                self.calcList.append([root])
                
            # loop through tree until all daughters are conected
            while len(viz) != 0:                
                # get the mother vessel (already added) and add its daughters
                motherVessel = viz.pop(0)
                calcListTemp=[motherVessel]
                # find left daughter
                leftDaughter = self.vessels[motherVessel].leftDaughter
                calcListTemp.append(leftDaughter)
                # set start node
                self.vessels[leftDaughter].start = self.vessels[motherVessel].end
                #set end node
                self.vessels[leftDaughter].end = count
                #increase count
                count = count + 1
                
                calcListTemp2= None
                # check if leftDaughter has also daughters which should be visualized
                if self.vessels[leftDaughter].leftDaughter != None:
                    viz.append(leftDaughter)
                else:
                    calcListTemp2 = [leftDaughter]
                    
                # find right daughter
                if self.vessels[motherVessel].rightDaughter != None:
                    rightDaughter = self.vessels[motherVessel].rightDaughter
                    calcListTemp.append(rightDaughter)
                    # set start node
                    self.vessels[rightDaughter].start = self.vessels[motherVessel].end
                    #set end node
                    self.vessels[rightDaughter].end = count
                    #increase count
                    count = count + 1
                    
                    # check if rightDaughter has also daughters which should be visualized
                    if self.vessels[rightDaughter].leftDaughter != None:
                        viz.append(rightDaughter)
                    else:
                        calcListTemp2 = [rightDaughter]
                        
                self.calcList.append(calcListTemp)
                if calcListTemp2 != None:self.calcList.append(calcListTemp2)
        print self.calcList
        
        if self.boundaryConditions != {}:
            ## Set position of boundary conditions
            # check position if one Vessel
            if len(self.vessels) == 1:
                id = self.boundaryConditions.keys()[0]
                if id != self.root[0]: print "Error Wrong Root found"
                for bc in self.boundaryConditions[id]:
                    if '_' not in bc.name: bc.setPosition(0)
                    else: bc.setPosition(-1)
            else:
                for id, bcs in self.boundaryConditions.iteritems():
                    if id in self.root:
                        for bc in bcs:
                            bc.setPosition(0)
                    elif id in self.ends:
                        for bc in bcs:
                            bc.setPosition(-1)
                            
        
    
    def calculateNetworkResistance(self):
        '''
        This function travers the network tree and calculates the 
        cumultative system resistances Rcum for each vessel in the Network.
        '''
        ## firstTraverse the tree and create a back traverse list
        calcList = []
        viz = []
        #find root
        root = self.root[0]        
        # add root to the viz vessels if root has daughters:
        if self.vessels[root].leftDaughter != None:
            viz.append(root)
            self.Rcum[root] = None
        else:
            self.Rcum[root] = self.vessels[root].R
        # loop through tree until all daughters are conected
        while len(viz) != 0:
            # get the mother vessel (already added) and add its daughters
            motherVessel = viz.pop(0)
            
            # find left daughter
            leftDaughter  = self.vessels[motherVessel].leftDaughter
            rightDaughter = self.vessels[motherVessel].rightDaughter
            #append the mother to the calc list
            curCalcList = [motherVessel]
            
            if leftDaughter != None:
                #append the left daughter to the calc list
                curCalcList.append(leftDaughter)
                self.Rcum[leftDaughter] = None
                if rightDaughter != None:
                    #append the left daughter to the calc list
                    curCalcList.append(rightDaughter)
                    self.Rcum[rightDaughter] = None
                else:
                    curCalcList.append(None)
            
                # check if leftDaughter has also daughters 
                if self.vessels[leftDaughter].leftDaughter != None:
                    viz.append(leftDaughter)
                else:
                    #calcList.extend([leftDaughter,-1,-1])
                    self.Rcum[leftDaughter] = self.vessels[leftDaughter].R
                    
                if rightDaughter != None:
                    # check if rightDaughter has also daughters 
                    if self.vessels[rightDaughter].leftDaughter != None:
                        viz.append(rightDaughter)
                    else:
                        #calcList.extend([rightDaughter,-1,-1])
                        self.Rcum[rightDaughter] = self.vessels[rightDaughter].R 
            calcList.append(curCalcList)

        ## Loop through the caclList and calculate the resistance
        for toCalc in reversed(calcList):
            
            if toCalc[2] != None:
                ## Add calculation Method
                #Rcum[toCalc[0]] = ''.join([self.vessels[toCalc[0]].R,'+(',Rcum[toCalc[1]],'^(-1)+',Rcum[toCalc[2]],'^(-1))^(-1)'])
                self.Rcum[toCalc[0]] = ''.join(['[(',self.Rcum[toCalc[1]],'+',self.Rcum[toCalc[2]],')]'])
            else:
                ## Add calculation Method
                self.Rcum[toCalc[0]] = ''.join(['[',self.Rcum[toCalc[1]],']'])
                #Rcum[toCalc[0]] = ''.join([self.vessels[toCalc[0]].R,'+',Rcum[toCalc[1]]])
        pprint.pprint(self.Rcum[self.root[0]])   
        

        