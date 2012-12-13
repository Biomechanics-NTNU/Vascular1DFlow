import numpy as np
import scipy as sc
from itertools import islice

import sys,os

# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath( __file__ ) )
sys.path.append(cur+'/NetworkLib')

from classBoundarys import Boundary_02

from classDerivatives import *
from classSystemEquations import *
from classConnections import *

class FlowSolver(object):
    
    def __init__(self,vascularNetwork, quiet=False):
        '''
        Constructor       
        '''
        if vascularNetwork == None: raise Error("No vascularNetwork given!")
        
        # the vascular network to solve
        self.vascularNetwork = vascularNetwork
        # the boundarys of the network { vesselID : [<instance>::classBoundary_02(Characteristics.py), .. ]}
        # 1 boundary for each start/end-vessel except if only 1 vessel in the network
        self.boundarys = {}
        # the system Equations of each vessel { vesselID : <instance>::classSystemEquations(SystemEquations) }
        self.systemEquations = {}
        # the connections of the Network { motherVesselID : <instance>::classConnections, ...}
        self.connections  = {}                
        # time step
        self.dt = None
        # number of Timesteps
        self.Tsteps = None
        # total Simulation time ??? need this?! is saved in vascularNetwork??
        self.totalTime = None
        
        # Set div output 
        self.output = {}
        
        # solution of the system over time {vesselID: [ [solution at N nodes]<-one array for each timePoint , ...  ]  }
        self.P  = {}
        self.Q  = {}
        self.A  = {}
        self.d0 = {} # not transient
        # solution of the predictorStep // earlier called _P
        # {vesselID: [solution at N nodes]  }
        self.P_pre  = {}
        self.Q_pre  = {}
        self.A_pre  = {}
        self.d0_pre = {}
        
        # Set derivatives functions
        self.df = DeltaForward
        self.db = DeltaBackward
        
        # characteristic system variable (0/1)
        self.charSys = self.vascularNetwork.networkSolver['CharSystem']
        # euqSys A = A0 or A = A(P)
        self.equSys = self.vascularNetwork.networkSolver['EqSystem']
        print 'self.equSys',self.equSys
         # 4. Define the output of A, dependend on the characteristic system 0.1
        if self.equSys == '02':
            self.AFunction = self.AFunctionSys0
        elif self.equSys == '2': 
            self.AFunction = self.AFunctionSys1
        else:
            raise Error("CharSys not properly defined!")
        # initialize system
        # refresh/update vascularNetwork (just to be sure it is done)
        vascularNetwork.evaluateConnections()
        vascularNetwork.initialize()
        
        self.initializeTimeVariables()
        # initialize boundarys
        self.initializeBoundarys()
        # initialize system Equations
        self.initializeSystemEquations()
        # initialize Connections
        self.initializeConnections()
        # initialize solution
        self.initializeSolutionMatrices()
        # feedback
        if quiet==False: self.initOutput()

        #def solving sheme
        self.solvingScheme = self.vascularNetwork.networkSolver['NumScheme']
        ## ADD EVALUATE solver schmes! each solver scheme get an one function
        
        #self.solve = self.MacCormack_2_tree
        self.solve = self.MacCormack_2
        #self.solve = self.MacCormack_2InB
    '''       
    ########################################################################################
    # initialisation Methods
    ########################################################################################
    '''
    def calcTimeStep(self,dz,c,CFL):
        return (CFL*dz)/c

    def initializeTimeVariables(self):
        '''
        initialize time variable dt and Tstep
        '''
        self.totalTime = self.vascularNetwork.simulationContext['totalTime']
        
        dt_min,dz_min,c_max,gridNodens = [],[],[],[]
        #create waveSpeed Log file
        logfileData = {}
        
        dt = self.totalTime
        for vessel_i in self.vascularNetwork.vessels.itervalues():
        # Calculate time variables
            
            A0_max = max(vessel_i.A0)
            c_high = vessel_i.c(A0_max,vessel_i.pref)
            dz_low = min(vessel_i.dz)
            dt = self.calcTimeStep(dz_low,c_high,self.vascularNetwork.simulationContext['CFL']) 
            
            c_max = np.append(c_max,c_high)
            dt_min = np.append(dt_min,dt)
            dz_min = np.append(dz_min,dz_low)
            gridNodens = np.append(gridNodens,vessel_i.N)
            
            logfileData[vessel_i.Id] = [max(c_high),min(c_high),min(dt),vessel_i.dz,vessel_i.N]
            
        # Set time variables 
        self.dt = min(dt_min)
        self.Tsteps = int(round(self.totalTime/self.dt))
        
        
        self.output['dz_min']     = min(dz_min)
        self.output['c_min']      = min(c_max)
        self.output['c_max']      = max(c_max)
        self.output['gridNodens'] = sum(gridNodens)
        
        self.output['CFLcorrect'] = []
        
        logfile = open(str(cur+'/../'+'LOGcurrentWaveSpeed.txt'),'wb')
        logfile2 = open(str(cur+'/../'+'LOGproposedGrid.txt'),'wb')
        CFL = self.vascularNetwork.simulationContext['CFL']
        for vesselT,data in logfileData.iteritems():
            #number of deltaX
            Nnew = int((sum(data[3])*CFL/(self.dt*data[0])))
            
            L = sum(data[3])
            dz_new = L/Nnew
            dt_new = self.calcTimeStep(dz_new,data[1],CFL)
            while dt_new < self.dt:
                Nnew = Nnew-1
                dz_new = L/Nnew
                dt_new = self.calcTimeStep(dz_new,data[1],CFL)
            if Nnew == data[4]: logfile.write(''.join(['vessel ',str(vesselT), ' c_max: %2.3f'%(data[0]),'   dt (ms) %2.6f'%(data[2]*1.0E3),' | res CFL: %2.3f'%(data[0]*self.dt/min(data[3])),' || already best N', '\n']))
            else: 
                logfile.write(''.join(['vessel ',str(vesselT), ' c_max: %2.3f'%(data[0]),'   dt (ms) %2.6f'%(data[2]*1.0E3),' | res CFL: %2.3f'%(data[0]*self.dt/min(data[3])),' || prop: N%2.0f'%Nnew,'    dN %.0f'%(Nnew-data[4]),'    dtNew %.4f'%(dt_new*1.e3),'   CFL: %2.3f'%(data[0]*self.dt/dz_new), '\n']))
                self.output['CFLcorrect'].append(' '.join([str(vesselT),str(data[4]),str(Nnew)]))
            logfile2.write(''.join([str(int(Nnew+1)),'\n']))
        logfile.close()
        logfile2.close()    
         
    def initializeBoundarys(self):
        '''
        initialize boundarys
        '''
        if len(self.vascularNetwork.vessels) == 1:
            id = self.vascularNetwork.root[0]
            bcList0 = []
            bcList1 = []
            for bc in self.vascularNetwork.boundaryConditions[id]:
                print 'here', bc
                if bc.position == 0:     bcList0.append(bc)
                elif bc.position == -1:  bcList1.append(bc)
            self.boundarys[id] = [Boundary_02(bcList0,self.vascularNetwork.vessels[id].A_nID,self.equSys),
                                  Boundary_02(bcList1,self.vascularNetwork.vessels[id].A_nID,self.equSys)]  
            self.output['BndrNR'] = 2
        else:
            for id,boundaryConditions in self.vascularNetwork.boundaryConditions.iteritems():
                self.boundarys[id] = [Boundary_02(boundaryConditions,self.vascularNetwork.vessels[id].A_nID,self.equSys)]
            self.output['BndrNR'] = len(self.boundarys)

    def initializeSystemEquations(self):
        '''
        initialize system Equations
        '''
        for id,vessel in self.vascularNetwork.vessels.iteritems():
            self.systemEquations[id] = System(vessel,self.charSys)
    
    def initializeConnections(self):
        '''
        initialize Connections of the network
        by traversing the network tree
        '''
        vcalc = []
        root = self.vascularNetwork.root[0]
        if self.vascularNetwork.vessels[root].leftDaughter != None:
            vcalc.append(root)
        # loop through tree until the end of it is reached
        while len(vcalc) != 0:
            motherVessel = vcalc.pop(0)
            leftDaughter = self.vascularNetwork.vessels[motherVessel].leftDaughter
            if self.vascularNetwork.vessels[leftDaughter].leftDaughter != None:
                vcalc.append(leftDaughter)
            if self.vascularNetwork.vessels[motherVessel].rightDaughter != None:
                rightDaughter = self.vascularNetwork.vessels[motherVessel].rightDaughter
                if self.vascularNetwork.vessels[rightDaughter].leftDaughter != None:
                    vcalc.append(rightDaughter)
                # initialize bifucration between mother and left&right daughters
                self.connections[motherVessel] = Connection([self.vascularNetwork.vessels[motherVessel]],
                                                             [self.systemEquations[motherVessel]],
                                                             [self.vascularNetwork.vessels[leftDaughter],
                                                              self.vascularNetwork.vessels[rightDaughter]],
                                                             [self.systemEquations[leftDaughter],
                                                              self.systemEquations[rightDaughter]],
                                                              self.equSys)
                
            else:
                # initialize connection between mother and left&right daughters
                self.connections[motherVessel] = Connection([self.vascularNetwork.vessels[motherVessel]],
                                                             [self.systemEquations[motherVessel]],
                                                             [self.vascularNetwork.vessels[leftDaughter]],
                                                             [self.systemEquations[leftDaughter]],
                                                             self.equSys)
    
    def initializeSolutionMatrices(self):
        '''
        initialize solution matrices
        '''
        for id,vessel in self.vascularNetwork.vessels.iteritems():
            
            self.P[id]    = np.zeros((self.Tsteps,vessel.N))
            self.P[id][0] = np.ones((1,vessel.N))*vessel.pref
            self.Q[id]    = np.zeros((self.Tsteps,vessel.N))
            self.A[id]    = np.zeros((self.Tsteps,vessel.N))
            self.A[id][0] = np.ones((1,vessel.N))*vessel.A0
            self.d0[id]   = [np.zeros((2)),np.zeros((2))]
            # do we need this below here?
            self.P_pre[id]  = np.ones((vessel.N))
            self.Q_pre[id]  = np.ones((vessel.N))
            self.A_pre[id]  = np.ones((vessel.N))
            self.d0_pre[id] = [np.ones((2)),np.ones((2))]
        
    def initOutput(self):
        '''
        initialize solution matrices
        '''
        print '====================================='
        print '___________Time variables ____________'
        print '%-20s %2.1f' % ('totaltime (sec)',self.totalTime)
        print '%-20s %2.3f' % ('dt (ms)',self.dt*1.0E3)
        print '%-20s %4d' % ('Tsteps',self.Tsteps)
        print '___________Div variables _____________'
        print '%-20s %2.1f' % ('CFL max',self.vascularNetwork.simulationContext['CFL'])
        print '%-20s %2.1f' % ('dz min (mm)',self.output['dz_min']*1.0E3)
        print '%-20s %2.1f' % ('c min (m/s)',self.output['c_min'])
        print '%-20s %2.1f' % ('c max (m/s)',self.output['c_max'])
        print '%-20s %4d' % ('Grid nodens',self.output['gridNodens'])
        print '%-20s %4d' % ('BC',self.output['BndrNR'])
        print '___________CFL-correction_____________'
        for CFLcorr in self.output['CFLcorrect']:
            print '%-20s %4s' % ('id N-old N-prop.',CFLcorr)
            
    '''
    ########################################################################################
    # Solver Methods:
    #
    #    MacCormack_2
    #    MacCormack_2InB
    #    MacCormack_2_tree
    #
    ########################################################################################
    '''
    
    def MacCormack_2(self):
        '''
        MacCormack_2 solver method with forward-euler time steping,
        Using either Characteristic system 0 or 1 as defined in the XML-file.
        
        This method is solving the system by looping through all vessels, then boundarys, conncetions, bifucations?!
        '''
        print "Solving system"
        for n in range(self.Tsteps-1):
            '''Predictor Step'''
            for id,vessel in self.vascularNetwork.vessels.iteritems():
                
                self.systemEquations[id].updateSystem(self.P[id][n],self.Q[id][n],self.A[id][n])
                dPdz  = self.df(self.P[id][n],self.vascularNetwork.vessels[id].dz)
                dQdz  = self.df(self.Q[id][n],self.vascularNetwork.vessels[id].dz) 
                dQ2dz = self.df(pow(self.Q[id][n],2.)/self.A[id][n],self.vascularNetwork.vessels[id].dz)
                            
                
                self.P_pre[id] = self.P[id][n] - (self.systemEquations[id].M[0][1]*dQdz)*self.dt
                self.Q_pre[id] = self.Q[id][n] - (self.systemEquations[id].M[1][0]*dPdz + self.systemEquations[id].M[1][1]*dQ2dz - self.systemEquations[id].B[1] )*self.dt
                self.A_pre[id] = self.AFunction(id,self.A[id][n],self.P_pre[id])
                
            for id,boundarys in self.boundarys.iteritems():
                for boundary in boundarys:
                    position = boundary.position
                    '''Solve start - boundary Condition of the vessel ''' 
                    self.P_pre[id][position],self.Q_pre[id][position],self.A_pre[id][position],self.d0_pre[id] =\
                                         boundary(self.P[id][n],
                                                  self.Q[id][n],
                                                  self.A[id][n],
                                                  self.d0[id],
                                                  self.vascularNetwork.vessels[id].z,
                                                  n,
                                                  self.dt,
                                                  self.systemEquations[id].L,
                                                  self.systemEquations[id].R,
                                                  self.systemEquations[id].LAMBDA)
                                         
            for id,connection in self.connections.iteritems():
                leftDaughter = self.vascularNetwork.vessels[id].leftDaughter
                rightDaughter = self.vascularNetwork.vessels[id].rightDaughter
                
                if rightDaughter == None: 
                
                    self.P_pre[id][-1],self.Q_pre[id][-1],self.A_pre[id][-1],\
                    self.P_pre[leftDaughter][0],self.Q_pre[leftDaughter][0],self.A_pre[leftDaughter][0] = connection([self.P[id][n],self.P[leftDaughter][n]],
                                                                                                                     [self.Q[id][n],self.Q[leftDaughter][n]],
                                                                                                                     [self.A[id][n],self.A[leftDaughter][n]],
                                                                                                                      self.dt)
                else:
                    '''solve bifurcation/connection between mother and left&right daughters'''
                    self.P_pre[id][-1],self.Q_pre[id][-1],self.A_pre[id][-1],\
                    self.P_pre[leftDaughter][0],self.Q_pre[leftDaughter][0],self.A_pre[leftDaughter][0],\
                    self.P_pre[rightDaughter][0],self.Q_pre[rightDaughter][0],self.A_pre[rightDaughter][0] = connection([self.P[id][n],self.P[leftDaughter][n],self.P[rightDaughter][n]],
                                                                                                                        [self.Q[id][n],self.Q[leftDaughter][n],self.Q[rightDaughter][n]],
                                                                                                                        [self.A[id][n],self.A[leftDaughter][n],self.A[rightDaughter][n]],
                                                                                                                         self.dt)
            '''Corrector Step'''
            for id,vessel in self.vascularNetwork.vessels.iteritems():
                
                self.systemEquations[id].updateSystem(self.P_pre[id],self.Q_pre[id],self.A_pre[id])
            
                dPdz  = self.db(self.P_pre[id],self.vascularNetwork.vessels[id].dz)
                dQdz  = self.db(self.Q_pre[id],self.vascularNetwork.vessels[id].dz)     
                dQ2dz = self.db(pow(self.Q_pre[id],2.)/self.A_pre[id],self.vascularNetwork.vessels[id].dz)
                
                self.P[id][n+1] = (self.P[id][n] + self.P_pre[id] - (self.systemEquations[id].M[0][1]*dQdz)*self.dt )/2.0
                self.Q[id][n+1] = (self.Q[id][n] + self.Q_pre[id] - (self.systemEquations[id].M[1][0]*dPdz + self.systemEquations[id].M[1][1]*dQ2dz - self.systemEquations[id].B[1] )*self.dt )/2.0
                self.A[id][n+1] = self.AFunction(id,self.A_pre[id],self.P[id][n+1])
            
            for id,boundarys in self.boundarys.iteritems():
                for boundary in boundarys:
                    position = boundary.position
                    '''Solve start - boundary Condition of the vessel ''' 
                    P_corr,Q_corr,self.A[id][n+1][position],self.d0[id] =\
                                         boundary(self.P_pre[id],
                                                  self.Q_pre[id],
                                                  self.A_pre[id],
                                                  self.d0_pre[id],
                                                  self.vascularNetwork.vessels[id].z,
                                                  n,
                                                  self.dt,
                                                  self.systemEquations[id].L,
                                                  self.systemEquations[id].R,
                                                  self.systemEquations[id].LAMBDA)
                                         
                    self.P[id][n+1][position] = P_corr
                    self.Q[id][n+1][position] = Q_corr
            
            for id,connection in self.connections.iteritems():
                    leftDaughter = self.vascularNetwork.vessels[id].leftDaughter
                    rightDaughter = self.vascularNetwork.vessels[id].rightDaughter
                    if rightDaughter == None: 
                        '''solve bifurcation/connection between mother and left&right daughters'''
                        self.P[id][n+1][-1],self.Q[id][n+1][-1],self.A[id][n+1][-1],\
                        self.P[leftDaughter][n+1][0],self.Q[leftDaughter][n+1][0],self.A[leftDaughter][n+1][0] = connection([self.P_pre[id],self.P_pre[leftDaughter]],
                                                                                                                            [self.Q_pre[id],self.Q_pre[leftDaughter]],
                                                                                                                            [self.A_pre[id],self.A_pre[leftDaughter]],
                                                                                                                             self.dt)
                    else:
                        '''solve bifurcation/connection between mother and left&right daughters'''
                        self.P[id][n+1][-1],self.Q[id][n+1][-1],self.A[id][n+1][-1],\
                        self.P[leftDaughter][n+1][0],self.Q[leftDaughter][n+1][0],self.A[leftDaughter][n+1][0],\
                        self.P[rightDaughter][n+1][0],self.Q[rightDaughter][n+1][0],self.A[rightDaughter][n+1][0] = connection([self.P_pre[id],self.P_pre[leftDaughter],self.P_pre[rightDaughter]],
                                                                                                                               [self.Q_pre[id],self.Q_pre[leftDaughter],self.Q_pre[rightDaughter]],
                                                                                                                               [self.A_pre[id],self.A_pre[leftDaughter],self.A_pre[rightDaughter]],
                                                                                                                                self.dt)
        return self.P,self.Q,self.A
    
    
    
    def MacCormack_2InB(self):
        '''
        MacCormack_2 solver method with forward-euler time steping,
        Using either Characteristic system 0 or 1 as defined in the XML-file.
        
        This method is solving the system by looping through all vessels, then boundarys, conncetions, bifucations?!
        '''
        print "Solving system"
        for n in range(self.Tsteps-1):
            '''Predictor Step'''
            for id,vessel in self.vascularNetwork.vessels.iteritems():
                
                dz = vessel.dz
                self.systemEquations[id].updateSystem(self.P[id][n],self.Q[id][n],self.A[id][n])
                dPdz  = DeltaForwardInnerPart(self.P[id][n],dz)
                dQdz  = DeltaForwardInnerPart(self.Q[id][n],dz) 
                dQ2dz = DeltaForwardInnerPart(pow(self.Q[id][n],2.)/self.A[id][n],dz)
                         
                self.P_pre[id][1:-1] = self.P[id][n][1:-1] - (self.systemEquations[id].M[0][1][1:-1]*dQdz)*self.dt
                self.Q_pre[id][1:-1] = self.Q[id][n][1:-1] - (self.systemEquations[id].M[1][0][1:-1]*dPdz + self.systemEquations[id].M[1][1][1:-1]*dQ2dz - self.systemEquations[id].B[1][1:-1] )*self.dt
                self.A_pre[id][1:-1] = self.AFunction(id,self.A[id][n],self.P_pre[id])[1:-1]
                
                if id in self.boundarys:
                    for boundary in self.boundarys[id]:
                        position = boundary.position
                        '''Solve start - boundary Condition of the vessel ''' 
                        self.P_pre[id][position],self.Q_pre[id][position],self.A_pre[id][position],self.d0_pre[id] =\
                                             boundary(self.P[id][n],
                                                      self.Q[id][n],
                                                      self.A[id][n],
                                                      self.d0[id],
                                                      vessel.z,
                                                      n,
                                                      self.dt,
                                                      self.systemEquations[id].L,
                                                      self.systemEquations[id].R,
                                                      self.systemEquations[id].LAMBDA)
                                         
                if id in self.connections:
                    leftDaughter = vessel.leftDaughter
                    rightDaughter = vessel.rightDaughter
                    
                    if rightDaughter == None: 
                    
                        self.P_pre[id][-1],self.Q_pre[id][-1],self.A_pre[id][-1],\
                        self.P_pre[leftDaughter][0],self.Q_pre[leftDaughter][0],self.A_pre[leftDaughter][0] = self.connections[id]([self.P[id][n],self.P[leftDaughter][n]],
                                                                                                                         [self.Q[id][n],self.Q[leftDaughter][n]],
                                                                                                                         [self.A[id][n],self.A[leftDaughter][n]],
                                                                                                                          self.dt)
                    else:
                        '''solve bifurcation/connection between mother and left&right daughters'''
                        self.P_pre[id][-1],self.Q_pre[id][-1],self.A_pre[id][-1],\
                        self.P_pre[leftDaughter][0],self.Q_pre[leftDaughter][0],self.A_pre[leftDaughter][0],\
                        self.P_pre[rightDaughter][0],self.Q_pre[rightDaughter][0],self.A_pre[rightDaughter][0] = self.connections[id]([self.P[id][n],self.P[leftDaughter][n],self.P[rightDaughter][n]],
                                                                                                                            [self.Q[id][n],self.Q[leftDaughter][n],self.Q[rightDaughter][n]],
                                                                                                                            [self.A[id][n],self.A[leftDaughter][n],self.A[rightDaughter][n]],
                                                                                                                             self.dt)
            '''Corrector Step'''
            for id,vessel in self.vascularNetwork.vessels.iteritems():
                
                self.systemEquations[id].updateSystem(self.P_pre[id],self.Q_pre[id],self.A_pre[id])
            
                dz = vessel.dz
                dPdz  = DeltaBackwardInnerPart(self.P_pre[id],dz)
                dQdz  = DeltaBackwardInnerPart(self.Q_pre[id],dz)     
                dQ2dz = DeltaBackwardInnerPart(pow(self.Q_pre[id],2.)/self.A_pre[id],dz)
                
                self.P[id][n+1][1:-1] = (self.P[id][n][1:-1] + self.P_pre[id][1:-1] - (self.systemEquations[id].M[0][1][1:-1]*dQdz)*self.dt )/2.0
                self.Q[id][n+1][1:-1] = (self.Q[id][n][1:-1] + self.Q_pre[id][1:-1] - (self.systemEquations[id].M[1][0][1:-1]*dPdz + self.systemEquations[id].M[1][1][1:-1]*dQ2dz - self.systemEquations[id].B[1][1:-1] )*self.dt )/2.0
                self.A[id][n+1][1:-1] = self.AFunction(id,self.A_pre[id],self.P[id][n+1])[1:-1]
            
                if id in self.boundarys:
                    for boundary in self.boundarys[id]:
                        position = boundary.position
                        '''Solve start - boundary Condition of the vessel ''' 
                        self.P[id][n+1][position],self.Q[id][n+1][position],self.A[id][n+1][position],self.d0[id] =\
                                             boundary(self.P_pre[id],
                                                      self.Q_pre[id],
                                                      self.A_pre[id],
                                                      self.d0_pre[id],
                                                      vessel.z,
                                                      n,
                                                      self.dt,
                                                      self.systemEquations[id].L,
                                                      self.systemEquations[id].R,
                                                      self.systemEquations[id].LAMBDA)
            
            
                if id in self.connections:
                        leftDaughter = vessel.leftDaughter
                        rightDaughter = vessel.rightDaughter
                        if rightDaughter == None: 
                            '''solve bifurcation/connection between mother and left&right daughters'''
                            self.P[id][n+1][-1],self.Q[id][n+1][-1],self.A[id][n+1][-1],\
                            self.P[leftDaughter][n+1][0],self.Q[leftDaughter][n+1][0],self.A[leftDaughter][n+1][0] = self.connections[id]([self.P_pre[id],self.P_pre[leftDaughter]],
                                                                                                                                [self.Q_pre[id],self.Q_pre[leftDaughter]],
                                                                                                                                [self.A_pre[id],self.A_pre[leftDaughter]],
                                                                                                                                 self.dt)
                        else:
                            '''solve bifurcation/connection between mother and left&right daughters'''
                            self.P[id][n+1][-1],self.Q[id][n+1][-1],self.A[id][n+1][-1],\
                            self.P[leftDaughter][n+1][0],self.Q[leftDaughter][n+1][0],self.A[leftDaughter][n+1][0],\
                            self.P[rightDaughter][n+1][0],self.Q[rightDaughter][n+1][0],self.A[rightDaughter][n+1][0] = self.connections[id]([self.P_pre[id],self.P_pre[leftDaughter],self.P_pre[rightDaughter]],
                                                                                                                                   [self.Q_pre[id],self.Q_pre[leftDaughter],self.Q_pre[rightDaughter]],
                                                                                                                                   [self.A_pre[id],self.A_pre[leftDaughter],self.A_pre[rightDaughter]],
                                                                                                                                    self.dt)
        return self.P,self.Q,self.A
    
    
    
    def AFunctionSys0(self,vesselID,A,p):
        return A
    def AFunctionSys1(self,vesselID,A,p):
        return self.vascularNetwork.vessels[vesselID].A(p)
     
    def MacCormack_2_tree(self):
        '''
        MacCormack_2 solver method with forward-euler time steping,
        Using either Characteristic system 0 or 1 as defined in the XML-file.
        
        This method is solving the system by traversing the binary tree of the network
        '''
        print "Solving system \n"
        root = self.vascularNetwork.root[0]
        
        #for n in range(3): 
        for n in range(self.Tsteps-1): #:
            
            '''Predictor Step'''
            vcalc = []
            
            '''Solve system of the root vessel'''
            self.systemEquations[root].updateSystem(self.P[root][n],self.Q[root][n],self.A[root][n])
            
            dPdz  = self.df(self.P[root][n],self.vascularNetwork.vessels[root].dz)
            dQdz  = self.df(self.Q[root][n],self.vascularNetwork.vessels[root].dz) 
            dQ2dz = self.df(pow(self.Q[root][n],2.)/self.A[root][n],self.vascularNetwork.vessels[root].dz)
                        
            self.P_pre[root] = self.P[root][n] - (self.systemEquations[root].M[0][1]*dQdz)*self.dt
            self.Q_pre[root] = self.Q[root][n] - (self.systemEquations[root].M[1][0]*dPdz + self.systemEquations[root].M[1][1]*dQ2dz - self.systemEquations[root].B[1] )*self.dt
            self.A_pre[root] = self.AFunction(root,self.A[root][n],self.P_pre[root])
            
            '''Solve start - boundary Condition of the vessel ''' 
            self.P_pre[root][0],self.Q_pre[root][0],self.A_pre[root][0],self.d0_pre[root] =\
                 self.boundarys[root][0](self.P[root][n],
                                         self.Q[root][n],
                                         self.A[root][n],
                                         self.d0[root],
                                         self.vascularNetwork.vessels[root].z,
                                         n,
                                         self.dt,
                                         self.systemEquations[root].L,
                                         self.systemEquations[root].R,
                                         self.systemEquations[root].LAMBDA)
            
            if self.vascularNetwork.vessels[root].leftDaughter != None:
                vcalc.append(root)
            else:
                '''Solve end - boundary Condition of the vessel if existence'''
                self.P_pre[root][-1],self.Q_pre[root][-1],self.A_pre[root][-1],self.d0_pre[root] =\
                 self.boundarys[root][-1](self.P[root][n],
                                         self.Q[root][n],
                                         self.A[root][n],
                                         self.d0[root],
                                         self.vascularNetwork.vessels[root].z,
                                         n,
                                         self.dt,
                                         self.systemEquations[root].L,
                                         self.systemEquations[root].R,
                                         self.systemEquations[root].LAMBDA)
            
            '''loop through tree until the end of it is reached'''
            while len(vcalc) != 0:
                motherVessel = vcalc.pop(0)
                leftDaughter = self.vascularNetwork.vessels[motherVessel].leftDaughter
                '''Solve system of the left daughter vessel'''
                self.systemEquations[leftDaughter].updateSystem(self.P[leftDaughter][n],self.Q[leftDaughter][n],self.A[leftDaughter][n])
            
                dPdz  = self.df(self.P[leftDaughter][n],self.vascularNetwork.vessels[leftDaughter].dz)
                dQdz  = self.df(self.Q[leftDaughter][n],self.vascularNetwork.vessels[leftDaughter].dz) 
                dQ2dz = self.df(pow(self.Q[leftDaughter][n],2.)/self.A[leftDaughter][n],self.vascularNetwork.vessels[leftDaughter].dz)
                            
                self.P_pre[leftDaughter] = self.P[leftDaughter][n] - (self.systemEquations[leftDaughter].M[0][1]*dQdz)*self.dt
                self.Q_pre[leftDaughter] = self.Q[leftDaughter][n] - (self.systemEquations[leftDaughter].M[1][0]*dPdz + self.systemEquations[leftDaughter].M[1][1]*dQ2dz - self.systemEquations[leftDaughter].B[1] )*self.dt
                self.A_pre[leftDaughter] = self.AFunction(leftDaughter,self.A[leftDaughter][n],self.P_pre[leftDaughter])
            

                if self.vascularNetwork.vessels[leftDaughter].leftDaughter != None:
                    vcalc.append(leftDaughter)
                else:
                    '''Solve end - boundary Condition of the vessel'''
                    self.P_pre[leftDaughter][-1],self.Q_pre[leftDaughter][-1],self.A_pre[leftDaughter][-1],self.d0_pre[leftDaughter] =\
                         self.boundarys[leftDaughter][0](   self.P[leftDaughter][n],
                                                             self.Q[leftDaughter][n],
                                                             self.A[leftDaughter][n],
                                                             self.d0[leftDaughter],
                                                             self.vascularNetwork.vessels[leftDaughter].z,
                                                             n,
                                                             self.dt,
                                                             self.systemEquations[leftDaughter].L,
                                                             self.systemEquations[leftDaughter].R,
                                                             self.systemEquations[leftDaughter].LAMBDA)
                    
                # find right daughter
                if self.vascularNetwork.vessels[motherVessel].rightDaughter != None:
                    rightDaughter = self.vascularNetwork.vessels[motherVessel].rightDaughter
                    
                    ''' Solve system of the right daughter '''
                    self.systemEquations[rightDaughter].updateSystem(self.P[rightDaughter][n],self.Q[rightDaughter][n],self.A[rightDaughter][n])
            
                    dPdz  = self.df(self.P[rightDaughter][n],self.vascularNetwork.vessels[rightDaughter].dz)
                    dQdz  = self.df(self.Q[rightDaughter][n],self.vascularNetwork.vessels[rightDaughter].dz) 
                    dQ2dz = self.df(pow(self.Q[rightDaughter][n],2.)/self.A[rightDaughter][n],self.vascularNetwork.vessels[rightDaughter].dz)
                                
                    self.P_pre[rightDaughter] = self.P[rightDaughter][n] - (self.systemEquations[rightDaughter].M[0][1]*dQdz)*self.dt
                    self.Q_pre[rightDaughter] = self.Q[rightDaughter][n] - (self.systemEquations[rightDaughter].M[1][0]*dPdz + self.systemEquations[rightDaughter].M[1][1]*dQ2dz - self.systemEquations[rightDaughter].B[1] )*self.dt
                    self.A_pre[rightDaughter] = self.AFunction(rightDaughter,self.A[rightDaughter][n],self.P_pre[rightDaughter])
                    
                    
                    if self.vascularNetwork.vessels[rightDaughter].leftDaughter != None:
                        vcalc.append(rightDaughter)
                        
                    else:
                        '''Solve end - boundary Condition of the vessel'''
                        self.P_pre[rightDaughter][-1],self.Q_pre[rightDaughter][-1],self.A_pre[rightDaughter][-1],self.d0_pre[rightDaughter] =\
                         self.boundarys[rightDaughter][0](   self.P[rightDaughter][n],
                                                             self.Q[rightDaughter][n],
                                                             self.A[rightDaughter][n],
                                                             self.d0[rightDaughter],
                                                             self.vascularNetwork.vessels[rightDaughter].z,
                                                             n,
                                                             self.dt,
                                                             self.systemEquations[rightDaughter].L,
                                                             self.systemEquations[rightDaughter].R,
                                                             self.systemEquations[rightDaughter].LAMBDA)
                        
                    '''solve bifurcation between mother and left&right daughters'''
                    
                    self.P_pre[motherVessel][-1],self.Q_pre[motherVessel][-1],self.A_pre[motherVessel][-1],\
                    self.P_pre[leftDaughter][0],self.Q_pre[leftDaughter][0],self.A_pre[leftDaughter][0],\
                    self.P_pre[rightDaughter][0],self.Q_pre[rightDaughter][0],self.A_pre[rightDaughter][0] = self.connections[motherVessel]([self.P[motherVessel][n],self.P[leftDaughter][n],self.P[rightDaughter][n]],
                                                                                                                                             [self.Q[motherVessel][n],self.Q[leftDaughter][n],self.Q[rightDaughter][n]],
                                                                                                                                             [self.A[motherVessel][n],self.A[leftDaughter][n],self.A[rightDaughter][n]],
                                                                                                                                             self.dt)
                        
                else:
                    
                    ''' solve Connection between mother and left daughter'''
                    self.P_pre[motherVessel][-1],self.Q_pre[motherVessel][-1],self.A_pre[motherVessel][-1],\
                    self.P_pre[leftDaughter][0],self.Q_pre[leftDaughter][0],self.A_pre[leftDaughter][0] = self.connections[motherVessel]([self.P[motherVessel][n],self.P[leftDaughter][n]],
                                                                                                                                         [self.Q[motherVessel][n],self.Q[leftDaughter][n]],
                                                                                                                                        [self.A[motherVessel][n],self.A[leftDaughter][n]],
                                                                                                                                             self.dt)
            ########################
                                                                                                                                        
            '''Corrector Step'''
            vcalc = []
            #find root
                        
            self.systemEquations[root].updateSystem(self.P_pre[root],self.Q_pre[root],self.A_pre[root])
            
            dPdz  = self.db(self.P_pre[root],self.vascularNetwork.vessels[root].dz)
            dQdz  = self.db(self.Q_pre[root],self.vascularNetwork.vessels[root].dz)     
            dQ2dz = self.db(pow(self.Q_pre[root],2.)/self.A_pre[root],self.vascularNetwork.vessels[root].dz)
            #print 'Solve system of the root vessel',root
            self.P[root][n+1] = (self.P[root][n] + self.P_pre[root] - (self.systemEquations[root].M[0][1]*dQdz)*self.dt )/2.0
            self.Q[root][n+1] = (self.Q[root][n] + self.Q_pre[root] - (self.systemEquations[root].M[1][0]*dPdz + self.systemEquations[root].M[1][1]*dQ2dz - self.systemEquations[root].B[1] )*self.dt )/2.0
            self.A[root][n+1] = self.AFunction(root,self.A_pre[root],self.P[root][n+1])
            
            #print 'Solve start - boundary Condition of the vessel ',root 
            self.P[root][n+1][0],self.Q[root][n+1][0],self.A[root][n+1][0],self.d0[root] =\
                 self.boundarys[root][0](self.P_pre[root],
                                         self.Q_pre[root],
                                         self.A_pre[root],
                                         self.d0_pre[root],
                                         self.vascularNetwork.vessels[root].z,
                                         n,
                                         self.dt,
                                         self.systemEquations[root].L,
                                         self.systemEquations[root].R,
                                         self.systemEquations[root].LAMBDA)
            
            # add vcalc to the viz vessels if root has daughters:
            if self.vascularNetwork.vessels[root].leftDaughter != None:
                vcalc.append(root)
            else:
                #print 'Solve end - boundary Condition of the vessel if existence',root
                self.P[root][n+1][-1],self.Q[root][n+1][-1],self.A[root][n+1][-1],self.d0[root] =\
                 self.boundarys[root][-1](self.P_pre[root],
                                         self.Q_pre[root],
                                         self.A_pre[root],
                                         self.d0_pre[root],
                                         self.vascularNetwork.vessels[root].z,
                                         n,
                                         self.dt,
                                         self.systemEquations[root].L,
                                         self.systemEquations[root].R,
                                         self.systemEquations[root].LAMBDA)
            '''loop through tree until the end of it is reached'''
            while len(vcalc) != 0:
                motherVessel = vcalc.pop(0)
                leftDaughter = self.vascularNetwork.vessels[motherVessel].leftDaughter
                '''Solve system of the left daughter vessel'''
                self.systemEquations[leftDaughter].updateSystem(self.P_pre[leftDaughter],self.Q_pre[leftDaughter],self.A_pre[leftDaughter])
            
                dPdz  = self.db(self.P_pre[leftDaughter],self.vascularNetwork.vessels[leftDaughter].dz)
                dQdz  = self.db(self.Q_pre[leftDaughter],self.vascularNetwork.vessels[leftDaughter].dz)     
                dQ2dz = self.db(pow(self.Q_pre[leftDaughter],2.)/self.A_pre[leftDaughter],self.vascularNetwork.vessels[leftDaughter].dz)
                #print 'Solve system of the root vessel',root
                self.P[leftDaughter][n+1] = (self.P[leftDaughter][n] + self.P_pre[leftDaughter] - (self.systemEquations[leftDaughter].M[0][1]*dQdz)*self.dt )/2.0
                self.Q[leftDaughter][n+1] = (self.Q[leftDaughter][n] + self.Q_pre[leftDaughter] - (self.systemEquations[leftDaughter].M[1][0]*dPdz + self.systemEquations[leftDaughter].M[1][1]*dQ2dz - self.systemEquations[leftDaughter].B[1] )*self.dt )/2.0
                self.A[leftDaughter][n+1] = self.AFunction(leftDaughter,self.A_pre[leftDaughter],self.P[leftDaughter][n+1])
            

                if self.vascularNetwork.vessels[leftDaughter].leftDaughter != None:
                    vcalc.append(leftDaughter)
                else:
                    '''Solve end - boundary Condition of the vessel'''
                    self.P[leftDaughter][n+1][-1],self.Q[leftDaughter][n+1][-1],self.A[leftDaughter][n+1][-1],self.d0[leftDaughter] =\
                                    self.boundarys[leftDaughter][0](self.P_pre[leftDaughter],
                                                         self.Q_pre[leftDaughter],
                                                         self.A_pre[leftDaughter],
                                                         self.d0_pre[leftDaughter],
                                                         self.vascularNetwork.vessels[leftDaughter].z,
                                                         n,
                                                         self.dt,
                                                         self.systemEquations[leftDaughter].L,
                                                         self.systemEquations[leftDaughter].R,
                                                         self.systemEquations[leftDaughter].LAMBDA)
                    
                # find right daughter
                if self.vascularNetwork.vessels[motherVessel].rightDaughter != None:
                    rightDaughter = self.vascularNetwork.vessels[motherVessel].rightDaughter
                    
                    ''' Solve system of the right daughter '''
                    self.systemEquations[rightDaughter].updateSystem(self.P_pre[rightDaughter],self.Q_pre[rightDaughter],self.A_pre[rightDaughter])
            
                    dPdz  = self.db(self.P_pre[rightDaughter],self.vascularNetwork.vessels[rightDaughter].dz)
                    dQdz  = self.db(self.Q_pre[rightDaughter],self.vascularNetwork.vessels[rightDaughter].dz)     
                    dQ2dz = self.db(pow(self.Q_pre[rightDaughter],2.)/self.A_pre[rightDaughter],self.vascularNetwork.vessels[rightDaughter].dz)
                    #print 'Solve system of the root vessel',root
                    self.P[rightDaughter][n+1] = (self.P[rightDaughter][n] + self.P_pre[rightDaughter] - (self.systemEquations[rightDaughter].M[0][1]*dQdz)*self.dt )/2.0
                    self.Q[rightDaughter][n+1] = (self.Q[rightDaughter][n] + self.Q_pre[rightDaughter] - (self.systemEquations[rightDaughter].M[1][0]*dPdz + self.systemEquations[rightDaughter].M[1][1]*dQ2dz - self.systemEquations[rightDaughter].B[1] )*self.dt )/2.0
                    self.A[rightDaughter][n+1] = self.AFunction(rightDaughter,self.A_pre[rightDaughter],self.P[rightDaughter][n+1])
                
                    if self.vascularNetwork.vessels[rightDaughter].leftDaughter != None:
                        vcalc.append(rightDaughter)
                        
                    else:
                        '''Solve end - boundary Condition of the vessel'''
                        self.P[rightDaughter][n+1][-1],self.Q[rightDaughter][n+1][-1],self.A[rightDaughter][n+1][-1],self.d0[rightDaughter] =\
                                    self.boundarys[rightDaughter][0](self.P_pre[rightDaughter],
                                                         self.Q_pre[rightDaughter],
                                                         self.A_pre[rightDaughter],
                                                         self.d0_pre[rightDaughter],
                                                         self.vascularNetwork.vessels[rightDaughter].z,
                                                         n,
                                                         self.dt,
                                                         self.systemEquations[rightDaughter].L,
                                                         self.systemEquations[rightDaughter].R,
                                                         self.systemEquations[rightDaughter].LAMBDA)
                        
                    '''solve bifurcation between mother and left&right daughters'''
                    
                    self.P[motherVessel][n+1][-1],self.Q[motherVessel][n+1][-1],self.A[motherVessel][n+1][-1],\
                    self.P[leftDaughter][n+1][0],self.Q[leftDaughter][n+1][0],self.A[leftDaughter][n+1][0],\
                    self.P[rightDaughter][n+1][0],self.Q[rightDaughter][n+1][0],self.A[rightDaughter][n+1][0] =\
                     self.connections[motherVessel]([self.P_pre[motherVessel],self.P_pre[leftDaughter],self.P_pre[rightDaughter]],
                                                     [self.Q_pre[motherVessel],self.Q_pre[leftDaughter],self.Q_pre[rightDaughter]],
                                                     [self.A_pre[motherVessel],self.A_pre[leftDaughter],self.A_pre[rightDaughter]],
                                                     self.dt)

                else:
                    
                    ''' solve Connection between mother and left daughter'''
                    self.P[motherVessel][n+1][-1],self.Q[motherVessel][n+1][-1],self.A[motherVessel][n+1][-1],\
                    self.P[leftDaughter][n+1][0],self.Q[leftDaughter][n+1][0],self.A[leftDaughter][n+1][0] =\
                     self.connections[motherVessel]([self.P_pre[motherVessel],self.P_pre[leftDaughter]],
                                                     [self.Q_pre[motherVessel],self.Q_pre[leftDaughter]],
                                                     [self.A_pre[motherVessel],self.A_pre[leftDaughter]],
                                                     self.dt)
        return self.P,self.Q,self.A
        
    def MacCormack_2_treeOpt(self):
        '''
        MacCormack_2 solver method with forward-euler time steping,
        Using either Characteristic system 0 or 1 as defined in the XML-file.
        
        This method is solving the system by traversing the binary tree of the network
        '''
        print "Solving system \n"
        root = self.vascularNetwork.root[0]
        #for n in range(3): 
        for n in range(self.Tsteps-1): #:
            '''Solve system of the root vessel'''
            self.systemEquations[root].updateSystem(self.P[root][n],self.Q[root][n],self.A[root][n])
            
            dPdz  = self.df(self.P[root][n],self.vascularNetwork.vessels[root].dz)
            dQdz  = self.df(self.Q[root][n],self.vascularNetwork.vessels[root].dz) 
            dQ2dz = self.df(pow(self.Q[root][n],2.)/self.A[root][n],self.vascularNetwork.vessels[root].dz)
                        
            self.P_pre[root] = self.P[root][n] - (self.systemEquations[root].M[0][1]*dQdz)*self.dt
            self.Q_pre[root] = self.Q[root][n] - (self.systemEquations[root].M[1][0]*dPdz + self.systemEquations[root].M[1][1]*dQ2dz - self.systemEquations[root].B[1] )*self.dt
            self.A_pre[root] = self.AFunction(root,self.A[root][n],self.P_pre[root])
            
            '''Solve start - boundary Condition of the vessel ''' 
            self.P_pre[root][0],self.Q_pre[root][0],self.A_pre[root][0],self.d0_pre[root] =\
                 self.boundarys[root][0](self.P[root][n],
                                         self.Q[root][n],
                                         self.A[root][n],
                                         self.d0[root],
                                         self.vascularNetwork.vessels[root].z,
                                         n,
                                         self.dt,
                                         self.systemEquations[root].L,
                                         self.systemEquations[root].R,
                                         self.systemEquations[root].LAMBDA)
            
            if self.vascularNetwork.vessels[root].leftDaughter != None:
                vcalc.append(root)
            else:
                '''Solve end - boundary Condition of the vessel if existence'''
                self.P_pre[root][-1],self.Q_pre[root][-1],self.A_pre[root][-1],self.d0_pre[root] =\
                 self.boundarys[root][-1](self.P[root][n],
                                         self.Q[root][n],
                                         self.A[root][n],
                                         self.d0[root],
                                         self.vascularNetwork.vessels[root].z,
                                         n,
                                         self.dt,
                                         self.systemEquations[root].L,
                                         self.systemEquations[root].R,
                                         self.systemEquations[root].LAMBDA)
        
        return self.P,self.Q,self.A
        