import numpy as np 
from numpy.linalg import solve
from scipy.optimize import fsolve
from numpy.linalg import inv

import sys,os
# set the path relative to THIS file not the executing file!
cur = os.path.dirname( os.path.realpath( __file__ ) )
sys.path.append(cur+'/NetworkLib')
from classBoundaryConditions import *

warn = False
maxfevfsolve = 10
warningActivated = False

class Boundary_02:
    def __init__(self,boundaryConditions,A_nID,equSys):
        '''
        Constructor of Boundary_02
        
        Initializes one Boundary of a vessel
        input: - boundaryConditions:
                  a list of all boundaryConditions of the vessel at one position i.g 0 or -1 

        Note: if network consits of a single vessel make sure the function is called twice 
              once for the start and once for the end boundary at 0 and -1!
              
        variables (most important):
        self.position        =    position of the boundary
        self.duFunction        =    function which gives du-vector back
        self.bcType1        =    list of all type1 functions found in the given boundaryConditions
        self.omegaInput        =    function which gives the _omega of omega_ back depending on self.position
        self.omegaFunctio     =    function with boundaryCondtion of type2 calculating the omega-vector
        
        __call__:
        returns p,q, omega
        '''
        self.position = None

        #Function which gives du-vector back
        self.duFunction = None
        # list of all type1 functions found in the given boundaryConditions, which determine 
        # the du-Function
        self.bcType1 = []
        #Function which gives the _omega of omega_ back depending on self.position
        self.omegaInput = None
        #Function with boundaryCondtion of type2 calculating the omega-vector
        self.omegaFunction = None
        #list of all type2 functions found in the given boundaryConditions
        self.bcType2 = []
        
        posTemp = []
        for bC in boundaryConditions:            
            if bC.type == 1:
                self.bcType1.append(bC)
            elif bC.type == 2:
                self.bcType2.append(bC)
            posTemp.append(bC.position)
            
        #find out position of the Boundary + check if all conditions are defined at the same side
        if sum(posTemp) == 0: self.position = 0
        elif sum(np.array(posTemp)+1) == 0: self.position = -1
        else: raise ValueError("One position of one boundaryCondition is not correct")
        
        # 1. Set duFunction 
        if len(self.bcType1) == 0:    
            # set calculation type for Standard Boundary
            self.duFunction = self.duFunctionZero
        elif len(self.bcType1) == 1:
            self.bcType1 = self.bcType1[0]
            self.duFunction = self.duFunctionSingle
        elif len(self.bcType1) > 0:
            self.duFunction = self.duFunctionMulti
        else:
            raise Error("Too many type2-boundary Conditions defined!") 
        
        # 2. omegaInput Function
        if self.position == 0:
            self.omegaInput = self.omegaInputPos0
        elif self.position == -1:
            self.omegaInput = self.omegaInputPos1
        
        # 3. Set Condition there should only one! if None apply the one condition depending on the side
        if len(self.bcType2) == 1:
            self.omegaFunction = self.bcType2[0]
            
        elif len(self.bcType2) == 0:    
            # set calculation type for Standard Boundary
            self.omegaFunction = standardType2()
            self.omegaFunction.setPosition(self.position)
        else:
            raise Error("Too many type2-boundary Conditions defined!")
        
        # 4. Define the output of A, dependend on the characteristic system 0.1
        if equSys == '02':
            self.AFunction = self.AFunctionSys0
        elif equSys == '2': 
            self.A_nID = A_nID
            self.AFunction = self.AFunctionSys1
        else:
            raise Error("EquSys not properly defined!")
    ### Function which calculated du
    def duFunctionZero(self,n,dt):
        '''
        Determine the du-vector = [0,0] if no type1 boundaryConditions are given
        '''
        return np.array([0,0])
    
    def duFunctionSingle(self,n,dt):
        '''
        Determine the du-vector with the values of the given
        type1 boundaryCondition 
        '''
        return self.bcType1(n+1,dt)-self.bcType1(n,dt)
    
    def duFunctionMulti(self,n,dt):
        '''
        Determine the summized du-vector with the values of all given
        type1 boundaryConditions 
        '''
        du = np.array([0,0])
        for bc in self.bcType1:
            du = du + (bc(n+1,dt)-bc(n,dt))
        return du
    
    ### Functions for omegaInput determined by the position of the boundary
    def omegaInputPos0(self,lmbd,L,dt,Q,P,z,position):
        '''
        calculate the outgoing charakteristic _omega (w2) at the boundary
        and return it
        '''
        zi = z[position]-lmbd[1]*dt   
        _du = np.array([np.interp(zi,z,P),np.interp(zi,z,Q)])-np.array([P[position],Q[position]])
        
        _domega = np.dot(L[1],_du)
        return _domega
    
    def omegaInputPos1(self,lmbd,L,dt,Q,P,z,position):
        '''
        calculate the ingoing charakteristic _omega (w1) at the boundary
        and return it
        '''
        zi = z[position]-lmbd[0]*dt
        du_ = np.array([np.interp(zi,z,P),np.interp(zi,z,Q)])-np.array([P[position],Q[position]])
        
        domega_ = np.dot(L[0],du_)
        return domega_
    
    ## Function to define the output of A, dependend on the characteristic system 0.1
    
    def AFunctionSys0(self,A,p):
        return A
    
    def AFunctionSys1(self,A,p):
        return self.A_nID(p,self.position)
    
    def __call__(self,P,Q,A,d0,z,n,dt,L,R,LMBD):    
        '''
        P
        Q
        A
        d0
        z
        n
        dt
        systemEquation.L
        systemEquation.R
        systemEquation.LMBD
        '''
        position = self.position
        #calculate the du_vector using given boundaryConditions of type 2
        du_temp = self.duFunction(n,dt)
        
        #'get omega_ or _omega depending on the position'
        _omega_temp = self.omegaInput(LMBD[position],L[position],dt,Q,P,z,position)
        
        #calculate the omgea using given boundaryConditions of type 2
        omega = self.omegaFunction(_omega_temp,d0[position],du_temp,R[position],L[position],n,dt)
        
        p,q = np.dot(R[position],omega) + np.array([P[position],Q[position]])
        
        a = self.AFunction(A[position],[p])
        
        d0[position] = omega
        
        return p,q,a,d0
    
