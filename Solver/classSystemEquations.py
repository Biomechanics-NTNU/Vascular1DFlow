import numpy as np 


class System(object):
    
    def __init__(self,vessel,charSys):
        
        # compliance properties
        self.C  = vessel.C
        self.c  = vessel.c
        
        self.C_nID = vessel.C_nID
        self.c_nID = vessel.c_nID
        
        #Fluid properties
        self.my     = vessel.my
        self.rho    = vessel.rho
        self.dlt    = vessel.dlt
        self.gamma  =  vessel.gamma
                
        # system matrices
        self.M      = None
        self.B      = None
        self.LAMBDA = [None,None]  # list with 2 objects one for each boundary
        self.R      = [None,None]  # list with 2 objects one for each boundary
        self.L      = [None,None]  # list with 2 objects one for each boundary
        
        # eigen values
        if charSys == 0:
            self.updateLARL = self.updateLARLSys0
        elif charSys == 1:
            self.updateLARL = self.updateLARLSys1
    
    def updateSystem(self,P,Q,A):
        
        # calculate needed values
        C = self.C(P) 
        c = self.c(A,P)
        v = Q/A
        # set up M matrix
        m12 = 1./C
        m21 = C*(c**2-self.dlt*v**2)
        m22 = 2*self.dlt*v
        self.M = [[0.0,m12],[m21,m22]]
        # set up B matrix
        b2 = -(1.0/self.rho)*2*np.pi*v*(self.gamma+2)*self.my 
        self.B = [0.0,b2]
        # set up LAMBDA, L and R matrices
        self.updateLARL(P,Q,A,Ct=C,ct=c)
        
    def updateLARLSys0(self,P,Q,A,idArray=[0,-1],bool='all',Ct=None,ct=None):
        
        for id in idArray:       
            
            if Ct == None: C = self.C_nID(P,id)
            else:  C = Ct[id]
            if ct == None: c = self.c_nID(A,P,id)
            else:  c = ct[id]
                        
            #l1 = c
            #l2 = -c
            Z1 =  1.0/(C*c)
            Z2 = -1.0/(C*-c)
            
            if bool != 'L': 
                self.LAMBDA[id] = np.array([c,-c])
            
           
            #if bool != 'L': self.R[id] = np.array([[Z1,-Z2],[1.0,1.0]])
            #self.L[id] = (1.0/(Z1+Z2))*np.array([[1.0,Z2],[-1.0,Z1]])
            
            #print 1./(C*c), Z1
            #if bool != 'L': self.R[id] = np.array([[1. ,1.  ],[C*c,-C*c]])
            #self.L[id] = np.array([[1., 1./(C*c)],[1.,-1./(C*c)]])/2.

                self.R[id] = np.array([[1. ,1.  ],[1.0/Z1,-1.0/Z2]])
            #self.L[id] = (Z1*Z2)/(Z1+Z2)*np.array([[1/Z2 , 1.],[1./Z1,-1.]])
            self.L[id] = np.array([[Z1/(Z1+Z2),(Z1*Z2)/(Z1+Z2)],[Z2/(Z1+Z2),-(Z1*Z2)/(Z1+Z2)]])
            

    def updateLARLSys1(self,P,Q,A,idArray=[0,-1],bool='all'):
        
        for id in idArray:       
            
            C = self.C_nID(P,id)
            c = self.c_nID(A,P,id)
            
            v = Q[id]/A[id] 
            
            l1 = self.dlt*v+np.sqrt(c**2.+self.dlt*(self.dlt-1.)*(v)**2.0)
            l2 = self.dlt*v-np.sqrt(c**2.+self.dlt*(self.dlt-1.)*(v)**2.0)
            
            Z1 =  1.0/(C*l1)
            Z2 = -1.0/(C*l2)
            
            if bool != 'L': 
                self.LAMBDA[id] = np.array([l1,l2])
            
            
            
            #if bool != 'L': self.R[id] = np.array([[Z1,-Z2],[1.0,1.0]])
            #self.L[id] = (1.0/(Z1+Z2))*np.array([[1.0,Z2],[-1.0,Z1]])
            
            #print 1./(C*c), Z1
            #if bool != 'L': self.R[id] = np.array([[1. ,1.  ],[C*c,-C*c]])
            #self.L[id] = np.array([[1., 1./(C*c)],[1.,-1./(C*c)]])/2.
            
            
                self.R[id] = np.array([[1. ,1.  ],[1.0/Z1,-1.0/Z2]])
            #self.L[id] = (Z1*Z2)/(Z1+Z2)*np.array([[1/Z2 , 1.],[1./Z1,-1.]])
            self.L[id] = np.array([[Z1/(Z1+Z2),(Z1*Z2)/(Z1+Z2)],[Z2/(Z1+Z2),-(Z1*Z2)/(Z1+Z2)]])



    
    
        
    