import numpy as np
from numpy.linalg import solve
from scipy.optimize import fsolve

class Connection(object):
    
    def __init__(self,mothers,mothersSys,daughters,daugthersSys,equSys):
        # vessel variables initially set, constant through simulation
        self.rho = []
        self.systemEquations = []
        self.z = []
        self.c_func = []
        self.A_func = []
        
        self.positions =[]
        self.vz = []
        
        self.names = []
        self.number = 0
        
        # calculated in call, constant in fsolve
        self.zi    = []
        self.P     = []
        self.Q     = []
        self.PQ    = None
        self.A     = []
        self.domega = []
        # changed in fsolve
        self.du    = []
            
        # equations to solve in f solve
        self.fsolveFunction = None
        
        #initialize
        for mother,motherSys in zip(mothers,mothersSys):
            self.rho.append(mother.rho)
            self.z.append(mother.z)
            self.systemEquations.append(motherSys)
            self.positions.append(-1)
            self.names.append(mother.Id)
            self.c_func.append(mother.c_nID)
            self.A_func.append(mother.A_nID)
            self.vz.append(-1)
        
        for daughter,daughterSys in zip(daughters,daugthersSys):
            self.rho.append(daughter.rho)
            self.z.append(daughter.z)
            self.systemEquations.append(daughterSys)
            self.positions.append(0)
            self.names.append(daughter.Id)
            self.c_func.append(daughter.c_nID)
            self.A_func.append(daughter.A_nID)
            self.vz.append(1)
            
        if len(self.rho) == len(self.z) == len(self.systemEquations) == len(self.positions) == len(self.names):
            self.number = len(self.rho)
        else:
            raise Error("classConnections: Indifferent length of connection variables")
        
        if self.number == 2:
            
            if   equSys == '02':  self.fsolveFunction = self.fsolveConnectionSys0
            elif equSys == '2':  self.fsolveFunction = self.fsolveConnectionSys1
            else: raise Error("classConnections: EquSys not properly defined!")
                     
        elif self.number == 3:
            
            if   equSys == '02':  self.fsolveFunction = self.fsolveConnectionSysMulti
            elif equSys == '2':  self.fsolveFunction = self.fsolveBifurcationSys1
            else: raise Error("classConnections: EquSys not properly defined!")
        else:
            raise Error("classConnections: Too many vessels given")
        
        
        a = range(self.number)
        self.P  = [None for i in a]
        self.PQ = [None for i in a]
        self.Q  = [None for i in a]
        self.A  = [None for i in a]
        self.domega= [None for i in a]
        # changed in fsolve
        self.du= [None for i in a]
            
    
    def __call__(self,P,Q,A,dt):
        '''
        P,Q,A must be list in the form [P1 (mother) ,P2 (LD) ..(,P3 (RD))  ]
        in the same form as the bifucrtaion/connection was initialized
        '''
        if len(A) != len(Q) != len(P) != self.number:
            raise Error(" Wrong length of input variables")
        # solution vector
        x = []
                
        for i in range(self.number):
            
            pos = self.positions[i]
            self.P[i] = P[i][pos]
            self.Q[i] = Q[i][pos]
            self.A[i] = A[i][pos]
            
            zi = self.z[i][pos] + self.vz[i]  * self.c_func[i](A[i],P[i],pos) * dt
            
            self.du[i] = np.array([np.interp(zi,self.z[i],P[i]),np.interp(zi,self.z[i],Q[i])])-np.array([self.P[i],self.Q[i]])  
                        
            self.systemEquations[i].updateLARL(P[i],Q[i],A[i],idArray=[pos],bool='L') #
            
            self.domega[i] =  np.dot(self.systemEquations[i].L[pos][pos+1],self.du[i])

            x.extend([self.P[i],self.Q[i]])
               
        sol = fsolve(self.fsolveFunction,x)
        sol = sol.tolist()
        for i in range(self.number):
            sol.insert(i*3+2,self.A[i]) ## change to pSolution
        
        return sol
    
    def fsolveConnectionSysMulti(self,x):
        
        #self.count = self.count+1
        
        out = [0 for i in x]
        #outt = ['' for i in x]
        for i in range(self.number):
           
            Ai = self.A[i]
            out[0]=out[0] + self.vz[i]*x[i*2+1]
            if i != 0: out[i]= (x[0]+(self.rho[0]/2.)*(x[1]/self.A[0])**2.) - (x[i*2]+(self.rho[i]/2.)*(x[i*2+1]/Ai)**2.)
                                 
            du_t = np.array([x[i*2],x[i*2+1]])-np.array([self.P[i],self.Q[i]])  
            self.systemEquations[i].updateLARL([x[i*2]],[x[i*2+1]],[Ai],idArray=[self.positions[i]],bool='L') 
            
            out[i+self.number]= np.dot(self.systemEquations[i].L[self.positions[i]][self.positions[i]+1],du_t)/self.domega[i] - 1
            
            #outt[0] = outt[0]+str('mass eq'+str(i)+'+')
            #if i != 0: outt[i]= 'P0 - P'+str(i)
            #outt[i+self.number] = 'domega '+str(i)
            
        #if self.count == 8:print outt
        return out
    
    def fsolveConnectionSys0(self,x):
        
        P1,Q1,P2,Q2 = x
        #P1,Q1,P2,Q2,P3,Q3
        A1 = self.A[0]
        A2 = self.A[1]
      
        res1 = self.vz[0]*Q1+self.vz[1]*Q2
        
        res2 = -P1-(self.rho[0]/2.)*(Q1/A1)**2.+P2+(self.rho[1]/2.)*(Q2/A2)**2.
        
        du1 = np.array([P1,Q1])-np.array([self.P[0],self.Q[0]])  
        du2 = np.array([P2,Q2])-np.array([self.P[1],self.Q[1]]) 
        
        self.systemEquations[0].updateLARL([P1],[Q1],[A1],idArray=[self.positions[0]],bool='L') 
        self.systemEquations[1].updateLARL([P2],[Q2],[A2],idArray=[self.positions[1]],bool='L') 
        
        res3 = np.dot(self.systemEquations[0].L[self.positions[0]][self.positions[0]+1],du1)-self.domega[0]
        res4 = np.dot(self.systemEquations[1].L[self.positions[1]][self.positions[1]+1],du2)-self.domega[1]
        
        return [res1,res2,res3,res4]
    
    def fsolveConnectionSys1(self,x):
        
        P1,Q1,P2,Q2 = x
        A1 = self.A_func[self.positions[0]]([P1],self.positions[0])
        A2 = self.A_func[self.positions[1]]([P2],self.positions[1])
                
        out = [ self.vz[0]*Q1+self.vz[1]*Q2 ]
        out.append(P1+(self.rho[0]/2.)*(Q1/A1)**2.-P2-(self.rho[1]/2.)*(Q2/A2)**2.)
        
        # into for loop
        du1 = np.array([P1,Q1])-np.array([self.P[0],self.Q[0]])  
        du2 = np.array([P2,Q2])-np.array([self.P[1],self.Q[1]])  
                
        self.systemEquations[0].updateLARL([P1],[Q1],[A1],idArray=[self.positions[0]])
        self.systemEquations[1].updateLARL([P2],[Q2],[A2],idArray=[self.positions[1]])
        
        out.append(np.dot(self.systemEquations[0].L[self.positions[0]+1],du1)-self.omega[0])
        out.append(np.dot(self.systemEquations[1].L[self.positions[0]+1],du2)-self.omega[1])
        
        return out
    
    
    def fsolveBifurcation(self):
        pass 