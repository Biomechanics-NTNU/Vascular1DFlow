import numpy as np 
from scipy.optimize import fsolve
from math import exp

class TimeStep:
	def __init__(self,dz,c,CFL):
		self.dz,self.c,self.CFL = \
		dz,c,CFL
	def __call__(self):
		return (self.CFL*self.dz)/self.c
	
class BoundaryCondition(object):
	
	position = None
	name = None
	
	def update(self,bcDict):
		'''
        updates the updateBoundaryDict data using a dictionary in from of 
	    bcDict = {'variableName': value}
        '''
		for key,value in bcDict.iteritems():
			try:
				self.__getattribute__(key)
				self.__setattr__(key,value)
			except: 
				print 'ValueError: wrong key: %s, could not set up boundaryCondition' %key
				
	def setPosition(self, position):
		'''
		Set the position of the boundaryCondition
		'''
		self.position = position
		
########################################################################################
# Type 1 boundary Conditions
########################################################################################

class Mean(BoundaryCondition):
	"""	
	Boundary profile - type 1
	
	ramps to a mean amplitude
	
	call function input (n,dt)
	gives amplitude value for pressure/flow back
	for the function delta g
	"""
	def __init__(self):		
		
		self.type = 1
			
		self.mean = 0
		self.Tmean = 1
		self.Traise = 1
		
		self.duMatrix = np.array([0,0])
	
	def __call__(self,n,dt):
		t = n*dt
		
		if t <= self.Tmean:
			mean = 0.0
		else:
			if t <= self.Tmean + self.Traise:
				mean = self.mean*pow(np.sin(np.pi*(t - self.Tmean)/(self.Traise*2.0)),2.0)
			else:
				mean = self.mean
				
		return mean*self.duMatrix
		
class Sinus(BoundaryCondition):
	"""	
	Boundary profile - type 1
	
	creates in a periodic sinus signal
	
	call function input (n,dt)
	gives amplitude value for pressure/flow back
	for the function delta g
	"""
	
	def __init__(self):		
		
		self.type = 1
		
		self.amp 	= 0
		self.Npulse = 1
		self.freq 	= 1
		
		self.duMatrix = np.array([0,0])
		
	def __call__(self,n,dt):
		t = n*dt
		if t <= self.Npulse/self.freq:
			amp = self.amp*np.sin(np.pi*self.freq*t)
		else:
			amp = 0.0
		
		return amp*self.duMatrix

class Sinus2(BoundaryCondition):
	"""	
	Boundary profile - type 1
	
	creates in a periodic sinus-squared signal
	
	call function input (n,dt)
	gives amplitude value for pressure/flow back
	for the function delta g
	"""
	def __init__(self):
		
		self.type = 1
		
		self.amp 	= 0
		self.Npulse = 1
		self.freq 	= 1
		self.Tpulse = 1
		
		self.duMatrix = np.array([0,0])
		
	def __call__(self,n,dt):
		t = n*dt
		if t <= self.Tpulse:
			ampT = 0.0
		else:
			if t <= self.Tpulse + self.Npulse/self.freq:
				ampT = self.amp*pow(np.sin(np.pi*self.freq*(t - self.Tpulse)),2.0)
			else:
				ampT = 0.0
				
		return ampT*self.duMatrix
    
class PhysiologicalFunction(BoundaryCondition):
	"""	
	Boundary profile - type 1
	
	creates in a similar heart-outflow signal as found in
	Stergiopulos et al. 1992
	
	call function input (n,dt)
	gives amplitude value for flow back
	for the function delta g
	"""
	def __init__(self):
		
		self.type = 1
		
		# physiological pulse
		self.amp    = 0
		self.Npulse = 1
		self.freq   = 1
		self.Tpulse = 1
		self.Tspace = 1

		self.lowPoint = 0.1739

		self.fracSin2 = 0.36
		self.fracCos = 0.43
		self.fracRes = 1.0 - (self.fracSin2+self.fracCos)
		self.pulseTime = []
	
	def update(self,bcDict):
		'''
        updates the updateBoundaryDict data using a dictionary in from of 
	    bcDict = {'variableName': value}
        '''
		for key,value in bcDict.iteritems():
			try:
				self.__getattribute__(key)
				self.__setattr__(key,value)
			except: 
				print 'ValueError: wrong key: %s, could not set up boundaryCondition' %key	
				
		for puls in np.arange(self.Npulse):
			print "i am in here"
			self.pulseTime.append([self.Tpulse+ (self.Tspace+1.0/self.freq)*puls, self.Tpulse+1.0 / self.freq+(self.Tspace+1.0/self.freq)*puls])
    
	def __call__(self,n,dt):
		t = n*dt
		ampT = 0.0
		for pTA in self.pulseTime:
			if pTA[0] < t and t < pTA[1]:
				pulsNum = self.pulseTime.index(pTA)
				if t < self.Tpulse+pulsNum*(self.Tspace+1.0/self.freq)+self.fracSin2*1/self.freq:
					ampT = self.amp*pow(np.sin(np.pi*self.freq/(2*self.fracSin2)*(t - (self.Tpulse+pulsNum*(self.Tspace+1.0/self.freq)))),2.0)
                
                		elif t < self.Tpulse+pulsNum*(self.Tspace+1.0/self.freq)+(self.fracSin2+self.fracCos)*1/self.freq:
                			ampT = self.amp*(np.cos((np.pi*self.freq/(2.0*self.fracCos)*(t-self.fracSin2/self.freq- (self.Tpulse+pulsNum*(self.Tspace+1.0/self.freq))))) )
                	
        			elif t < self.Tpulse+pulsNum*(self.Tspace+1.0/self.freq)+(self.fracSin2+self.fracCos+self.fracRes*3.0/8.0)*1/self.freq:
                    			ampT =self.amp*(pow(np.sin((t-(self.fracSin2+self.fracCos)/self.freq-(self.Tpulse+pulsNum*(self.Tspace+1.0/self.freq))-(1.0/(self.freq)*self.fracRes)*3./8.)*np.pi*self.freq/(1.5*self.fracRes)),2.0)*2*self.lowPoint-self.lowPoint)
        			else:
                			ampT =self.amp*(pow(np.sin((t-(self.fracSin2+self.fracCos)/self.freq-(self.Tpulse+pulsNum*(self.Tspace+1.0/self.freq))-(1.0/(self.freq)*self.fracRes)*3./8.0)*np.pi*self.freq/(self.fracRes*1.25)),2.0)*self.lowPoint-self.lowPoint)
		
		return ampT*np.array([0,1])
    
class expVelocity(BoundaryCondition):
	"""	
	Boundary profile - type 1
	
	creates a single gaussian peak signal as found
	in D.Xiu and S.P.Sherwin 2007
	
	call function input (n,dt)
	gives amplitude value for flow back
	for the function delta g
	"""
	def __init__(self):
		
		self.type = 1
		
		self.amp 	= 0
		self.C 		= 1
		self.Tpulse = 0
		
		self.area = 0

	def __call__(self,n,dt):
		t = n*dt
		
		ampT =  self.amp * exp(-self.C*(t-self.Tpulse)**2)
		
		return ampT*np.array([0,1])*self.area[self.position]

class Fourier(BoundaryCondition):
	"""	
	Boundary profile - type 1
	
	creates fourier signal similar to the flow of the heart
	
	call function input (n,dt)
	gives amplitude value for flow back
	for the function delta g
	"""
	def __init__(self):
		
		self.type = 1
				
		self.scale  = 0
		self.Npulse = 1
		self.Tpulse = 1
		
		self.harm = [[0.2465,0.1975,0.2624,0.2132,0.0424,0.0722,0.0411,0.0093,0.0141,0.0044],
				 [0.0,2.2861,4.5723,6.8584,9.1445,11.4307,13.7168,16.0029,18.2891,20.5752],
				 [0.0,2.5010,-2.9986,-0.2689,1.4904,-2.9854,-0.0317,1.5292,-2.5394,0.5771]]
		
	def __call__(self,n,dt):
		t = n*dt
		if t <= self.Tpulse:
			amp = 0.0
		else:
			if t <= self.Tpulse + self.Npulse/self.harm[1][1]:
				amp = 0.0
				for i in range(len(self.harm[0])):
					amp += self.harm[0][i]*np.sin(2.0*np.pi*self.harm[1][i]*(t - self.Tpulse) + self.harm[2][i] + np.pi/2.0)
			else:
				amp = 0.0
			amp *= self.scale

		return amp*np.array([0,1])

class Physiological(BoundaryCondition):
	"""	
	Boundary profile - type 1
	
	creates signal based on measured physiological values
	(source unknown: values in NetworkLib/physiologicalData.py)
	
	call function input (n,dt)
	gives amplitude value for pressure/flow back
	for the function delta g
	"""
	def __init__(self):
		
		self.type = 1
		
		from physiologicalData import aorticFlowPressure
		self.data = aorticFlowPressure()
		self.tmax = max(self.data.t)
	
	def __call__(self,n,dt):
		t = n*dt
		pulseNumber = np.floor(t/self.tmax)
		P = np.interp(t-pulseNumber*self.tmax,self.data.t,self.data.P)
		Q = np.interp(t-pulseNumber*self.tmax,self.data.t,self.data.Q)
		return Q*np.array([0,1])
		#return P*np.array([1,0])
		
	
		
########################################################################################
# Type 2 boundary Conditions
########################################################################################
class standardType2(BoundaryCondition):
	"""	
	Boundary profile - type 2
	
	calculates domega if only type1 function(s) are given
	
	call function input:
	 _domega_,dO,du,R,L,n,dt
	returns the domega-vector with (domega_ , _domega) based on the input values
	and its returnFunction
	"""
	def __init__(self):
		self.type = 2
		self.name = "standardType2"

		self.returnFunction = None

	def __call__(self,_domega_,dO,du,R,L,n,dt):
		if self.returnFunction != None: return self.returnFunction(_domega_,dO,du,R,L,n,dt)
		else: raise Error("Retrun function of standardType2 is not proper defined")
	
	def setPosition(self, position):
		'''
		Set the position of the boundaryCondition
		and determines the return function
		based on the position of the bC (0,-1)
		'''
		self.position = position
		if self.position == 0:
			self.returnFunction = self.funcPos0
		elif self.position ==-1:
			self.returnFunction = self.funcPos1
		else:
			self.returnFunction = None
	
	def funcPos0(self,_domega_,dO,du,R,L,n,dt):
		'''
		return function for position 0 at the start
		of the vessel
		'''
		domega_prescribed = np.dot(L[0],du)
		
		return np.array([domega_prescribed , _domega_])
	
	def funcPos1(self,_domega_,dO,du,R,L,n,dt):
		'''
		return function for position -1 at the end
		of the vessel
		'''
		prescribed_domega = np.dot(L[1],du*np.array([1,-1]))
		
		return np.array([_domega_ , prescribed_domega])


class Terminal(BoundaryCondition):
	"""	
	Boundary profile - type 2
	
	Terminal reflection
	
	call function input:
	 _domega_,dO,du,R,L,n,dt
	returns the domega-vector with (domega_ , _domega) based on the input values
	and its returnFunction
	
	"""
	def __init__(self):
		self.type = 2
		self.Rt = 0
	
		self.returnFunction = None
	
	def __call__(self,_domega_,dO,du,R,L,n,dt):
		if self.returnFunction != None: return self.returnFunction(_domega_,dO,du,R,L,n,dt)
		else: raise Error("Retrun function of TerminalReflection is not proper defined")
	
	def setPosition(self, position):
		'''
		Set the position of the boundaryCondition
		and determines the return function
		based on the position of the bC (0,-1)
		'''
		self.position = position
		if self.position == 0:
			self.returnFunction = self.funcPos0
		elif self.position ==-1:
			self.returnFunction = self.funcPos1
		else:
			self.returnFunction = None
	
	def funcPos0(self,_domega_,dO,du,R,L,n,dt):
		'''
		return function for position 0 at the start
		of the vessel
		'''
		domega_prescribed = np.dot(L[0],du)
		domega_ = -_domega_*self.Rt + domega_prescribed
		return np.array([domega_ , _domega_])
	
	def funcPos1(self,_domega_,dO,du,R,L,n,dt):
		'''
		return function for position -1 at the end
		of the vessel
		'''
		prescribed_omgega = np.dot(L[1],du*np.array([1,-1]))
		
		_domega = -_domega_*self.Rt + prescribed_omgega
		return np.array([_domega_ , _domega])
	

class Resistance(BoundaryCondition):
	"""	
	Boundary profile - type 2
	
	signle Resistance element
	
	call function input:
	 _domega_,dO,du,R,L,n,dt
	returns the domega-vector with (domega_ , _domega) based on the input values
	and its returnFunction
	"""
	def __init__(self):
		self.type = 2
		self.Rf = 1
		
		self.returnFunction = None
		
	def __call__(self,_domega_,dO,du,R,L,n,dt):
		if self.returnFunction != None: return self.returnFunction(_domega_,dO,du,R,L,n,dt)
		else: raise Error("Retrun function of Resistance is not proper defined")
	
	def setPosition(self, position):
		'''
		Set the position of the boundaryCondition
		and determines the return function
		based on the position of the bC (0,-1)
		'''
		self.position = position
		if self.position == 0:
			self.returnFunction = self.funcPos0
		elif self.position ==-1:
			self.returnFunction = self.funcPos1
		else:
			self.returnFunction = None
	
	def funcPos0(self,_domega_,dO,du,R,L,n,dt):
		'''
		return function for position 0 at the start
		of the vessel
		'''
		P0 = du[0]/2. ##??
		
		r21,r22 = R[1][0],R[1][1]
		if self.Rf == None: self.Rf = -1./R[1][1]	
		denom = 1.+r21*self.Rf
		domega_ = ((-1.-r22*Rf)*_domega_ + P0)/denom
		return np.array([domega_ , _domega_])
	
	def funcPos1(self,_domega_,dO,du,R,L,n,dt):
		'''
		return function for position -1 at the end
		of the vessel
		'''
		P0 = du[0]/2. ##??
		if self.Rf == None: self.Rf = 1./R[1][0]	
		denom = 1.-r22*self.Rf
		_domega = ((-1.+r21*self.Rf)*_domega_ + P0)/denom
			
		return np.array([_domega_ , _domega])
	

class Windkessel2(BoundaryCondition):
	"""	
	Boundary profile - type 2
	
	2 Element Windkessel 
	
	call function input:
	 _domega_,dO,du,R,L,n,dt
	returns the domega-vector with (domega_ , _domega) based on the input values
	and its returnFunction
	"""
	def __init__(self):
		self.type = 2
		
		self.Rc = 1
		self.C  = 0

		self.returnFunction = None

	def __call__(self,_domega_,dO,du,R,L,n,dt):
		if self.returnFunction != None: return self.returnFunction(_domega_,dO,du,R,L,n,dt)
		else: raise Error("Retrun function of 2-Element Windkessel is not proper defined")
	
	def setPosition(self, position):
		'''
		Set the position of the boundaryCondition
		and determines the return function
		based on the position of the bC (0,-1)
		'''
		self.position = position
		if self.position == 0:
			self.returnFunction = self.funcPos0
		elif self.position ==-1:
			self.returnFunction = self.funcPos1
		else:
			self.returnFunction = None
	
	def funcPos0(self,_domega_,dO,du,R,L,n,dt):
		'''
		return function for position 0 at the start
		of the vessel
		'''
		P0 = du[0]/2.
		r21,r22 = R[1][0],R[1][1]
		dw_,_dw = dO[self.position][0],dO[self.position][1]
		
		if self.Rc == None: self.Rc = -1./R[1][1]	
		tau = 2.*self.Rc*self.C
		taudt = tau/dt
		denom = 1.+taudt+r21*self.Rc
		domega_ = ((-1.-taudt-r22*self.Rc)*_domega_ + taudt*(dw_+_dw) + P0)/denom
		return np.array([domega_ , _domega_])
	
	def funcPos1(self,_domega_,dO,du,R,L,n,dt):
		'''
		return function for position -1 at the end
		of the vessel
		'''
		P0 = du[0]/2.
		r21,r22 = R[1][0],R[1][1]
		dw_,_dw = dO[self.position][0],dO[self.position][1]
		
		if self.Rc == None: self.Rc = 1./R[1][0]
		tau = 2.*self.Rc*self.C
		taudt = tau/dt
		denom = 1.+taudt-r22*self.Rc
		_domega = ((-1.-taudt+r21*self.Rc)*_domega_ + taudt*(dw_+_dw) + P0)/denom
		return np.array([_domega_ , _domega])



class Windkessel3(BoundaryCondition):
	"""	
	Boundary profile - type 2
	
	3 Element Windkessel 
	
	call function input:
	 _domega_,dO,du,R,L,n,dt
	returns the domega-vector with (domega_ , _domega) based on the input values
	and its returnFunction
	"""
	def __init__(self):
		self.type = 2
		
		self.Rc = 1
		self.C  = 0
		self.R1 = None
		self.RT = None
		
		self.returnFunction = None
		
	def __call__(self,_domega_,dO,du,R,L,n,dt):
		if self.returnFunction != None: return self.returnFunction(_domega_,dO,du,R,L,n,dt)
		else: raise Error("Retrun function of Resistance is not proper defined")
	
	def setPosition(self, position):
		'''
		Set the position of the boundaryCondition
		and determines the return function
		based on the position of the bC (0,-1)
		'''
		self.position = position
		if self.position == 0:
			self.returnFunction = self.funcPos0
		elif self.position ==-1:
			self.returnFunction = self.funcPos1
		else:
			self.returnFunction = None
	
	def funcPos0(self,_domega_,dO,du,R,L,n,dt):
		'''
		return function for position 0 at the start
		of the vessel
		'''
		P0 = du[0]/2.

		r21,r22 = R[1][0],R[1][1]
		dw_,_dw = dO[self.position][0],dO[self.position][1]
		
		if self.R1 == None: self.R1 = -1./r22
		if self.Rc == None: self.Rc = self.RT - self.R1 
		tau = 2.*self.Rc*self.C
		taudt = tau/dt
		denom = taudt + 1. + r21*(self.R1*taudt + self.Rc + self.R1)
		domega_ = ((1.+r21*self.R1)*taudt*dw_ + (1.+r22*self.R1)*taudt*_dw + (-1.-taudt-r22*(self.R1*taudt + self.Rc + self.R1))*_domega_ + P0)/denom
		return np.array([domega_ , _domega_])
	
	def funcPos1(self,_domega_,dO,du,R,L,n,dt):
		'''
		return function for position -1 at the end
		of the vessel
		'''
		P0 = du[0]/2.

		r21,r22 = R[1][0],R[1][1]
		dw_,_dw = dO[self.position][0],dO[self.position][1]

		if self.R1 == None: self.R1 = 1./r21
		if self.Rc == None: self.Rc = self.RT - self.R1 
		tau = 2*self.Rc*self.C
		taudt = tau/dt	
		denom =  taudt + 1. - r22*(self.R1*taudt + self.Rc + self.R1)
		_domega = ((1.-r21*self.R1)*taudt*dw_ + (1.-r22*self.R1)*taudt*_dw + (-1.-taudt+r21*(self.R1*taudt + self.Rc + self.R1))*_domega_ + P0)/denom
		return np.array([_domega_ , _domega])

		print "3 element windkessel is not correct defined"

class L_network(BoundaryCondition):
	"""	
	Boundary profile - type 2
	
	L-network
	
	call function input:
	 _domega_,dO,du,R,L,n,dt
	returns the domega-vector with (domega_ , _domega) based on the input values
	and its returnFunction
	"""
	def __init__(self):
		self.type = 2
		
		self.C  = 0
		self.R1 = 1

		self.returnFunction = None

	def __call__(self,_domega_,dO,du,R,L,n,dt):
		if self.returnFunction != None: return self.returnFunction(_domega_,dO,du,R,L,n,dt)
		else: raise Error("Retrun function of Resistance is not proper defined")
	
	def setPosition(self, position):
		'''
		Set the position of the boundaryCondition
		and determines the return function
		based on the position of the bC (0,-1)
		'''
		self.position = position
		if self.position == 0:
			self.returnFunction = self.funcPos0
		elif self.position ==-1:
			self.returnFunction = self.funcPos1
		else:
			self.returnFunction = None
	
	def funcPos0(self,_domega_,dO,du,R,L,n,dt):
		'''
		return function for position 0 at the start
		of the vessel
		'''
		dQ0 = du[1]/2.0

		r21,r22 = R[1][0],R[1][1]
		dw_,_dw = dO[self.position][0],dO[self.position][1]
		if self.position == 0:
			if self.R1 == None: self.R1 = -1./R[1][1]	
			tau = 2.*self.R1*self.C
			taudt = tau/dt	
			denom = taudt + taudt*r21*self.R1 + r21*self.R1
			domega_ = ((1.+r21*self.R1)*taudt*dw_ + (1.+r22*self.R1)*taudt*_dw + (-taudt - taudt*r22*self.R1 - r22*self.R1)*_domega_ + self.R1*dQ0)/denom
			return np.array([domega_ , _domega_])
	
	def funcPos1(self,_domega_,dO,du,R,L,n,dt):
		'''
		return function for position -1 at the end
		of the vessel
		'''
		dQ0 = du[1]/2.0
		
		if self.R1 == None: self.R1 = 1./R[1][0]	
		tau = 2.*self.R1*self.C
		taudt = tau/dt
		denom = -taudt + taudt*r22*self.R1 + r22*self.R1
		_domega = ((-1.+r21*self.R1)*taudt*dw_ + (-1.+r22*self.R1)*taudt*_dw + (taudt - taudt*r21*self.R1 - r21*self.R1)*_domega_ + self.R1*dQ0)/denom
		return np.array([_domega_ , _domega])
	





