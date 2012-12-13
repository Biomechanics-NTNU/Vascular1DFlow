import numpy as np 

class Uni:
	""" 1D uniform grid """
	def __init__(self,args):
		self.l = args['l']
		self.R = args['R']
		self.N = args['N']
		self.A = np.pi*self.R**2.
			
	def __call__(self):
		l,A,N = self.l,self.A,self.N
		z  = np.linspace(0.,l,N)
		dz = np.diff(z)
		A0  = np.ones(N)*A
		return [z,dz,A0]

class Cons:
	""" 1D grid with constriction"""
	def __init__(self,args):
		self.l,self.R,self.r,self.N = args['l'],args['R'],args['r'],args['N']
		self.A = np.pi*self.R**2.
		self.a = np.pi*self.r**2.

	def __call__(self):
		l,A,a,N = self.l,self.A,self.a,self.N
		z  = np.linspace(0.,l,N)
		dz = np.diff(z)
		z_cons  = np.linspace(-1.,1.,N)
		D  = np.sqrt((4./np.pi)*A)
		d  = np.sqrt((4./np.pi)*a)
		D0 = (d-D)*(1.0 + np.cos(z_cons*np.pi))/2.
		D0 = D + D0
		A0 = (np.pi/4.)*D0**2.
		return [z,dz,A0]

class Cone:
	""" 1D cone grid """
	def __init__(self,args):
		self.l,self.R,self.r,self.N = args['l'],args['R'],args['r'],args['N']
		self.A = np.pi*self.R**2.
		self.a = np.pi*self.r**2.

	def __call__(self):
		l,A,a,N = self.l,self.A,self.a,self.N
		z  = np.linspace(0.,l,N)
		dz = np.diff(z)
		D  = np.sqrt((4./np.pi)*A)
		d  = np.sqrt((4./np.pi)*a)
		D0 = D + ((d-D)/l*z)
		A0 = (np.pi/4.)*D0**2.
		return [z,dz,A0]

