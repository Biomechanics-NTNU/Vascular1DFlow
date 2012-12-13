import numpy as np 

class Compliance:		
	def c(self,A,P):						
		return np.sqrt(A/(self.rho*self.C(P)))
	def c0(self,A,P):						
		return np.sqrt(A/(self.rho*self.C0(P)))
	def c_(self,A,P,ID):			
		return np.sqrt(A[ID]/(self.rho*self.C_(P,ID)))
	def c0_(self,A,P,ID):
		return np.sqrt(A[ID]/(self.rho*self.C0_(P,ID)))

class Exponential(Compliance):
	def __init__(self,solid,rho):
		self.rho = rho					
		self.As = solid['As']
		self.beta = solid['beta']
		self.Ps = solid['Ps']
		self.P0 = solid['P0']

	def P00(self,A):				
		return self.Ps*np.exp(self.beta*(A/self.As-1.))
	def A(self,P):
		return self.As*(np.log((P+self.P00)/self.Ps)/self.beta + 1.)
	def C(self,P):
		return self.As/(self.beta*(P+self.P0))	
	def C0(self,P):
		return self.As/(self.beta*(self.Ps+self.P0))
	def A_(self,P,ID):
		return self.As[ID]*(np.log((P[ID]+self.P0[ID])/self.Ps)/self.beta[ID] + 1.)
	def C_(self,P,ID):
		return self.As[ID]/(self.beta[ID]*(P[ID]+self.P0[ID]))	
	def C0_(self,P,ID):
		return self.As[ID]/(self.beta[ID]*(self.Ps+self.P0[ID]))


class Laplace(Compliance):
	def __init__(self,solid,rho):							
		self.As,self.H,self.E,self.Ps = solid['As'],solid['H'],solid['E'],solid['Ps']
		self.beta = (np.sqrt(np.pi)*self.H*self.E)/self.As
		self.rho = rho	
	def A(self,P):
		return (P/self.beta + np.sqrt(self.As))**2.
	def C(self,P):
		return (2.*(P/self.beta + np.sqrt(self.As)))/self.beta
	def C0(self,P):
		return (2.*(self.Ps/self.beta + np.sqrt(self.As)))/self.beta
	def A_(self,P,ID):
		return (P[ID]/self.beta[ID] + np.sqrt(self.As[ID]))**2.
	def C_(self,P,ID):
		return (2.*(P[ID]/self.beta[ID] + np.sqrt(self.As[ID])))/self.beta[ID]
	def C0_(self,P,ID):
		return (2.*(self.Ps/self.beta[ID] + np.sqrt(self.As[ID])))/self.beta[ID]

	

class Laplace2(Compliance):
	def __init__(self,solid,rho):							
		self.As,self.beta,self.Ps = solid['As'],solid['beta'],solid['Ps']
		self.rho = rho	
	def A(self,P):
		return (P/self.beta + np.sqrt(self.As))**2.
	def C(self,P):
		return (2.*(P/self.beta + np.sqrt(self.As)))/self.beta
	def C0(self,P):
		return (2.*(self.Ps/self.beta + np.sqrt(self.As)))/self.beta
	def A_(self,P,ID):
		return (P[ID]/self.beta[ID] + np.sqrt(self.As[ID]))**2.
	def C_(self,P,ID):
		return (2.*(P[ID]/self.beta[ID] + np.sqrt(self.As[ID])))/self.beta[ID]
	def C0_(self,P,ID):
		return (2.*(self.Ps/self.beta[ID] + np.sqrt(self.As[ID])))/self.beta[ID]

	