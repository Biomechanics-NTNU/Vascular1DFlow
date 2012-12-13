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

class Connection:
	def __init__(self,P1,Q1,A1,P2,Q2,A2):
		self.P1,self.Q1,self.A1,self.P2,self.Q2,self.A2 = P1,Q1,A1,P2,Q2,A2 

def connect_02(x,*args):
	pq1,pq2,A01,A02,domega1,domega2,rho1,rho2,S1,S2,I1,I2 = args
	P1,Q1,P2,Q2 = x
	
	du1 = np.array([P1,Q1])-pq1
	du2 = np.array([P2,Q2])-pq2
	out = [I1[2]*Q1+I2[2]*Q2]
	out.append(P1+(rho1/2.)*(Q1/A01)**2.-P2-(rho2/2.)*(Q2/A02)**2.)
	out.append(np.dot(S1.L(np.array([P1,Q1,A01]),I1[0])[I1[1]],du1)-domega1)
	out.append(np.dot(S2.L(np.array([P2,Q2,A02]),I2[0])[I2[1]],du2)-domega2)
	return out
 
class Connection_02(Connection):
	def __call__(self,n,dt,Vessel1,Vessel2,ID):			
		P1,Q1,A1,P2,Q2,A2 = self.P1,self.Q1,self.A1,self.P2,self.Q2,self.A2		
	
 		I1,I2 = ID[1:4],ID[5:8]
		z1,dz1,rho1,S1 = Vessel1['S'].Segment
		z2,dz2,rho2,S2 = Vessel2['S'].Segment
		zi1 = z1[I1[0]] + I1[2]*S1.c_(rho1,A1[I1[0]],P1[I1[0]],I1[0])*dt
		zi2 = z2[I2[0]] + I2[2]*S2.c_(rho2,A2[I2[0]],P2[I2[0]],I2[0])*dt
		PQ1 = np.array([P1[I1[0]],Q1[I1[0]]])
		PQ2 = np.array([P2[I2[0]],Q2[I2[0]]])
		du1 = np.array([np.interp(zi1,z1,P1),np.interp(zi1,z1,Q1)])-PQ1             	
		du2 = np.array([np.interp(zi2,z2,P2),np.interp(zi2,z2,Q2)])-PQ2
		domega1 = np.dot(S1.L(np.array([P1[I1[0]],Q1[I1[0]],A1[I1[0]]]),I1[0])[I1[1]],du1)
		domega2 = np.dot(S2.L(np.array([P2[I2[0]],Q2[I2[0]],A2[I2[0]]]),I2[0])[I2[1]],du2)

		x = [P1[I1[0]],Q1[I1[0]],P2[I2[0]],Q2[I2[0]]]
		if warningActivated == True:
			P1,Q1,P2,Q2 = fsolve(connect_02,x,warning=warn,maxfev=maxfevfsolve,args=(PQ1,PQ2,A1[I1[0]],A2[I2[0]],domega1,domega2,rho1,rho2,S1,S2,I1,I2))
		else:
			P1,Q1,P2,Q2 = fsolve(connect_02,x,args=(PQ1,PQ2,A1[I1[0]],A2[I2[0]],domega1,domega2,rho1,rho2,S1,S2,I1,I2))
		return P1,Q1,A1[I1[0]],P2,Q2,A2[I2[0]]
	
	

def connect_2(x,*args):
	pq1,pq2,domega1,domega2,rho1,rho2,S1,S2,I1,I2 = args
	P1,Q1,A1,P2,Q2,A2 = x

	du1 = np.array([P1,Q1])-pq1
	du2 = np.array([P2,Q2])-pq2
	out = [I1[2]*Q1+I2[2]*Q2]
	out.append(P1+(rho1/2.)*(Q1/A1)**2.-P2-(rho2/2.)*(Q2/A2)**2.)
	out.append(S1.A_(P1,I1[0])-A1)
	out.append(S2.A_(P2,I2[0])-A2)
	out.append(np.dot(S1.L(np.array([P1,Q1,A1]),I1[0])[I1[1]],du1)-domega1)
	out.append(np.dot(S2.L(np.array([P2,Q2,A2]),I2[0])[I2[1]],du2)-domega2)
	return out

class Connection_2(Connection):
	def __call__(self,n,dt,Vessel1,Vessel2,ID):			
		P1,Q1,A1,P2,Q2,A2 = self.P1,self.Q1,self.A1,self.P2,self.Q2,self.A2	
	
 		I1,I2 = ID[1:4],ID[5:8]
		z1,dz1,rho1,S1 = Vessel1['S'].Segment
		z2,dz2,rho2,S2 = Vessel2['S'].Segment
		zi1 = z1[I1[0]] + I1[2]*S1.c_(rho1,A1[I1[0]],P1[I1[0]],I1[0])*dt
		zi2 = z2[I2[0]] + I2[2]*S2.c_(rho2,A2[I2[0]],P2[I2[0]],I2[0])*dt
		PQ1 = np.array([P1[I1[0]],Q1[I1[0]]])
		PQ2 = np.array([P2[I2[0]],Q2[I2[0]]])
		du1 = np.array([np.interp(zi1,z1,P1),np.interp(zi1,z1,Q1)])-PQ1             	
		du2 = np.array([np.interp(zi2,z2,P2),np.interp(zi2,z2,Q2)])-PQ2
		domega1 = np.dot(S1.L(np.array([P1[I1[0]],Q1[I1[0]],A1[I1[0]]]),I1[0])[I1[1]],du1)
		domega2 = np.dot(S2.L(np.array([P2[I2[0]],Q2[I2[0]],A2[I2[0]]]),I2[0])[I2[1]],du2) 

		x = [P1[I1[0]],Q1[I1[0]],A1[I1[0]],P2[I2[0]],Q2[I2[0]],A2[I2[0]]]	
		if warningActivated == True:
			P1,Q1,A1,P2,Q2,A2 = fsolve(connect_2,x,warning=warn,maxfev=maxfevfsolve,args=(PQ1,PQ2,domega1,domega2,rho1,rho2,S1,S2,I1,I2))
		else:
			P1,Q1,A1,P2,Q2,A2 = fsolve(connect_2,x,args=(PQ1,PQ2,domega1,domega2,rho1,rho2,S1,S2,I1,I2))
		return P1,Q1,A1,P2,Q2,A2

def connect_3(x,*args):
	pqa1,pqa2,domega1,domega2,rho1,rho2,S1,S2,I1,I2,dt = args
	P1,Q1,A1,P2,Q2,A2 = x
	PQA1,PQA2 = np.array([P1,Q1,A1]),np.array([P2,Q2,A2])
	du1,du2 = PQA1-pqa1,PQA2-pqa2
	S1L,S2L = S1.L(PQA1,I1[0]),S2.L(PQA2,I2[0])

	out = [I1[2]*Q1+I2[2]*Q2]
	out.append(P1+(rho1/2.)*(Q1/A1)**2.-P2-(rho2/2.)*(Q2/A2)**2.)
	out.append(np.dot(S1L[2],du1))
	out.append(np.dot(S2L[2],du2))
	out.append(np.dot(S1L[I1[1]],du1)-domega1)
	out.append(np.dot(S2L[I2[1]],du2)-domega2)
	return out

class Connection_3(Connection):
	def __call__(self,n,dt,Vessel1,Vessel2,ID):			
		P1,Q1,A1,P2,Q2,A2 = self.P1,self.Q1,self.A1,self.P2,self.Q2,self.A2	
	
 		I1,I2 = ID[1:4],ID[5:8]
		z1,dz1,rho1,S1 = Vessel1['S'].Segment
		z2,dz2,rho2,S2 = Vessel2['S'].Segment
		zi1 = z1[I1[0]] + I1[2]*S1.c_(rho1,A1[I1[0]],P1[I1[0]],I1[0])*dt
		zi2 = z2[I2[0]] + I2[2]*S2.c_(rho2,A2[I2[0]],P2[I2[0]],I2[0])*dt
		PQA1 = np.array([P1[I1[0]],Q1[I1[0]],A1[I1[0]]])
		PQA2 = np.array([P2[I2[0]],Q2[I2[0]],A2[I2[0]]])
		du1 = np.array([np.interp(zi1,z1,P1),np.interp(zi1,z1,Q1),np.interp(zi1,z1,A1)])-PQA1             	
		du2 = np.array([np.interp(zi2,z2,P2),np.interp(zi2,z2,Q2),np.interp(zi2,z2,A2)])-PQA2
		domega1 = np.dot(S1.L(PQA1,I1[0])[I1[1]],du1) 
		domega2 = np.dot(S2.L(PQA2,I2[0])[I2[1]],du2) 

		x = [P1[I1[0]],Q1[I1[0]],A1[I1[0]],P2[I2[0]],Q2[I2[0]],A2[I2[0]]]		
		if warningActivated == True:
			P1,Q1,A1,P2,Q2,A2 = fsolve(connect_3,x,warning=warn,maxfev=maxfevfsolve,args=(PQA1,PQA2,domega1,domega2,rho1,rho2,S1,S2,I1,I2,dt))
		else:
			P1,Q1,A1,P2,Q2,A2 = fsolve(connect_3,x,args=(PQA1,PQA2,domega1,domega2,rho1,rho2,S1,S2,I1,I2,dt))
		return P1,Q1,A1,P2,Q2,A2

class Bifurcation:
	def __init__(self,P1,Q1,A1,P2,Q2,A2,P3,Q3,A3):
		self.P1,self.Q1,self.A1,self.P2,self.Q2,self.A2,self.P3,self.Q3,self.A3 = P1,Q1,A1,P2,Q2,A2,P3,Q3,A3  

def bifurcation_02(x,*args):
	pq1,pq2,pq3,A01,A02,A03,domega1,domega2,domega3,rho1,rho2,rho3,S1,S2,S3,I1,I2,I3 = args
	P1,Q1,P2,Q2,P3,Q3 = x
	du1 = np.array([P1,Q1])-pq1
	du2 = np.array([P2,Q2])-pq2
	du3 = np.array([P3,Q3])-pq3
	out = [I1[2]*Q1+I2[2]*Q2+I3[2]*Q3]
	out.append(P1+(rho1/2.)*(Q1/A01)**2.-P2-(rho2/2.)*(Q2/A02)**2.)
	out.append(P1+(rho1/2.)*(Q1/A01)**2.-P3-(rho3/2.)*(Q3/A03)**2.)
	out.append(np.dot(S1.L(np.array([P1,Q1,A01]),I1[0])[I1[1]],du1)-domega1)
	out.append(np.dot(S2.L(np.array([P2,Q2,A02]),I2[0])[I2[1]],du2)-domega2)
	out.append(np.dot(S3.L(np.array([P3,Q3,A03]),I3[0])[I3[1]],du3)-domega3)
	return out

class Bifurcation_02(Bifurcation):
	def __call__(self,n,dt,Vessel1,Vessel2,Vessel3,ID):			
		P1,Q1,A1,P2,Q2,A2,P3,Q3,A3 = self.P1,self.Q1,self.A1,self.P2,self.Q2,self.A2,self.P3,self.Q3,self.A3		
	
 		I1,I2,I3= ID[1:4],ID[5:8],ID[9:12]
		z1,dz1,rho1,S1 = Vessel1['S'].Segment
		z2,dz2,rho2,S2 = Vessel2['S'].Segment
		z3,dz3,rho3,S3 = Vessel3['S'].Segment
		zi1 = z1[I1[0]] + I1[2]*S1.c_(rho1,A1[I1[0]],P1[I1[0]],I1[0])*dt
		zi2 = z2[I2[0]] + I2[2]*S2.c_(rho2,A2[I2[0]],P2[I2[0]],I2[0])*dt
		zi3 = z3[I3[0]] + I3[2]*S3.c_(rho3,A3[I3[0]],P3[I3[0]],I3[0])*dt
		PQ1 = np.array([P1[I1[0]],Q1[I1[0]]])
		PQ2 = np.array([P2[I2[0]],Q2[I2[0]]])
		PQ3 = np.array([P3[I3[0]],Q3[I3[0]]])
		du1 = np.array([np.interp(zi1,z1,P1),np.interp(zi1,z1,Q1)])-PQ1             	
		du2 = np.array([np.interp(zi2,z2,P2),np.interp(zi2,z2,Q2)])-PQ2
		du3 = np.array([np.interp(zi3,z3,P3),np.interp(zi3,z3,Q3)])-PQ3
		domega1 = np.dot(S1.L(np.array([P1[I1[0]],Q1[I1[0]],A1[I1[0]]]),I1[0])[I1[1]],du1)
		domega2 = np.dot(S2.L(np.array([P2[I2[0]],Q2[I2[0]],A2[I2[0]]]),I2[0])[I2[1]],du2)
		domega3 = np.dot(S3.L(np.array([P3[I3[0]],Q3[I3[0]],A3[I3[0]]]),I3[0])[I3[1]],du3)

		x = [P1[I1[0]],Q1[I1[0]],P2[I2[0]],Q2[I2[0]],P3[I3[0]],Q3[I3[0]]]	
		if warningActivated == True:
			P1,Q1,P2,Q2,P3,Q3 = fsolve(bifurcation_02,x,warning=warn,maxfev=maxfevfsolve,args=(PQ1,PQ2,PQ3,\
								A1[I1[0]],A2[I2[0]],A3[I3[0]],domega1,domega2,domega3,rho1,rho2,rho3,S1,S2,S3,I1,I2,I3))
		else:
			P1,Q1,P2,Q2,P3,Q3 = fsolve(bifurcation_02,x,args=(PQ1,PQ2,PQ3,\
								A1[I1[0]],A2[I2[0]],A3[I3[0]],domega1,domega2,domega3,rho1,rho2,rho3,S1,S2,S3,I1,I2,I3))
		return P1,Q1,A1[I1[0]],P2,Q2,A2[I2[0]],P3,Q3,A3[I3[0]]

#def bifurcation_2(x,*args):
#	pq1,pq2,pq3,domega1,domega2,domega3,rho1,rho2,rho3,S1,S2,S3,I1,I2,I3 = args
#	P1,Q1,A1,P2,Q2,A2,P3,Q3,A3 = x
#	du1 = np.array([P1,Q1])-pq1
#	du2 = np.array([P2,Q2])-pq2
#	du3 = np.array([P3,Q3])-pq3
#	out = [I1[2]*Q1+I2[2]*Q2+I3[2]*Q3]
#	out.append(P1+(rho1/2.)*(Q1/A1)**2.-P2-(rho2/2.)*(Q2/A2)**2.)
#	out.append(P1+(rho1/2.)*(Q1/A1)**2.-P3-(rho3/2.)*(Q3/A3)**2.)
#	out.append(S1.A_(P1,I1[0])-A1)
#	out.append(S2.A_(P2,I2[0])-A2)
#	out.append(S3.A_(P3,I3[0])-A3)
#	out.append(np.dot(S1.L(np.array([P1,Q1,A1]),I1[0])[I1[1]],du1)-domega1)
#	out.append(np.dot(S2.L(np.array([P2,Q2,A2]),I2[0])[I2[1]],du2)-domega2)
#	out.append(np.dot(S3.L(np.array([P3,Q3,A3]),I3[0])[I3[1]],du3)-domega3)
#	return out

def bifurcation_2(x,*args):
	pq1,pq2,pq3,domega1,domega2,domega3,rho1,rho2,rho3,S1,S2,S3,I1,I2,I3 = args
	P1,Q1,P2,Q2,P3,Q3 = x
	du1 = np.array([P1,Q1])-pq1
	du2 = np.array([P2,Q2])-pq2
	du3 = np.array([P3,Q3])-pq3
	out = [I1[2]*Q1+I2[2]*Q2+I3[2]*Q3]
	out.append(P1+(rho1/2.)*(Q1/S1.A_(P1,I1[0]))**2.-P2-(rho2/2.)*(Q2/S2.A_(P2,I2[0]))**2.)
	out.append(P1+(rho1/2.)*(Q1/S1.A_(P1,I1[0]))**2.-P3-(rho3/2.)*(Q3/S3.A_(P3,I3[0]))**2.)
	out.append(np.dot(S1.L(np.array([P1,Q1,S1.A_(P1,I1[0])]),I1[0])[I1[1]],du1)-domega1)
	out.append(np.dot(S2.L(np.array([P2,Q2,S2.A_(P2,I2[0])]),I2[0])[I2[1]],du2)-domega2)
	out.append(np.dot(S3.L(np.array([P3,Q3,S3.A_(P3,I3[0])]),I3[0])[I3[1]],du3)-domega3)
	return out

class Bifurcation_2(Bifurcation):
	def __call__(self,n,dt,Vessel1,Vessel2,Vessel3,ID):			
		P1,Q1,A1,P2,Q2,A2,P3,Q3,A3 = self.P1,self.Q1,self.A1,self.P2,self.Q2,self.A2,self.P3,self.Q3,self.A3		
	
 		I1,I2,I3= ID[1:4],ID[5:8],ID[9:12]
		z1,dz1,rho1,S1 = Vessel1['S'].Segment
		z2,dz2,rho2,S2 = Vessel2['S'].Segment
		z3,dz3,rho3,S3 = Vessel3['S'].Segment
		PQ1 = np.array([P1[I1[0]],Q1[I1[0]]])
		PQ2 = np.array([P2[I2[0]],Q2[I2[0]]])
		PQ3 = np.array([P3[I3[0]],Q3[I3[0]]])
		zi1 = z1[I1[0]] + I1[2]*S1.c_(rho1,A1[I1[0]],P1[I1[0]],I1[0])*dt
		zi2 = z2[I2[0]] + I2[2]*S2.c_(rho2,A2[I2[0]],P2[I2[0]],I2[0])*dt
		zi3 = z3[I3[0]] + I3[2]*S3.c_(rho3,A3[I3[0]],P3[I3[0]],I3[0])*dt
		du1 = np.array([np.interp(zi1,z1,P1),np.interp(zi1,z1,Q1)])-PQ1             	
		du2 = np.array([np.interp(zi2,z2,P2),np.interp(zi2,z2,Q2)])-PQ2
		du3 = np.array([np.interp(zi3,z3,P3),np.interp(zi3,z3,Q3)])-PQ3
		domega1 = np.dot(S1.L(np.array([P1[I1[0]],Q1[I1[0]],A1[I1[0]]]),I1[0])[I1[1]],du1)
		domega2 = np.dot(S2.L(np.array([P2[I2[0]],Q2[I2[0]],A2[I2[0]]]),I2[0])[I2[1]],du2)
		domega3 = np.dot(S3.L(np.array([P3[I3[0]],Q3[I3[0]],A3[I3[0]]]),I3[0])[I3[1]],du3)
#		x = [P1[I1[0]],Q1[I1[0]],A1[I1[0]],P2[I2[0]],Q2[I2[0]],A2[I2[0]],P3[I3[0]],Q3[I3[0]],A3[I3[0]]]
#		P1,Q1,A1,P2,Q2,A2,P3,Q3,A3 = fsolve(bifurcation_2,x,args=(PQ1,PQ2,PQ3,\
#			domega1,domega2,domega3,rho1,rho2,rho3,S1,S2,S3,I1,I2,I3))
		x = [P1[I1[0]],Q1[I1[0]],P2[I2[0]],Q2[I2[0]],P3[I3[0]],Q3[I3[0]]]
		P1,Q1,P2,Q2,P3,Q3 = fsolve(bifurcation_2,x,args=(PQ1,PQ2,PQ3,\
			domega1,domega2,domega3,rho1,rho2,rho3,S1,S2,S3,I1,I2,I3))
		A1,A2,A3 = S1.A_(P1,I1[0]),S2.A_(P2,I2[0]),S3.A_(P3,I3[0])
		return P1,Q1,A1,P2,Q2,A2,P3,Q3,A3

def bifurcation_3(x,*args):
	pqa1,pqa2,pqa3,domega1,domega2,domega3,rho1,rho2,rho3,S1,S2,S3,I1,I2,I3 = args
	P1,Q1,A1,P2,Q2,A2,P3,Q3,A3 = x
	PQA1,PQA2,PQA3 = np.array([P1,Q1,A1]),np.array([P2,Q2,A2]),np.array([P3,Q3,A3])
	du1,du2,du3 = PQA1-pqa1,PQA2-pqa2,PQA3-pqa3
	S1L,S2L,S3L = S1.L(PQA1,I1[0]),S2.L(PQA2,I2[0]),S3.L(PQA3,I3[0])

	out = [I1[2]*Q1+I2[2]*Q2+I3[2]*Q3]
	out.append(P1+(rho1/2.)*(Q1/A1)**2.-P2-(rho2/2.)*(Q2/A2)**2.)
	out.append(P1+(rho1/2.)*(Q1/A1)**2.-P3-(rho3/2.)*(Q3/A3)**2.)
	out.append(np.dot(S1L[2],du1))
	out.append(np.dot(S2L[2],du2))
	out.append(np.dot(S3L[2],du3))
	out.append(np.dot(S1L[I1[1]],du1)-domega1)
	out.append(np.dot(S2L[I2[1]],du2)-domega2)
	out.append(np.dot(S3L[I3[1]],du3)-domega3)
	return out

class Bifurcation_3(Bifurcation):
	def __call__(self,n,dt,Vessel1,Vessel2,Vessel3,ID):			
		P1,Q1,A1,P2,Q2,A2,P3,Q3,A3 = self.P1,self.Q1,self.A1,self.P2,self.Q2,self.A2,self.P3,self.Q3,self.A3		
	
 		I1,I2,I3= ID[1:4],ID[5:8],ID[9:12]
		z1,dz1,rho1,S1 = Vessel1['S'].Segment
		z2,dz2,rho2,S2 = Vessel2['S'].Segment
		z3,dz3,rho3,S3 = Vessel3['S'].Segment	
		PQA1 = np.array([P1[I1[0]],Q1[I1[0]],A1[I1[0]]])
		PQA2 = np.array([P2[I2[0]],Q2[I2[0]],A2[I2[0]]])
		PQA3 = np.array([P3[I3[0]],Q3[I3[0]],A3[I3[0]]])
		zi1 = z1[I1[0]] + I1[2]*S1.c_(rho1,A1[I1[0]],P1[I1[0]],I1[0])*dt
		zi2 = z2[I2[0]] + I2[2]*S2.c_(rho2,A2[I2[0]],P2[I2[0]],I2[0])*dt
		zi3 = z3[I3[0]] + I3[2]*S3.c_(rho3,A3[I3[0]],P3[I3[0]],I3[0])*dt
		du1 = np.array([np.interp(zi1,z1,P1),np.interp(zi1,z1,Q1),np.interp(zi1,z1,A1)])-PQA1 
		du2 = np.array([np.interp(zi2,z2,P2),np.interp(zi2,z2,Q2),np.interp(zi2,z2,A2)])-PQA2
		du3 = np.array([np.interp(zi3,z3,P3),np.interp(zi3,z3,Q3),np.interp(zi3,z3,A3)])-PQA3
		domega1 = np.dot(S1.L(PQA1,I1[0])[I1[1]],du1)
		domega2 = np.dot(S2.L(PQA2,I2[0])[I2[1]],du2)
		domega3 = np.dot(S3.L(PQA3,I3[0])[I3[1]],du3)

		x = [P1[I1[0]],Q1[I1[0]],A1[I1[0]],P2[I2[0]],Q2[I2[0]],A2[I2[0]],P3[I3[0]],Q3[I3[0]],A3[I3[0]]]
		if warningActivated == True:
			P1,Q1,A1,P2,Q2,A2,P3,Q3,A3 = fsolve(bifurcation_3,x,warning=warn,maxfev=maxfevfsolve,args=(PQA1,PQA2,PQA3,\
										domega1,domega2,domega3,rho1,rho2,rho3,S1,S2,S3,I1,I2,I3))
		else:
			P1,Q1,A1,P2,Q2,A2,P3,Q3,A3 = fsolve(bifurcation_3,x,args=(PQA1,PQA2,PQA3,\
										domega1,domega2,domega3,rho1,rho2,rho3,S1,S2,S3,I1,I2,I3))
		return P1,Q1,A1,P2,Q2,A2,P3,Q3,A3


