import numpy as np 


def DeltaForward(f,h): 
	""" 1st-order forward diff array """
	dfdz_t = (f[2::] - f[1:-1])/h[1::]
	
	dfdz = np.empty_like(f)

	dfdz[1:-1] = dfdz_t
	dfdz[0] = 0
	dfdz[-1] = 0
	
	return dfdz

def DeltaBackward(f,h): 
	""" 1st-order backward diff array"""
	dbdz_t = (f[1:-1] - f[0:-2])/h[0:-1]
	
	dbdz = np.empty_like(f)
	dbdz[1:-1] = dbdz_t
	dbdz[0] = 0
	dbdz[-1] = 0
	return dbdz



def DeltaForwardInnerPart(f,h): 
	""" 1st-order forward diff array """
	return (f[2::] - f[1:-1])/h[1::]

def DeltaBackwardInnerPart(f,h): 
	""" 1st-order backward diff array"""
	return (f[1:-1] - f[0:-2])/h[0:-1]
