import numpy as N
def aeroflux(flux):
	flux=N.sum(flux,(1,2,3))
	flux=flux[12]
	flux=N.mean(flux)
	return flux
