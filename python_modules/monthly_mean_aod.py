import numpy as N
def mmean(arrayin):
	a=N.shape(arrayin)
	yrs=a[1]/12
	arrayout=N.empty([6,yrs,12,a[2],a[3]])
        stdout=N.empty([6,yrs,12,a[2],a[3]])
	for i in range(0,a[1]):
		mo=i % 12
		yr=(i-mo)/12
		arrayout[:,yr,mo,:,:]=arrayin[:,i,:,:]
	stdout=N.std(arrayout,axis=1)
	arrayout=N.mean(arrayout,axis=1)
	return (arrayout, stdout)
