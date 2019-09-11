import numpy as N
def mmean(arrayin):
	a=N.shape(arrayin)
	print a
	yrs=a[0]/12
	arrayout=N.empty([yrs,12,a[1],a[2],a[3]])
        stdout=N.empty([yrs,12,a[1],a[2],a[3]])
	for i in range(0,a[0]):
		mo=i % 12
		yr=(i-mo)/12
		arrayout[yr,mo,:,:,:]=arrayin[i,:,:,:]
	stdout=N.std(arrayout,axis=0)
	arrayout=N.mean(arrayout,axis=0)
	return (arrayout, stdout)
