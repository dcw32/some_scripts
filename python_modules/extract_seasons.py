import sys
#
#
def seas_counter(seas):
	if seas=='DJF':
		mmm=[0,1,-1]
	elif seas=='MAM':
		mmm=[2,3,4]
	elif seas=='JJA':
		mmm=[5,6,7]
	elif seas=='SON':
		mmm=[8,9,10]
	else:
		print "SEASON NOT DEFINED in ~/python_modules/"
		sys.exit()
	return mmm
def extract_season2d(times,field,field_s,seas):
	nyrs=len(times)/12
	yr=0
	for i in range(12):
		while yr<nyrs:
			field_s[i,:,:]=field[i+12*yr,:,:]+field_s[i,:,:]
			yr=yr+1
		yr=0
	field_s=field_s/nyrs
	mmm=seas_counter(seas)
	print mmm
	field_s=(field_s[mmm[0],:,:]+field_s[mmm[1],:,:]+field_s[mmm[2],:,:])/3.	
	return field_s
