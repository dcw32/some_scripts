''' Convert the time of the UM output form hours since forecast reference time to 
yyyy-mm-dd-hh-mm-ss
output variable is called "time"
'''
import scipy
import numpy as np
import calendar
import datetime

def convert_time(word,array):	
    start = word
    time = array 
# convert the string to a list of integers
    start2=[]				# year,month,day,hour,minute,second
    start2.append(int(start[0:4]))
    start2.append(int(start[5:7]))
    start2.append(int(start[8:10]))
    start2.append(int(start[11:13]))
    start2.append(int(start[14:16]))
    start2.append(int(start[17:]))

# add the hours passed since the forecast reference time
    start4=[]
    for i in time:
        y=i/(24*360)			# year: 360 days per year
        m=(y-int(y))*12			# month: 12 months per year
        d=(m-int(m))*30			# day: 30 days per month
        mn=(d-int(d))*24		# hour: 24 hours per day
        hr=(mn-int(mn))*60		# minute: 60 mins per hour
        sec=(hr-int(hr))*60		# seconds: 60 secs per minute
        start3=start2[:]
        start3[0]+=int(y)			
        start3[1]+=int(m)
        start3[2]+=int(d)
        start3[3]+=int(mn)
        start3[4]+=int(hr)
        start3[5]+=int(sec)
        start4.append(start3)

# convert to a time format
    time=[]
    for i in start4:
        time.append(datetime.datetime(i[0],i[1],i[2],i[3],i[4],i[5]))
    del y,m,d,mn,hr,sec
    del start2, start3, start4
    return time

def convert_to_hours(word, array):
    start = word
    timebound = array 
    # convert the string to a list of integers
    start2=[]	# year,month,day,hour,minute,second
    start2.append(int(start[0:4]))
    start2.append(int(start[5:7]))
    start2.append(int(start[8:10]))
    start2.append(int(start[11:13]))
    start2.append(int(start[14:16]))
    start2.append(int(start[17:]))
    difference=[]
    for i in range(0,len(timebound)):
        difference.append(timebound[i]-start2[i])
    # add the hours passed since the forecast reference time
    hours=difference[0]*(24*360)+difference[1]*(12)+difference[2]*(30)+difference[3]*(24)+difference[4]*(60)+difference[5]*(60)
    # year: 360 days per year
    # month: 12 months per year
    # day: 30 days per month
    # hour: 24 hours per day
    # minute: 60 mins per hour
    # seconds: 60 secs per minute
    return hours

