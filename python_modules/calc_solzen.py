import numpy as np
def solar_zen(n,lon,lat):
    Jstar=n-lon/360.0
    M=(357.5291+0.98560028*Jstar)%360
    Mrads=M*np.pi/180.
    C=1.9148*np.sin(Mrads)+0.02*np.sin(2*Mrads)+0.0003*np.sin(3*Mrads)
    lamb=(M+C+180.+102.9372)%360
    lambrads=lamb*np.pi/180.
    Jtrans=2451545.5+Jstar+0.0053*np.sin(Mrads)-0.0069*np.sin(2*lambrads)
    sind=np.sin(lambrads)*np.sin(23.44*np.pi/180.)
    cosd=np.sqrt(1-sind**2)
    cosz=sind*np.sin(lat*np.pi/180.)+cosd*np.cos(lat*np.pi/180.)
    return np.real((cosz+0j)**2.42)
