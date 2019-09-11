# calculates the area of a regular lon-lat gridcell

def calc_area(lon,lat,radians=False):
    """
    Calculates the area of a gridcell on the surface of a sphere. Assumed to be in degrees.
    """
    # A = R^2 |sin(lat1)-sin(lat2)| |lon1-lon2| where R is earth radius (6371 kms)
    # A = R^2 * solid angle of latitude-longitude rectangle
    # polar rows (at +/-90) are treated differently.
    import numpy as np
    Pi=np.float128(3.141592653589793238462643383279)
    Earth_Radius=np.float128(6371.0*1.0E3)
    lat_bound = np.float128(89.999999999999999999999999)
    lon = np.float128(lon)
    lat = np.float128(lat)
    if (not radians):
        rlon = (lon[:]/np.float128(180.0))*Pi
        rlat = (lat[:]/np.float128(180.0))*Pi
    else:
        rlon = lon[:]
        rlat = lat[:]
    dlat = (rlat[1] - rlat[0])/2.0
    dlon = (rlon[1] - rlon[0])/2.0
    area = np.zeros((len(rlat),len(rlon)),np.float128)
    j=0
    while j < len(rlat):
        if (lat[j] >= lat_bound):
            lat1 = rlat[j]
            lat2 = rlat[j] - dlat/2.0
        elif (lat[j] <= -1.0*lat_bound):
            lat1 = rlat[j] + dlat/2.0            
            lat2 = rlat[j] 
        else:
            lat1 = rlat[j] + dlat
            lat2 = rlat[j] - dlat
        i=0
        while i < len(rlon):
            lon1 = rlon[i] - dlon
            lon2 = rlon[i] + dlon
            area[j,i] = (Earth_Radius**2)*(abs(np.sin(lat1)-np.sin(lat2))*abs(lon1-lon2))
            i += 1
# NEED TO MULTIPLY THE POLAR ROWS (I.E. WITH A VALUE AT +/-90) BY A FACTOR OF 4 TO 
# GET THE CORRECT VALUE... NOT SURE WHY
        if (lat[j] >= lat_bound) or (lat[j] <= -1.0*lat_bound):
            area[j,:] = area[j,:]*4.0
        j += 1
    ## sum area and print out
    #print area[1:len(rlat)-1,:].sum(),\
    #    area[:,:].sum(),4.0*Pi*(Earth_Radius**2),\
    #    ((4.0*Pi*(Earth_Radius**2))-area[:,:].sum())/(4.0*Pi*(Earth_Radius**2)),\
    #    (4.0*Pi*(Earth_Radius**2)-area[1:len(rlat)-1,:].sum())/(2.0*len(rlon)),\
    #    area[0,0],area[-1,0],\
    #    (4.0*Pi*(Earth_Radius**2)-area[1:len(rlat)-1,:].sum())/(2.0*len(rlon))/area[0,0]
    del Pi
    del Earth_Radius
    del lat_bound
    del lon
    del lat
    del rlon
    del rlat
    del dlon
    del dlat
    del i
    del j
    return area

def redistribute_integral(
    orig_grid,
    new_grid,
    orig_data,
    total=None):
    """
    Re-distributes data from one grid to the other by integrating along and then splitting that back off again.
    The data must be in kg or an equivalent unit.
    It is assumed that the orig_grid and new_grid have the '0' levels also defined prior to entering this routine,
    and that orig_data has also been correctly defined so that the 0th value has been correctly set up.
    """
    import numpy as np
    bottom_orig  = np.float128(orig_grid[0:len(orig_grid)-1])
    top_orig     = np.float128(orig_grid[1:len(orig_grid)])
    bottom_new   = np.float128(new_grid[0:len(new_grid)-1])
    top_new      = np.float128(new_grid[1:len(new_grid)])
    orig_data    = np.float128(orig_data)
    new_data     = np.zeros(len(new_grid)-1,dtype=np.float128)
    #print bottom_orig,top_orig
    # NOTE: new_data[0] will correspond with new_grid[1] etc. i.e. a +1 shift
    #       also len(orig_data) should equal len(orig_grid)
    if len(orig_data)+1 != len(orig_grid):
        print 'Grid = ',len(orig_grid),' ; Data = ',len(orig_data)
        raise Exception, 'Grid sizes incompatible'
    #print 'In redist'
    if total is None:
        total = orig_data.sum()
    k = 0
    while k < len(new_grid)-1:
        #print '----------------------------------'
        #print 'A',k,len(new_grid)-1,bottom_new[k],top_orig[-1],bottom_new[k] < top_orig[-1]
        if bottom_new[k] < top_orig[-1]:
            k1 = 0
            # make sure new grid is a little bit higher than the old grid, or same size
            while (k1 < len(top_orig) and top_orig[k1] <= bottom_new[k]):
                #print 'B',k1,len(top_orig),bottom_orig[k1],top_orig[k1],bottom_new[k],top_new[k],\
                #    (k1 < len(top_orig) and top_orig[k1] <= bottom_new[k])
                k1 += 1
            # now interpolate to the point where the grids cross
            new_data[k] = ((top_orig[k1] - bottom_new[k])/\
                               (top_orig[k1] - bottom_orig[k1]))*orig_data[k1]
            #print 'E',k1,k,new_data[k],orig_data[k1]
            # now integrate into new_data
            k2 = k1 + 1
            #print 'C',k,k1,k2,len(new_data),len(orig_data)
            #print 'G',k,k1,k2,len(new_data),len(orig_data),new_data[k],orig_data[k2],\
            #    k2,len(orig_data),bottom_orig[k2],top_new[-1],\
            #    k2 < len(orig_data),bottom_orig[k2] < top_new[-1],\
            #    (k2 < len(orig_data) and bottom_orig[k2] < top_new[-1])
            while (k2 < len(orig_data) and bottom_orig[k2] < top_new[-1]):
                new_data[k] = new_data[k] + orig_data[k2]
                #print 'F',k,k1,k2,len(new_data),len(orig_data),new_data[k],orig_data[k2],\
                #k2,len(orig_data),bottom_orig[k2],top_new[-1],\
                #k2 < len(orig_data),bottom_orig[k2] < top_new[-1],\
                #(k2 < len(orig_data) and bottom_orig[k2] < top_new[-1])
                k2 += 1
        k += 1
    # now split off integration, leaving final field
    #print '----------------------------------'
    k = 0
    while k < len(new_data)-1:
        #print 'D',k
        new_data[k] = new_data[k] - new_data[k+1]
        #print k,new_data[k],new_data[k+1],new_data[k]+new_data[k+1],len(new_data)
        k += 1
    # now re-scale to make sure mass is conserved
    if new_data.sum() != total:
        new_data = new_data*(total/new_data.sum())
    del k
    del k1
    del k2
    del bottom_orig
    del top_orig   
    del bottom_new 
    del top_new    
    del orig_data 
    del orig_grid
    del new_grid
    return new_data


def expand_longitude(
    test_lon,orig_lon,orig_data):
    """
    Makes a cyclic array using basemap and expands the longitudes to give the buffer aroud the edge needed 
    for the area weighted re-gridding.
    """
    import numpy as np
    from mpl_toolkits.basemap import addcyclic, shiftgrid
    import copy
    # shift the grid of the emissions so it fits with the test longitude grid
    # first, find the value of the orig long grid nearest to the test long grid
    if len(orig_data.shape) != 2:
        raise Exception, 'Data array must be two-dimensional'
    idx = (np.abs(orig_lon-test_lon[0])).argmin()
    start = orig_lon[idx]
    orig_data_temp,orig_lon_temp=addcyclic(orig_data,orig_lon)
    orig_data_temp,orig_lon_temp=shiftgrid(start,orig_data_temp,orig_lon_temp)
    test_dlon = test_lon[1] - test_lon[0]
    orig_dlon = orig_lon[1] - orig_lon[0]
    extra_cells=np.int(test_dlon/orig_dlon)+2
    # NOTE: DATA AND LOGITUDE HAVE BEEN MADE CYCLIC, ORIGIN IS THAT CLOSEST IN ORIG TO THE START OF TEST
    new_lon   = np.linspace(orig_lon_temp[0] - extra_cells*orig_dlon,\
                              orig_lon_temp[-2]+(extra_cells)*orig_dlon,\
                              len(orig_lon)+2*extra_cells)
    lon_index = np.arange(-1*extra_cells,len(new_lon)-extra_cells,1,dtype=np.int)
    new_data = np.zeros((orig_data.shape[0],len(new_lon)))
    i = 0
    while i < len(new_lon):
        # ignore cyclic grid here for length, since orig_data_temp[0] == orig_data_temp[1], but
        # we don't want to replicate that
        new_data[:,i] = copy.deepcopy(orig_data_temp[:,lon_index[i]%orig_data.shape[1]])
        i+=1
    del idx
    del i
    del start
    del orig_data_temp
    del orig_lon_temp
    del test_dlon
    del orig_dlon
    del test_lon
    del orig_lon
    del orig_data
    del extra_cells
    del lon_index
    return new_lon,new_data



