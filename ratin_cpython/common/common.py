import numpy as np
from math import factorial
import scipy.signal

#Gaussian filter with convolution - faster and easier to handle
## Degree is equal to the number of values left and right of the central value 
## of the gaussian window:
##  ie degree=3 yields a window of length 7
## It uses normalized weights (sum of weights = 1)
## Based on: 
##    http://en.wikipedia.org/wiki/Gaussian_filter
##    http://en.wikipedia.org/wiki/Standard_deviation
##    http://en.wikipedia.org/wiki/Window_function#Gaussian_window
def smooth(array_in, degree=5):
    '''
    Gaussian smooth line using a window of specified degree (=half-length)
    '''

    degree = int(degree) #make sure it is of integer type
    n = 2*degree+1

    if degree <= 0:
        return array_in
        
    if type(array_in) == type(np.array([])) and len(array_in.shape)>1:
        array_in = array_in.flatten()
    
    array_in = list(array_in)
    # If degree is larger than twice the original data, make it smaller
    if len(array_in) < n:
        degree = len(array_in)/2
        n = 2*degree+1
        print "Changed smoothing degree to:",degree    
    #extend the array's initial and ending values with equal ones, accordingly
    array_in = np.array( [array_in[0]]*degree + array_in + [array_in[-1]]*degree )
    
    #TODO: These parameters are subject to change - depends on the implementation
    # Gaussian parameters:
    x = np.linspace(-degree,degree,n) 
    sigma = np.sqrt( sum( (x-np.mean(x))**2 ) / n )
    alpha = 1.0 / (2.0 * sigma**2)
    weight = np.sqrt(alpha/np.pi) * np.exp(-alpha*x**2 )  #gaussian
    weights = weight / sum(weight)   #normalize
    
    return np.convolve(array_in, weights, 'valid')
    
    
#TODO: revise
#Gaussian 2D smoothing, anisotropic
## http://homepages.inf.ed.ac.uk/rbf/HIPR2/gsmooth.htm
def smooth2D(matrix_in, fill, degree=5, sigma=2.0, a=1.0, b=1.0):
    '''
    Gaussian smooth matrix using a window of specified degree
    '''
    kx, ky = np.arange(-degree,degree+1.0),np.arange(-degree,degree+1.0)
    kernel = np.zeros([kx.shape[0],ky.shape[0]])
    for i in range(len(kx)):
        for j in range(len(ky)):
            kernel[i,j] = 1./(2*np.pi*sigma**2) * np.exp( -(b*kx[i]**2+a*ky[j]**2)/(2*sigma**2) )
    kernel /= kernel.sum()
    
    matrix_out = scipy.signal.convolve2d(matrix_in, kernel, mode='same', fillvalue=fill)
    return matrix_out



def get_direction(x, y, smoothdegree=0, units='degrees'):
    '''
    Return direction (cartesian reference) of point
    The direction of each point is calculated as the mean of directions
    on both sides
    '''
    #Calculate direction in RADIANS
    direction = np.array([])
    #first point: Can determine direction only based on next point
    direction = np.append(direction,np.angle((x[1]-x[0])+(y[1]-y[0])*1j))
    for j in range(1, len(x)-1):
        # Base direction on points before and after current point
        direction = np.append(direction,np.angle((x[j+1]-x[j-1])+(y[j+1]-y[j-1])*1j))
    #last point: Can determine direction only based on previous point
    direction = np.append(direction,np.angle((x[-1]-x[-2])+(y[-1]-y[-2])*1j))

    #fix 'jumps' in data
    direction = fix_angle_vector(direction)
    
    #Smoothing - do not perform if input degree is equal/less than 0.0
    if smoothdegree <= 0.0:
        pass
    else:
        direction = smooth(direction, degree=smoothdegree)

    #TODO: Review! Do we need to confine it?
    #Limit the representation in the space of [0,2*pi]
    gaps = np.where(np.abs(direction) > np.radians(360.0))[0]
    direction[gaps] -= np.radians(360.0)
    
    if units=='radians':
        pass
    elif units == 'degrees':
        direction = np.degrees(direction)
        
    return direction
    
    
def distance(p1, p2):
    """
    Distance in between two points (given as tuples)
    """
    dist = np.sqrt( (p2[0]-p1[0])**2 + (p2[1]-p1[1])**2 )
    return dist
    
    
def distance_matrix(x0, y0, x1, y1, aniso):
    """
    Returns distances between points in a matrix formation.
    An anisotropy factor is set as input. If >1, the points in
    x direction shift closer. If <1, the points in x direction
    shift further apart. If =1, normal distances are computed.
    """
    aniso = float(aniso)
    x0 = np.array(x0).flatten()
    y0 = np.array(y0).flatten()
    x1 = np.array(x1).flatten()
    y1 = np.array(y1).flatten()
    #transpose observations
    vertical = np.vstack((x0, y0)).T
    horizontal = np.vstack((x1, y1)).T
    # Make a distance matrix between pairwise observations
    # Note: from <http://stackoverflow.com/questions/1871536>
    if aniso<=0.0:
        print "Warning: Anisotropy factor cannot be 0 or negative; set to 1.0."
        aniso = 1.0
    d0 = np.subtract.outer(vertical[:,0], horizontal[:,0]) * (1./aniso)
    d1 = np.subtract.outer(vertical[:,1], horizontal[:,1])
    return np.hypot(d0, d1)
    
    
#retrieve s values streamwise
def get_chainage(x, y):
    """
    Get chain distances for a set of continuous points
    """
    s = np.array([0.0]) #start
    for j in range(1,len(x)):
        s = np.append( s, s[j-1] + distance([x[j-1],y[j-1]], [x[j],y[j]]) )
    return s
    
    
def to_sn(Gx, Gy):
    """
    Transform (Gx,Gy) Cartesian coordinates to flow-oriented ones (Gs,Gn), 
    where Gx and Gy stand for gridded x and gridded y, and Gs and Gn are their
    transformed counterparts.
    Gx,Gy,Gs,Gn are all numpy arrays in the form of matrices.
    """
    rows, cols = Gx.shape
    #find s-direction coordinates
    midrow = int(rows/2)
    c_x = Gx[midrow,:]
    c_y = Gy[midrow,:]
    Salong = get_chainage(c_x,c_y)
    #all s-direction points have the same spacing
    Gs = np.tile(Salong, (rows,1)) #"stretch" all longitudinals
    
    #find n-direction coordinates
    Gn = np.zeros([rows,cols])
    for j in range(cols): #for each column
        Gn[midrow::-1,j] = -get_chainage(Gx[midrow::-1,j],Gy[midrow::-1,j])
        Gn[midrow:,j] = get_chainage(Gx[midrow:,j],Gy[midrow:,j])
    return Gs, Gn


def to_grid(data, rows, cols):
    """
    Transform a list of data to a grid-like (matrix) form of specified shape
    """
    data = np.array(data).flatten()
    return data.reshape(rows,cols)

    
##??['Brute-force' way but works correctly]    
def fix_angle_vector(theta):
    '''
    Fixes a vector of angles (in radians) that show 'jumps' because of changes
    between 360 and 0 degrees
    '''
    thetadiff = np.diff(theta)
    gaps = np.where(np.abs(thetadiff) > np.radians(180))[0]
    while len(gaps)>0:
        gap = gaps[0]
        if thetadiff[gap]<0:
            theta[gap+1:] += np.radians(360)
        else:
            theta[gap+1:] -= np.radians(360)
        thetadiff = np.diff(theta)
        gaps = np.where(np.abs(thetadiff) > np.radians(180))[0]

    return theta



def get_parallel_line(x, y, direction, distance, units = 'degrees'):
    '''
    Create parallel lines for representation of MAT path.
    '''
    if units == 'degrees':
        direction = np.radians(direction)
    perpendicular_direction = np.array(direction)+0.5*np.pi
    xn = np.array(x)+np.array(distance)*np.array(np.cos(perpendicular_direction))
    yn = np.array(y)+np.array(distance)*np.array(np.sin(perpendicular_direction))

    return xn, yn
    
    
#http://wiki.scipy.org/Cookbook/SavitzkyGolay
def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int", msg)
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')
    
 