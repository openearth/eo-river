# -*- coding: utf-8 -*-
"""
Created on Wed May 20 11:26:56 2015

@author: zervakis
"""

import numpy as np
import matplotlib.pyplot as plt
import shapefile
from scipy.interpolate import griddata, Rbf
from pykrige.ok import OrdinaryKriging
from pykrige.uk import UniversalKriging
#from scipy.spatial import cKDTree as tree
from copy import deepcopy
from ratin.common import common
from shapely.geometry import Polygon, LineString, Point
from scipy.spatial import Voronoi
import gc


class Interpol():
    
    def __init__(self, locations, samples, skipped_points=[0,0]) :
        """
        Interpol class requires the grid locations (x,y) and the sampled data
        depths as a numpy matrix with unsampled points set to -999.0.
        Optionally the [right,left] number of cells from the river banks are 
        given to exclude from interpolation. These cells are instead filled in
        with a linear bank profile.
        
        Args:
             - locations : <dic> of 'x' and 'y' numpy arrays (matrix notation)
             - samples : <numpy.array> (in s-by-n matrix notation)
        
        Kwargs:
             - skipped_points : <list> of 2 <int> ([left,right])
        """
        try:
            #Transform x,y coordinates to s,n flow-oriented coordinates:
            self.x = locations['x']
            self.y = locations['y']
            self.rows, self.cols = self.x.shape
            s,n = common.to_sn(self.x, self.y)
        except KeyError:
            print "Error: Input locations should be given as a dictionary of {\'x\',\'y\'}!"
            raise
        
        #check
        if self.rows%2 != 1: #number of rows must be odd
            print "Error: Input locations must have an odd number of transverse coordinates!"
            raise
        
        #samples check
        if type(samples) is not type(np.array([])):
            print "ERROR: Sample data are not of matrix notation."
            raise
        elif samples.shape != s.shape:
            print "ERROR: Sample data do not match the provided grid locations."
            raise
        else:
            samples = deepcopy(samples) #avoid overwrite
        
        #Flatten matrix notation to vector (holds structure)
        self.s = s.flatten()
        self.n = n.flatten()
        samples[samples==-999.0] = np.nan #exchange notation of "no value"
        self.z = samples.flatten()
        
        self.isdata = ~np.isnan(self.z) #where sampled data exists
        self.nodata = np.isnan(self.z)  #where data is missing
        
        #number of grid points to be filled in with bank profile
        self.skip = [abs(int(skipped_points[0])),abs(int(skipped_points[1]))]
        #TODO: assumed water level // maybe exclude it ?
        self.wl = np.nanmax(self.z) 
        
        self.z_interpol = None #to hold the interpolated values


    #Representation
    def plot(self, fignum, clim = 'auto', classes = 100, system = 'xy', hold = False):
        
        #append nan values
        Gz = deepcopy(self.z_interpol)
        Gz[Gz==-999.0] = np.nan
        
        if system == 'xy':
            Gx = self.x
            Gy = self.y
        elif system == 'sn':
            Gx = common.to_grid(self.s, self.rows, self.cols)
            Gy = common.to_grid(self.n, self.rows, self.cols)
        elif system == 'mn':
            Gx, Gy = np.meshgrid(self.cols, self.rows)
        else:
            print "Wrong coordinate system chosen. Please specify input system as \'xy\', \'sn\' or \'mn\'."

        #plt.close(fignum)
        plt.figure(fignum)
        #plt.clf()
        
        # set colorlimits
        minn= np.nanmin(Gz)
        maxx= np.nanmax(Gz)
        if clim == 'auto':
            mean = np.mean(np.ma.masked_array(Gz,np.isnan(Gz)))
            classes = int(classes)
            if classes%2==1:
                classes+=1
            dl = np.linspace(minn,mean,classes/2)
            dr = np.linspace(mean,maxx,classes/2)
            d = 10.0
            flag = True
            c = 0
            while flag:
                levels = np.round(np.append(dl,dr[1:])/d) * d
                if c%2 == 0:
                    d = d/2.0
                else:
                    d = d - 4*d/5.0
                c += 1
                flag = len(levels)!=len(set(levels))
        elif clim == 'minmax':
            levels = np.linspace(minn,maxx,classes)
        else:
            levels = np.linspace(clim[0],clim[1],classes)

        #If more than half the points are unknown, the dataset to plot logically must be a testing dataset
        if np.sum(np.isnan(Gz)) > Gz.shape[0]*Gz.shape[1]/2.0:
            print "-> Assuming plot of Testing Dataset."
            cs = plt.scatter(Gx, Gy, s=10, c=Gz, cmap=plt.cm.jet, vmin=clim[0], vmax=clim[1], edgecolors='none')
        else:
            cs = plt.contourf(Gx, Gy, Gz, levels, cmap=plt.cm.jet, extend="both")
        cs.cmap.set_under('k')
        cs.set_clim(np.floor(levels[0]), np.ceil(levels[-1]))
        cbar = plt.colorbar(cs)
        cbar.set_label('Bed Level (m)', labelpad=15, weight='bold', size='14')
        plt.axis('equal')
        plt.hold(hold) #used when the data needs to be "on top" of other data
        plt.ion()
        plt.show()        
    
    
#    #TODO: Possibly delete
#    #unused, the initital thought was to smoothen the values along a a longitudinal line
#    #making use of only the previous values (a window of half-size, reaching back)
#    def fwdsmooth(self, array_in, degree=5):
#        degree = int(degree) #make sure it is of integer type
#        n = degree+1
#        if degree <= 0:
#            return array_in
#        array_in = list(array_in)
#        # If degree is larger than twice the original data, make it smaller
#        if len(array_in) < n:
#            degree = len(array_in)/2
#            print "Changed smoothing degree to:",degree    
#        #extend the array's initial values with equal ones
#        array_in = np.array( [array_in[0]]*degree + array_in )
#        
#        # Gaussian parameters:
#        x = np.linspace(0.0,degree,n) 
#        sigma = np.sqrt( sum( (x-np.mean(x))**2 ) / n )
#        alpha = 1.0 / (2.0 * sigma**2)
#        weight = np.sqrt(alpha/np.pi) * np.exp(-alpha*x**2 )  #gaussian
#        weights = weight / sum(weight)   #normalize
#        
#        return np.convolve(array_in, weights, 'valid')    

        
    #EIDW - interpolator
    def eidw_interpol(self, x, y, z, xi, yi, n, ar, radius):
        gc.enable()
        dist = common.distance_matrix(x,y, xi,yi, ar)
        #global interpolator
        if radius == 'None' or radius <= 0.0:
            # inverse-distances = weights are 1 / distance
            weights = 1.0 / dist**n
            # Make weights sum to one
            weights /= weights.sum(axis=0)
            # Multiply the weights for each interpolated point by all observed Z-values
            zi = np.dot(weights.T, z)
        #local interpolator
        else:
            dist[dist > radius] = np.nan
            # weights are 1 / distance
            weights = 1.0 / dist**n
            zi = np.zeros(len(xi))
            for i in range(weights.shape[1]):
                w = weights[:,i] / np.nansum(weights[:,i])
                zi[i] = np.nansum(w * z)
        del dist
        return zi
        
        
    #Solver function for Radial Basis interpolation
    def rbf_solver(self, x, y, z, xi, yi, ar, method='linear'):
        x *= 1./ar
        xi *= 1./ar
        if method == 'linear': #Linear solution of rbf matrix
            dist = common.distance_matrix(x,y, xi,yi, ar)
            internal_dist = common.distance_matrix(x,y, x,y, ar)
            weights = np.linalg.solve(internal_dist, z)
            zi =  np.dot(dist.T, weights)
        elif method == 'tps': #Thin-Plate solution of rbf matrix
            rbfi = Rbf(x,y,z,function='thin_plate')
            zi = rbfi(xi,yi)
        else:
            print "ERROR: Wrong definition for method of RBF interpolation!"
            raise
        return zi



###############################################################################
################################ INTERPOLATION ################################
###############################################################################

    #Inverse Distance Weighting (IDW) [global interpolator]
    def idw(self, power = 2.0, smoothie = 0, plane = 'sn', plot = False):
        """
        Performs Inverse Distance Weighting (IDW) interpolation to the instance's data
        
        Kwargs: power <float>  : power parameter, should be above 2.0
                smoothie <int> : smoothing degree (Gaussian window)
                plane <str>    : chooses plane for interpolation; 'xy' for Cartesian, 'sn' for flow-oriented
                plot <boolean> : decides to plot or not the resulting IDW values
        """
        gc.enable()        
        print "Calculating: Inverse Distance Weighting"
        #check
        if power <2.0:
            print "Warning: Power parameter too small; set to 2.0."
            power = 2.0
            
        idw = deepcopy(self.z) #avoid overwrite
        
        #calculate IDW (when anisotropy ratio is 1, we have IDW; radius='None' executes global IDW)        
        if plane == 'sn':

            idw[self.nodata] = self.eidw_interpol(self.s[self.isdata],self.n[self.isdata],idw[self.isdata], 
                                                  self.s[self.nodata],self.n[self.nodata], power, 1.0, 'None')
        elif plane == 'xy':
            x = self.x.flatten()
            y = self.y.flatten()
            idw[self.nodata] = self.eidw_interpol(x[self.isdata],y[self.isdata],idw[self.isdata], 
                                                  x[self.nodata],y[self.nodata], power, 1.0, 'None')
        else:
            print "Error: Plane for interpolation not correctly specified. No interpolation performed."
            return
                                              
        #grid (matrix) notation
        idw = common.to_grid(idw, self.rows, self.cols)
        #add pre-defined banks
        idw = self.__add_banks(idw)
        
        #smoothen
        for i in range(self.rows):
            idw[i,:] = common.smooth(idw[i,:],smoothie)
        print "Finished [IDW]."
        self.z_interpol = idw
        del idw
        #plot
        if plot:
            self.plot(fignum='IDW-GLOBAL')
    
    
    #Elliptical Inverse Distance Weighting (EIDW) [global interpolator]
    def eidw(self, power = 2.0, anisotropy = 5.0, smoothie = 0, plane = 'sn', plot = False):
        """
        Performs Elliptical Inverse Distance Weighting (EIDW) interpolation to the instance's data
        
        Kwargs: power <float>      : power parameter, should be above 2.0
                anisotropy <float> : anisotropy parameter: if >1.0, it brings points closer in the
                                     longitudinal (s) direction, if <1.0 it brings points closer 
                                     in the transverse (n) direction
                smoothie <int>     : smoothing degree (Gaussian window)
                plane <str>    : chooses plane for interpolation; 'xy' for Cartesian, 'sn' for flow-oriented
                plot <boolean> : decides to plot or not the resulting EIDW values
        """
        gc.enable()
        print "Calculating: Elliptical Inverse Distance Weighting [EIDW]"
        #check
        if power <2.0:
            print "Warning: Power parameter too small; set to 2.0."
            power = 2.0
        
        eidw = deepcopy(self.z) #avoid overwrite
        
        #calculate EIDW (radius='None' executes global EIDW)
        if plane == 'sn':            
            eidw[self.nodata] = self.eidw_interpol(self.s[self.isdata],self.n[self.isdata],eidw[self.isdata], 
                                                   self.s[self.nodata],self.n[self.nodata], power, anisotropy, 'None')
        elif plane == 'xy':
            x = self.x.flatten()
            y = self.y.flatten()
            eidw[self.nodata] = self.eidw_interpol(x[self.isdata],y[self.isdata],eidw[self.isdata], 
                                                   x[self.nodata],y[self.nodata], power, anisotropy, 'None')
        else:
            print "Error: Plane for interpolation not correctly specified. No interpolation performed."
            return
        
        #grid (matrix) notation
        eidw = common.to_grid(eidw, self.rows, self.cols)
        #add pre-defined banks
        eidw = self.__add_banks(eidw)
        #smoothen
        for i in range(self.rows):
            eidw[i,:] = common.smooth(eidw[i,:],smoothie)
        print "Finished [EIDW]."
        self.z_interpol = eidw
        del eidw
        #plot
        if plot:
            self.plot(fignum='EIDW-GLOBAL')

    
    #LINEAR INTERPOLATION
    def linear(self, smoothie = 0, plane = 'sn', fill_nans = True, plot = False):
        """
        Performs Linear Barycentric interpolation to the instance's data
        
        Kwargs: smoothie <int>  : smoothing degree (Gaussian window)
                plane <str>     : chooses plane for interpolation; 'xy' for Cartesian, 'sn' for flow-oriented
                fill_nans <str> : decides to fill NaN values in a linear manner or keep NaN values
                plot <boolean>  : decides to plot or not the resulting LINEAR values
        """
        gc.enable()
        print "Calculating: Linear Barycentric Interpolation"
        linear = deepcopy(self.z) #avoid overwrite
        
        #choice of plane
        if plane == 'sn':
            pnts = np.array(zip(self.s,self.n))
        elif plane == 'xy':
            pnts = np.array(zip(self.x.flatten(),self.y.flatten()))
        else:
            print "Error: Plane for interpolation not correctly specified. No interpolation performed."
            return
            
        #linear interpolation
        linear[self.nodata] = griddata(pnts[self.isdata],linear[self.isdata],pnts[self.nodata],method='linear',fill_value=np.nan)
        linear = common.to_grid(linear,self.rows,self.cols) #griddify
        if fill_nans:
            linear = self.__fill_nans(linear) #linear fill of nan values
        
        #add pre-defined banks
        linear = self.__add_banks(linear)
        #smoothen
        for i in range(self.rows):
            linear[i,:] = common.smooth(linear[i,:],smoothie)
        print "Finished [LINEAR]."
        self.z_interpol = linear
        del linear
        #plot
        if plot:
            self.plot(fignum='LINEAR-BARYCENTRIC')

    
    #CUBIC INTERPOLATION
    #WARNING!CAN HAVE TOO MANY COMPUTATIONS IN (X,Y) PLANE (MEMORY ERROR)
    def cubic(self, smoothie = 0, plane = 'sn', fill_nans = True, plot = False):
        """
        Performs Cubic interpolation to the instance's data (piecewise cubic,
        continuously differentiable and approximately curvature-minimizing 
        polynomial surface).
        WARNING: High memory usage, usually in \'xy\' plane
        
        Kwargs: smoothie <int>  : smoothing degree (Gaussian window)
                plane <str>     : chooses plane for interpolation; 'xy' for Cartesian, 'sn' for flow-oriented
                fill_nans <str> : decides to fill NaN values in a linear manner or keep NaN values
                plot <boolean>  : decides to plot or not the resulting CUBIC values
        """
        gc.enable()
        print "Calculating: Cubic Interpolation"
        cubic = deepcopy(self.z) #avoid overwrite
        
        #choice of plane
        if plane == 'sn':
            pnts = np.array(zip(self.s,self.n))
        elif plane == 'xy':
            pnts = np.array(zip(self.x.flatten(),self.y.flatten()))
        else:
            print "Error: Plane for interpolation not correctly specified. No interpolation performed."
            return
        
        #cubic interpolation
        cubic[self.nodata] = griddata(pnts[self.isdata],cubic[self.isdata],pnts[self.nodata],method='cubic',fill_value=np.nan)
        cubic = common.to_grid(cubic,self.rows,self.cols) #griddify
        if fill_nans:
            cubic = self.__fill_nans(cubic) #linear fill of nan values
        
        #add pre-defined banks
        cubic = self.__add_banks(cubic)
        #smoothen
        for i in range(self.rows):
            cubic[i,:] = common.smooth(cubic[i,:],smoothie)
        print "Finished [CUBIC]."
        self.z_interpol = cubic
        del cubic
        #plot
        if plot:
            self.plot(fignum='CUBIC')
        
    
    #NEAREST-NEIGHBOUR INTERPOLATION (not really an interpolation, more like value assignment)
    def near(self, smoothie = 0, plane = 'sn', plot = False):
        """
        Fills values of unsampled data points with the values of the closest
        sampled data point of the instance's data.
        
        Kwargs: smoothie <int>  : smoothing degree (Gaussian window)
                plane <str>     : chooses plane for interpolation; 'xy' for Cartesian, 'sn' for flow-oriented
                plot <boolean>  : decides to plot or not the result NEAREST values
        """
        gc.enable()
        print "Calculating: Nearest-Neighbour \"Interpolation\""
        near = deepcopy(self.z) #avoid overwrite
        
        #choice of plane
        if plane == 'sn':
            pnts = np.array(zip(self.s,self.n))
        elif plane == 'xy':
            pnts = np.array(zip(self.x.flatten(),self.y.flatten()))
        else:
            print "Error: Plane for interpolation not correctly specified. No interpolation performed."
            return
        
        #find nearest values
        near[self.nodata] = griddata(pnts[self.isdata],near[self.isdata],pnts[self.nodata],method='nearest')
        near = common.to_grid(near,self.rows,self.cols) #griddify
        
        #add pre-defined banks
        near = self.__add_banks(near)
        #smoothen
        for i in range(self.rows):
            near[i,:] = common.smooth(near[i,:],smoothie)
        print "Finished [NEAR]."
        self.z_interpol = near
        del near
        #plot
        if plot:
            self.plot(fignum='NEAREST')
            
            
    #Radial Basis Function (RBF)
    def rbf(self, method = 'linear', anisotropy = 5.0, smoothie = 0, plane = 'sn', plot = False):
        """
        Performs Radial Basis Interpolation using the defined method (linear or
        thin plate spline)
        
        Kwargs: method <str>        : choice of RBF interpolation ['linear' or 'tps']
                anisotropy <float>  : anisotropy parameter: if >1.0, it brings points closer in the
                                      longitudinal (s) direction, if <1.0 it brings points closer 
                                      in the transverse (n) direction
                smoothie <int>      : smoothing degree (Gaussian window)
                plane <str>         : chooses plane for interpolation; 'xy' for Cartesian, 'sn' for flow-oriented
                plot <boolean>      : decides to plot or not the resulting RBF values
        """
        gc.enable()
        print "Calculating: Radial Basis Function (RBF) with {0} method".format(method)
        rbf = deepcopy(self.z) #avoid overwrite
        
        #choice of plane and rbf interpolation
        if plane == 'sn':
            rbf[self.nodata] = self.rbf_solver(self.s[self.isdata], self.n[self.isdata], rbf[self.isdata], self.s[self.nodata], self.n[self.nodata], anisotropy, method)
        elif plane == 'xy':
            x = self.x.flatten()
            y = self.y.flatten()
            rbf[self.nodata] = self.rbf_solver(x[self.isdata], y[self.isdata], rbf[self.isdata], x[self.nodata], y[self.nodata], anisotropy, method)
        else:
            print "Error: Plane for interpolation not correctly specified. No interpolation performed."
            return
            
        rbf = common.to_grid(rbf,self.rows,self.cols) #griddify
        
        #add pre-defined banks
        rbf = self.__add_banks(rbf)
        #smoothen
        for i in range(self.rows):
            rbf[i,:] = common.smooth(rbf[i,:],smoothie)
        print "Finished [RBF]."
        self.z_interpol = rbf
        del rbf
        #plot
        if plot:
            self.plot(fignum='RADIAL-BASIS[{0}]'.format(method.upper()))
            
    
    #Anisotropic Ordinary Kriging
    def ok(self, sill, vmodel = 'spherical', anisotropy = 1.0, smoothie = 0, plane = 'sn', plot = False):
        """
        Performs Ordinary Kriging (OK) interpolation with defined anisotropy 
        to the instance's data using the specified model for its semivariogram.
        The sill for the semivariogram needs to be defined. The range used is
        equal to the sill times the anisotropy defined.
        
        Args:
                sill <float>        : Sill value for semivariogram in OK
        
        Kwargs: vmodel <str>        : Semivariogram model (default: 'spherical')
                anisotropy <float>  : anisotropy parameter: if >1.0, it brings points closer in the
                                      longitudinal (s) direction, if <1.0 it brings points closer 
                                      in the transverse (n) direction
                smoothie <int>      : smoothing degree (Gaussian window)
                plane <str>         : chooses plane for interpolation; 'xy' for Cartesian, 'sn' for flow-oriented
                plot <boolean>      : decides to plot or not the resulting OK values
        """
        gc.enable()
        print "Calculating: Ordinary Kriging [OK]"
        
        ok = deepcopy(self.z) #avoid overwrite

        #calculate OK with anisotropy defined
        if plane == 'sn':            
            OK = OrdinaryKriging(self.s[self.isdata], self.n[self.isdata], ok[self.isdata], 
                                 anisotropy_scaling = anisotropy, anisotropy_angle=0.0, nlags=self.cols,
                                 variogram_model=vmodel, variogram_parameters = [sill,sill*anisotropy,0.0], 
                                 verbose=True, enable_plotting=True)
            ok[self.nodata], sigmas = OK.execute('points', self.s[self.nodata], self.n[self.nodata])
        elif plane == 'xy':
            x = self.x.flatten()
            y = self.y.flatten()
            OK = OrdinaryKriging(x[self.isdata], y[self.isdata], ok[self.isdata], 
                                 anisotropy_scaling = anisotropy, anisotropy_angle=0.0, nlags=self.cols,
                                 variogram_model=vmodel, variogram_parameters = [sill,sill*anisotropy,0.0], 
                                 verbose=False, enable_plotting=False)
            ok[self.nodata], sigmas = OK.execute('points', x[self.nodata], y[self.nodata])
        else:
            print "Error: Plane for interpolation not correctly specified. No interpolation performed."
            return
        
        #grid (matrix) notation
        ok = common.to_grid(ok, self.rows, self.cols)
        #add pre-defined banks
        ok = self.__add_banks(ok)
        #smoothen
        for i in range(self.rows):
            ok[i,:] = common.smooth(ok[i,:],smoothie)
        print "Finished [OK]."
        self.z_interpol = ok
        del ok
        #plot
        if plot:
            self.plot(fignum='OK')
            
    #TODO: REVIEW!!
    def uk_funct(self,x,y):
        return x+y**2
            
    #Anisotropic Universal Kriging with f(x,y) = x+y^2
    #TODO: Probably requires fixing
    def uk(self, sill, vmodel = 'spherical', anisotropy = 1.0, smoothie = 0, plane = 'sn', plot = False):
        """
        Performs Universal Kriging (UK) interpolation with defined anisotropy 
        to the instance's data using the specified model for its semivariogram.
        The sill for the semivariogram needs to be defined. The range used is
        equal to the sill times the anisotropy defined.
        
        Args:
                sill <float>        : Sill value for semivariogram in UK
        
        Kwargs: vmodel <str>        : Semivariogram model (default: 'spherical')
                anisotropy <float>  : anisotropy parameter: if >1.0, it brings points closer in the
                                      longitudinal (s) direction, if <1.0 it brings points closer 
                                      in the transverse (n) direction
                smoothie <int>      : smoothing degree (Gaussian window)
                plane <str>         : chooses plane for interpolation; 'xy' for Cartesian, 'sn' for flow-oriented
                plot <boolean>      : decides to plot or not the resulting UK values
        """
        gc.enable()
        print "Calculating: Universal Kriging [UK]"
        
        uk = deepcopy(self.z) #avoid overwrite

        #calculate UK with anisotropy defined
        if plane == 'sn':            
            UK = UniversalKriging(self.s[self.isdata], self.n[self.isdata], uk[self.isdata], 
                                 anisotropy_scaling = anisotropy, anisotropy_angle=0.0, nlags=self.cols,
                                 variogram_model=vmodel, variogram_parameters = [sill,sill*anisotropy,0.0], 
                                 drift_terms=['functional'], functional_drift=[self.uk_funct],
                                 verbose=True, enable_plotting=True)
            uk[self.nodata], sigmas = UK.execute('points', self.s[self.nodata], self.n[self.nodata])
        elif plane == 'xy':
            x = self.x.flatten()
            y = self.y.flatten()
            UK = UniversalKriging(x[self.isdata], y[self.isdata], uk[self.isdata], 
                                 anisotropy_scaling = anisotropy, anisotropy_angle=0.0, nlags=self.cols,
                                 variogram_model=vmodel, variogram_parameters = [sill,sill*anisotropy,0.0], 
                                 drift_terms=['regional_linear'], functional_drift=[self.uk_funct],
                                 verbose=False, enable_plotting=False)
            uk[self.nodata], sigmas = UK.execute('points', x[self.nodata], y[self.nodata])
        else:
            print "Error: Plane for interpolation not correctly specified. No interpolation performed."
            return
        
        #grid (matrix) notation
        uk = common.to_grid(uk, self.rows, self.cols)
        #add pre-defined banks
        uk = self.__add_banks(uk)
        #smoothen
        for i in range(self.rows):
            uk[i,:] = common.smooth(uk[i,:],smoothie)
        print "Finished [UK]."
        self.z_interpol = uk
        del uk
        #plot
        if plot:
            self.plot(fignum='UK')
            
            
    #TODO: PROBABLY ONLY POSSIBLE TO PERFORM IN SN system!
    #Natural Neighbour interpolation [brute-force, bounded by river polygon]
    #Warning: requires a lot of execution time
    def natneigh(self, anisotropy=1.0, smoothie = 0, plane = 'sn', plot = False):
        """
        Performs Natural Neighbor (NN) interpolation with defined anisotropy 
        to the instance's data. The algorithm assumes a brute-force way of 
        calculating the Voronoi cells every time an unsampled point's value is
        computed.
        
        Kwargs: anisotropy <float>  : anisotropy parameter: if >1.0, it brings points closer in the
                                      longitudinal (s) direction, if <1.0 it brings points closer 
                                      in the transverse (n) direction
                smoothie <int>      : smoothing degree (Gaussian window)
                plane <str>         : chooses plane for interpolation; 'xy' for Cartesian, 'sn' for flow-oriented
                plot <boolean>      : decides to plot or not the resulting NN values
        """
        
        from scipy.spatial import voronoi_plot_2d
        import matplotlib.pyplot as plt        
        
        gc.enable()
        print "Calculating: Natural Neighbour [NN]"
        
        nn = deepcopy(self.z) #avoid overwrite

        #choose plane
        if plane == 'sn':
            x = self.s * 1./anisotropy  #already flattened
            y = self.n                  #already flattened
        elif plane == 'xy':
            x = self.x.flatten() * 1./anisotropy
            y = self.y.flatten()
        else:
            print "Error: Plane for interpolation not correctly specified. No interpolation performed."
            return
        
        #DEFINE BOUNDARY FOR INTERPOLATION
        #griddify points to choose boundary easier:
        Gx = common.to_grid(x,self.rows,self.cols)
        Gy = common.to_grid(y,self.rows,self.cols)
        bx = np.hstack( (Gx[0,:],Gx[1:-1,-1],Gx[-1,:][::-1],Gx[1:-1,0][::-1]) )
        by = np.hstack( (Gy[0,:],Gy[1:-1,-1],Gy[-1,:][::-1],Gy[1:-1,0][::-1]) )
        #define boundary:
        boundary = np.array(zip(bx,by))
        
        #VORONOI OF SAMPLED DATA POINTS:
        vorpoints = np.array(zip(x[self.isdata],y[self.isdata]))
        #shift points around central point of the dataset (for precision purposes)
        center = vorpoints.mean(axis=0)
        vorpoints -= center
        #construct Voronoi diagram from (centered) sampled data points
        vor = Voronoi(vorpoints)
        vor.close()
        """
        #TODO: delete:
        voronoi_plot_2d(vor)
        plt.ion()
        plt.axis('equal')
        plt.show()
        """
        #calculate areas of sampled dataset Voronoi cells
        original_areas,vor = self.__find_areas(vor,boundary-center)
        
        """
        #TODO: delete:        
        # colorize
        for region in vor.regions[1:]:
            polygon = vor.vertices[region]
            plt.fill(*zip(*polygon), alpha=0.4)
        plt.plot(vorpoints[:,0], vorpoints[:,1], 'ko')
        plt.xlim(vor.min_bound[0] - 0.1, vor.max_bound[0] + 0.1)
        plt.ylim(vor.min_bound[1] - 0.1, vor.max_bound[1] + 0.1)
        plt.axis('equal')
        plt.ion()
        plt.show()
        """
        
        #ITERATE THROUGH UNSAMPLED DATA POINTS:
        # ~ For each unsampled point, construct the new Voronoi diagram consisting of the 
        # ~ sampled data and the current unsampled point. Then calculate the new areas and
        # ~ find the normalized weights based on how much of each area is "stolen away" from
        # ~ each sampled dataset Voronoi cell (https://en.wikipedia.org/wiki/Natural_neighbor).
        # ~ The areas are always bounded by the river polygon (based on the grid defined).
        unknown = []
        for i in range(len(x[self.nodata])):
            if i%1000==0:
                print i, "out of ", len(x[self.nodata])
            #add new point shifted around central point
            varys = np.vstack( (vorpoints,[x[self.nodata][i]-center[0],y[self.nodata][i]-center[1]]) )
            
            #calculate new Voronoi
            pntvor = Voronoi(varys)
            pntvor.close()
            
            #calculate areas
            new_areas,pntvor = self.__find_areas(pntvor,boundary-center)
            new_areas = new_areas[:-1] #exclude new point's area
            w = new_areas / original_areas
            w[w>1.0] = 1.0 #make sure that no area is larger than initial areas
            areaweight = 1.0 - w
            normalize = np.nansum(areaweight)
            if normalize == 0.0: #to avoid division by 0
                normalize = 1e-12
            areaweight /=  normalize
            unknown.append( np.nansum(areaweight*self.z[self.isdata]) )
        
        nn[self.nodata] = np.array(unknown)
        
        #grid (matrix) notation
        nn = common.to_grid(nn, self.rows, self.cols)
        #add pre-defined banks
        nn = self.__add_banks(nn)
        #smoothen
        for i in range(self.rows):
            nn[i,:] = common.smooth(nn[i,:],smoothie)
        print "Finished [NN]."
        self.z_interpol = nn
        del nn
        #plot
        if plot:
            self.plot(fignum='NN')
        
    
    
    
    
    """
###############################################################################
##### LOCAL INTERPOLATION ? #########################
###############################################################################
    
    #EIDW - Local    
    #TODO: Fix local interpolator by distance/neighbourhoods
    #radius must be given. Radius describes the circle neighbourhood that would normally be taken into account in IDW, 
    #but the anisotropy ratio describes the ratio of major/minor axes centered at that circle.
    #This means that anisotropy>1.0 transforms the circle into an ellipse longer in longitudinal direction (s)
    #Anisotropy<1.0 should not be used, because it results in ellipses longer to the traverse direction (n)
    def l_eidw(self, radius, pidw=2, anisotropy=5.0, smoothie=0, banks=False, plot=False):
        print "Calculating: Elliptical Inverse Distance Weighting (local)"
    
        print "Local EIDW anisotropy ratio:", anisotropy
        
        x_to_s, y_to_n = self._to_sn(self.td)
        l_eidw = copy(self.td)
        leidw_nodata = copy(self.nodata)
        
        while np.sum(leidw_nodata) > 0:
            print "Still to calculate:",np.sum(leidw_nodata)
            x0 = x_to_s[~leidw_nodata]
            y0 = y_to_n[~leidw_nodata]
            z0 = l_eidw[:,2][~leidw_nodata]
            x1 = x_to_s[leidw_nodata]
            y1 = y_to_n[leidw_nodata]
            z1 = l_eidw[:,2][leidw_nodata]    
            pnts = np.vstack((x0,y0)).T
            interpnts = np.vstack((x1,y1)).T
            kdT = tree(pnts)
            for p in range(interpnts.shape[0]):
                neighs = kdT.query_ball_point(interpnts[p], radius*anisotropy)
                if neighs:
                    z1[p] = self.eidw_interpol(kdT.data[neighs][:,0], kdT.data[neighs][:,1], z0[neighs], 
                                               [interpnts[p][0]], [interpnts[p][1]], pidw, anisotropy, radius)
                else:
                    pass
            l_eidw[:,2][leidw_nodata] = z1
            if np.sum(np.isnan(l_eidw[:,2])) == np.sum(leidw_nodata):
                print "Doubling the max distance traverse flow."
                radius = radius*2   #max distance traverse flow
                anisotropy = 0.5*anisotropy  #anisotropy ratio
                print "New Anisotropy Ratio:", anisotropy
            leidw_nodata = np.isnan(l_eidw[:,2])
        if banks:
            l_eidw[:,2] = self.__add_banks(l_eidw[:,2])
        #smoothen
        neweidw = np.array([smooth(l_eidw[:,2],smoothie)]).T
        l_eidw = np.hstack((l_eidw[:,0:2],neweidw))
        #plot
        if plot:
            self.plot(l_eidw,'EIDW-LOCAL')
        return l_eidw
       
    """
        
    #############################
    ###### Save Functions #######
    #############################
        
    #Normal structure of 12 columns: inefficient (O(n^3)), but correct
    def save_as_dep(self, outputfile):
        z = deepcopy(self.z_interpol)
        nanum = -999.0
        z[np.isnan(z)] = nanum
        S = z.shape
        cols = 12
        extra = (S[1]+1)%cols
        fits = (S[1]-extra+1)
        slices = fits / cols
        newrow = np.tile(nanum, S[1])
        newcolumn = np.tile(nanum, S[0]+1)
        dep = np.vstack([z, newrow])
        dep = np.insert(dep, S[1], newcolumn, axis=1)
        #small check
        outputfile = str(outputfile)
        if not(outputfile.endswith('.dep')):
            outputfile += '.dep'
        f = open(outputfile, 'w')
        for chunk in dep:
            temp = np.split(chunk[:fits],slices)
            for part in temp:
                for value in part:
                    if value<0:
                        spacing = "  "
                    else:
                        spacing = "   "
                    f.write(spacing+'%.17E'%value)
                f.write("\n")
            if extra > 0:
                temp = chunk[-extra:]
                for value in temp:
                    if value<0:
                        spacing = "  "
                    else:
                        spacing = "   "
                    f.write(spacing+'%.17E'%value)
                f.write("\n")
        f.close()
        
    
    
    def save_as_samp(self, outputfile):
        gc.enable()
        x = np.array(self.x.flatten(),ndmin=2)
        y = np.array(self.y.flatten(),ndmin=2)
        z = np.array(self.z_interpol.flatten(),ndmin=2)
        keep = ~np.isnan(z)
        dataset = np.hstack((x.T,y.T,z.T))
        #small check
        outputfile = str(outputfile)
        if not(outputfile.endswith('.xyz')):
            outputfile += '.xyz'
        np.savetxt(outputfile, dataset[keep,:], fmt='%.17E', delimiter='\t')
        del x,y,z,keep,dataset
        
        
        
    def save_as_shape(self, outputfile):
        z = self.z_interpol.flatten()
        keep = (z!=-999.0) and (~np.isnan(z))
        z = z[keep]
        x = self.x.flatten()[keep]
        y = self.y.flatten()[keep]
        
        w = shapefile.Writer(shapefile.POINT)
        
        sx = len(str(int(max(abs(x)))))+12
        sy = len(str(int(max(abs(y)))))+12
        sz = len(str(int(max(abs(z)))))+4
        
        w.field("X","N",sx,12)
        w.field("Y","N",sy,12)
        w.field("Z","N",sz,3)
        for i in range(len(z)):
            w.point(x[i],y[i])
            w.record(x[i],y[i],round(z[i],2))
        #small check
        outputfile = str(outputfile)
        if not(outputfile.endswith('.shp')):
            outputfile += '.shp'
        w.save(outputfile)
        
        
    #############################
    ##### Private Functions #####
    #############################
    
    #TODO: Probably exclude. Interpolation does what it does, no banks described.
    def __add_banks(self, data, smoothen=True):
        #add banks ONLY if the "skipped points" have
        if sum(self.skip) > 0:
            #CROSS-SECTIONAL SLOPES:
            
            self.skip[1] = data.shape[0] - self.skip[1] #for better management
                        
            for j in range(self.cols): #for each cross-section
                #linear slopes
                data[:self.skip[0]+1,j] = np.linspace( self.wl, data[self.skip[0],j], self.skip[0]+1)
                data[self.skip[1]-1:,j] = np.linspace( data[self.skip[1]-1,j], self.wl, data.shape[0]-self.skip[1]+1)
                
                #cmn.smoothing: #???
                #if smoothen:
                    #data[:,j] = cmn.smooth(data[:,j], (data.shape[0]-sum(self.skip))/2)
                    #data[:,j] = common.smooth(data[:,j], sum(self.skip))
                
            self.skip[1] = data.shape[0] - self.skip[1] #return
        else:
            pass
        
        return data
        
        
    def __fill_nans(self,gridded_data):
        
        for i in range(self.rows):
            line = gridded_data[i,:] #one "line" of points on grid
            lnan = np.isnan(line)   #find where data is still missing on the line
            chunks = []
            flagon = False
            start = 0
            
            for j in range(self.cols):
                if lnan[j] and not(flagon):
                    flagon = True
                    start = j
                elif not(lnan[j]) and flagon:
                    flagon = False
                    chunks.append([start,j])
            if flagon:
                chunks.append([start,self.cols])
            
            for j in range(len(chunks)):
                s = chunks[j][0]
                e = chunks[j][1]
                if e == self.cols:
                    line[s-1] = line[s-2] #simple minor fix for artefacts
                    fill = [ line[s-1] ] * (e-s)
                elif s == 0:
                    fill = [ line[e] ] * e                    
                else:
                    gap = np.linspace(0.,1.,e-s+2)[1:-1]
                    fill = line[s-1]*gap[::-1] + line[e]*gap
                line[s:e] = np.array(fill)
                
            gridded_data[i,:] = line
            
        return gridded_data
        
        
    def __find_areas(self, vor, boundary):
        
        import itertools
        
        if boundary is not list:
            boundary = boundary.tolist()
        
        if boundary[0] != boundary[-1]:
            boundary.append(boundary[0])
        bounds = Polygon(boundary)
        #diagonal of bounding box = safe maximum distance to add:
        radius = common.distance( (bounds.bounds[0],bounds.bounds[1]), 
                                  (bounds.bounds[2],bounds.bounds[3]) )
        vor = self.__voronoi_correction(vor,radius)
        
        regions, vertices = self.__voronoi_finite_polygons_2d(vor,radius)
        
        areas = []
        for region in regions:
            #check for polygon validity:
            vs = np.asarray([tuple(vertices[i]) for i in region])
            if len(vs)<3:
                print "NON-POLYGON=>no area record"
                areas.append(np.nan)
                continue
            p = Polygon( vs )
            
            if not p.is_valid:
                #brute-force validation of polygon
                c = vs.mean(axis=0)            
                angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
                cc_order = np.argsort(angles)
                print len(cc_order),": INVALID POLYGON! Start permuting! (",len(angles), "points)",
                #permutation works on last elements of list first. However,
                #the first elements of the Voronoi vertices for a region tend 
                #to be the ones that are wrongly defined
                perm = itertools.permutations(cc_order[::-1])
                while True:
                    try:
                        opa = np.array(perm.next())
                        p = Polygon( vs[opa[::-1]] )
                    except:
                        print "RAN OUT OF PERMUTATIONS!"
                        raise
                    if p.is_valid:
                        print "VALIDATED POLYGON! End permuting!"
                        break
                    
            #intersect with boundary (river polygon)
            inter = p.intersection(bounds)
            #calculate area
            if inter.type not in ['Polygon','MultiPolygon']:
                print inter.type, "NON-POLYGON=>no area record"
                areas.append(np.nan) #unknown area
            elif inter.type == 'MultiPolygon':
                a = 0.0
                for poly in inter:
                    a += abs(poly.area)
                areas.append(a)
            else:
                areas.append(abs(inter.area))
                
        #adjust voronoi data structure
        vor.regions = [[]] + regions #add initial empty region
        vor.vertices = vertices
        
        return np.array(areas),vor
        
        

    #help from: http://stackoverflow.com/questions/20515554/colorize-voronoi-diagram
    def __voronoi_finite_polygons_2d(self, vor, radius=None):
        """
        Reconstruct infinite voronoi regions in a 2D diagram to finite
        regions.
    
        Parameters
        ----------
        vor : Voronoi
            Input diagram
        radius : float, optional
            Distance to 'points at infinity'.
    
        Returns
        -------
        regions : list of tuples
            Indices of vertices in each revised Voronoi regions.
        vertices : list of tuples
            Coordinates for revised Voronoi vertices. Same as coordinates
            of input vertices, with 'points at infinity' appended to the
            end.
    
        """
    
        if vor.points.shape[1] != 2:
            raise ValueError("Requires 2D input")
    
        new_regions = []
        new_vertices = vor.vertices.tolist()
    
        center = vor.points.mean(axis=0)
        if radius is None:
            radius = vor.points.ptp().max()
            
        # Construct a map containing all ridges for a given point
        all_ridges = {}
        for p1 in range(vor.points.shape[0]):
            all_ridges.setdefault(p1, []) #initialize
        for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
            all_ridges[p1].append((p2, v1, v2))
            all_ridges[p2].append((p1, v1, v2))
        
        # Reconstruct infinite regions
        for p1, region in enumerate(vor.point_region):
            
            vertices = vor.regions[region]
    
            if all(v >= 0 for v in vertices):
                # finite region
                new_regions.append(vertices)
                continue
    
            # reconstruct a non-finite region
            ridges = all_ridges[p1]
            new_region = [v for v in vertices if v >= 0]
    
            for p2, v1, v2 in ridges:
                if v2 < 0:
                    v1, v2 = v2, v1
                if v1 >= 0:
                    # finite ridge: already in the region
                    continue
                
                # Compute the missing endpoint of an infinite ridge
                t = vor.points[p2] - vor.points[p1] # tangent
                t /= np.linalg.norm(t)
                n = np.array([-t[1], t[0]])  # normal
    
                midpoint = vor.points[[p1, p2]].mean(axis=0)
                direction = np.sign(np.dot(midpoint - center, n)) * n
                far_point = vor.vertices[v2] + direction * radius
    
                new_region.append(len(new_vertices))
                new_vertices.append(far_point.tolist())
    
            # sort region counterclockwise
            vs = np.asarray([new_vertices[v] for v in new_region])
            c = vs.mean(axis=0)            
            angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
            new_region = np.array(new_region)[np.argsort(angles)]
    
            # finish
            new_regions.append(new_region.tolist())
    
        return new_regions, np.asarray(new_vertices)


    def __voronoi_correction(self, vor, max_distance):
        '''
            Helper function that corrects the Voronoi diagram.
            Any vertice that lies further than the max_distance described 
            (relative to the dataset central point) is deleted. All ridges that
            end to such vertices are assumed to extend to infinity. 
            The Voronoi data structure is corrected to accommodate these changes.
        '''
        vertices_idx_map = np.array(np.arange(vor.vertices.shape[0]))
        new_center = vor.points.mean(axis=0)
        
        for v in range(vor.vertices.shape[0]):
            if common.distance( tuple(new_center), tuple(vor.vertices[v]) ) > max_distance:
                vertices_idx_map[v] = -1      #mark for deletion
                vertices_idx_map[v+1:] -= 1   #decrease ids
        vor.vertices = vor.vertices[np.where(vertices_idx_map!=-1)] #delete
        vertices_idx_map = np.hstack((vertices_idx_map,-1)) #add last id to associate "[-1]" with last element (=-1)
        
        #fix vor.regions:
        for r in range(len(vor.regions)):
            reg = vertices_idx_map[vor.regions[r]].tolist()
            vor.regions[r] = sorted(set(reg), key=reg.index)
            
        #fix vor.ridge_vertices and vor.ridge_dict:
        for r in range(len(vor.ridge_vertices)):
            vor.ridge_vertices[r] = vertices_idx_map[vor.ridge_vertices[r]].tolist()
            vor.ridge_dict[ tuple(vor.ridge_points[r]) ] = vor.ridge_vertices[r]
            
        #delete any ridge edges that have both vertices to infinity:
        rv = np.array(vor.ridge_vertices)
        ridges_to_keep = ~np.all(rv==[-1,-1],axis=1)
        rv = rv[ridges_to_keep]
        vor.ridge_vertices = rv.tolist()
        vor.ridge_points = vor.ridge_points[ridges_to_keep] #adjust
        vor_ridge_points_list = vor.ridge_points.tolist()
        for key in vor.ridge_dict.keys():
            if list(key) not in vor_ridge_points_list:
                del vor.ridge_dict[key] #adjust
        
        return vor
