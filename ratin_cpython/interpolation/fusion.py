# -*- coding: utf-8 -*-
"""
Created on Wed May 27 10:12:26 2015

@author: zervakis
Fusion method: http://repository.tudelft.nl/view/ir/uuid%3Af41ba0a0-acd5-462e-b5aa-7c1c5878e1ac/
"""

import gc
import numpy as np
import shapefile
import matplotlib.pyplot as plt
from copy import deepcopy
from ratin.common import common


class Fusion:
    
    def __init__(self, locations, samples, model_data, interpolated_data):
        """
            Fusion class that takes cross-sections or trackline sample data,
            the result of a physics-based model (or other dataset) and an
            interpolation dataset result to calculate a unified version of the 
            two, based on the distance of the unknown points to known points.
        """
        
        #GRID LOCATIONS
        try:
            self.x = locations['x']
            self.y = locations['y']            
            self.rows, self.cols = self.x.shape
            self.s, self.n = common.to_sn(self.x,self.y)
        except KeyError:
            print "Error: Input locations should be given as a dictionary of {\'x\',\'y\'}!"
            raise
        #check
        if self.rows%2 != 1: #number of rows must be odd
            print "Error: Input locations must have an odd number of transverse coordinates!"
            raise
        
        #SAMPLES check
        if type(samples) is not type(np.array([])):
            print "ERROR: Sample data are not of matrix notation."
            raise
        elif samples.shape != self.s.shape:
            print "ERROR: Sample data do not match the provided grid locations."
            raise
        else:
            self.z = deepcopy(samples) #avoid overwrite
            self.z[self.z==-999.0] = np.nan #exchange notation of "no value"            
            self.nodata = np.isnan(self.z)  #where data is missing
            self.isdata = ~self.nodata #where sampled data exists
            
        #MODEL DATA check
        if type(model_data) is not type(np.array([])):
            print "ERROR: Model data are not of matrix notation."
            raise
        elif model_data.shape != self.s.shape:
            print "ERROR: Model data do not match the provided grid locations."
            raise
        else:
            self.model = deepcopy(model_data) #avoid overwrite
               
        #INTERPOLATION DATA
        if type(interpolated_data) is not type(np.array([])):
            print "ERROR: Interpolation data are not of matrix notation."
            raise
        elif interpolated_data.shape != self.s.shape:
            print "ERROR: Interpolation data do not match the provided grid locations."
            raise
        else:
            self.inter = deepcopy(interpolated_data) #avoid overwrite
        
        self.fused = None
    
    
    
    #Fusion function that takes cross-sections or trackline data and another's
    #interpolation result to calculate a unified version of the two, based on
    #the distance of the unknown points to known points
    def fuse(self, ftype='free-form', dthres=0.0, anisotropy=1.0, smoothie=0):
        gc.enable()
        
        print "Calculating: Fusion method ("+ftype+")"
        if ftype not in ['free-form','cross-sections']:
            print "Warning: Please specify testing dataset type as: \'free-form\' for trackline or free-form data, or \'cross-sections\' for cross-sectional data."
            print "Assuming free-form data."
            ftype = 'free-form'
        
        s, n = self.s.flatten(), self.n.flatten()
        isdata = self.isdata.flatten()
        nodata = self.nodata.flatten()
        
        #if samples in free-form (ie tracklines)
        if ftype == 'free-form':
            #check
            dthres = abs(float(dthres))
            
            #distance matrix:
            s0, n0 = s[isdata], n[isdata] #sampled
            s1, n1 = s[nodata], n[nodata] #unsampled
            dm = common.distance_matrix(s0, n0, s1, n1, anisotropy) #matrix
        
            #if no threshold is defined, calculate the furthest possible point 
            # = this is where bathymetry has the most impact
            dmins = dm.min(axis=0)
            if dthres <= 0:
                dthres = dmins.max()
            print "Threshold Distance:",dthres
    
            #weights
            w_bth = dmins / dthres #weights for model dataset
            w_bth[w_bth>1.0] = 1.0 #adjust to range of [0,1]
            w_int = 1.0 - w_bth    #weights for interpolation dataset
            del dm, dmins, s0, n0, s1, n1 #clear
        
        #if samples in cross-sections
        else:
            #find cross-sections and the points where there is data (clustering of data):
            cs_where = []
            for j in range(self.cols):
                if np.nansum(self.isdata[:,j]) != 0:
                    cs_where.append(j)
                    print "Cross-section at: ",j
            
            #Find minimum distances from 2 different cross-sections
            ratios = np.ones( (self.rows*self.cols,2) ) * np.inf
            for j in cs_where:
                dist_from_cs = common.distance_matrix( self.s[:,j][self.isdata[:,j]], self.n[:,j][self.isdata[:,j]], s, n, anisotropy )
                cs_mindist = np.min(dist_from_cs,axis=0)
                change = cs_mindist < ratios[:,1]
                ratios[:,1][change] = cs_mindist[change]
                ratios = np.sort(ratios,axis=1) #holds min distances from 2 cs
                del dist_from_cs
            
            #weights
            w_int = []
            w_bth = []
            for i in range(len(ratios)):
                #real distances
                d1 = ratios[i,0] #distance from one cs
                d2 = ratios[i,1] #distance from another cs
                if d1==0.0:
                    d1 = 1e-20 #avoid division by 0
                if d2==0.0:
                    d2 = 1e-20 #avoid division by 0
                hd = (d1+d2)/2.0 #half-distance
                if d1 < hd:
                    d2 -= hd
                    r1 = 1. / (1. + (d2/d1)**2 )
                else:
                    d1 -= hd
                    r1 = 1. / (1. + (d1/d2)**2 )
                if r1 > 1.0:
                    r1 = 1.0
                r2 = 1.0 - r1
                w_bth.append(r1)
                w_int.append(r2)
            w_int = np.array(w_int)[nodata]
            w_bth = np.array(w_bth)[nodata]
        
            
        #In both cases, fused dataset is the sum of the weighted datasets
        I = self.inter.flatten()
        M = self.model.flatten()
        fusion = self.z.flatten() #initialize, so to keep measured data values 
        fusion[nodata] = I[nodata]*w_int + M[nodata]*w_bth
        fusion = common.to_grid(fusion, self.rows, self.cols) #grid (matrix) notation
        del w_bth, w_int #clear
        
        #smoothen the longitudinals, if required
        if smoothie>0:
            for i in range(self.rows):
                fusion[i,:] = common.smooth(fusion[i,:],smoothie)
            if smoothie>int(fusion.shape[0]/4):
                smoothie = int(fusion.shape[0]/4)
            for i in range(self.cols):
                fusion[:,i] = common.smooth(fusion[:,i],smoothie)
            
        
        print "Finished [FUSION]."
        self.fused = fusion
        del fusion
        
        
    #Representation
    def plot(self, fignum, clim = 'auto', classes = 100, system = 'xy', hold = False):
        
        #append nan values
        Gz = deepcopy(self.fused)
        Gz[Gz==-999.0] = np.nan
        
        if system == 'xy':
            Gx = self.x
            Gy = self.y
        elif system == 'sn':
            Gx = self.s
            Gy = self.n
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
    
    
    #############################
    ###### Save Functions #######
    #############################
        
    #Normal structure of 12 columns: inefficient (O(n^3)), but correct
    def save_as_dep(self, outputfile):
        z = deepcopy(self.fused)
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
        z = np.array(self.fused.flatten(),ndmin=2)
        keep = ~np.isnan(z)
        dataset = np.hstack((x.T,y.T,z.T))
        #small check
        outputfile = str(outputfile)
        if not(outputfile.endswith('.xyz')):
            outputfile += '.xyz'
        np.savetxt(outputfile, dataset[keep,:], fmt='%.17E', delimiter='\t')
        del x,y,z,keep,dataset
        
        
        
    def save_as_shape(self, outputfile):
        z = self.fused.flatten()
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
