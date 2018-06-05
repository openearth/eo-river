# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 11:29:09 2015

@author: zervakis
"""

import numpy as np
#from shapely.geometry import LineString, Polygon
from scipy.interpolate import griddata
from copy import deepcopy
from ratin.common import common


class Evaluator:
    
    def __init__(self, ground_truth, prediction):
        """
        Evaluator class is used to calculate RMSE, MAE and NHWS values and
        display the results. Input must be 2 datasets of which the first is
        considered to be the ground truth dataset and the second the calculated
        prediction dataset. They both must be in numpy matrix notation 
        (array of 2 dimensions) with equal shape (rows-x-columns).
        
        NOTE: The class assumes identical locations for input datasets.
                
        Args:
             - ground_truth : <numpy.array> 
             - prediction : <numpy.array>
        
        Kwargs:
             -
        """
        #check
        try:
            if ground_truth.shape != prediction.shape:
                print "ERROR: The two datasets do not match in size!"
                raise
        except:
            print "Mismatch in input datasets"
            raise
        
        self.GT = deepcopy(ground_truth)
        self.P = deepcopy(prediction)
        self.GT = self.__noval_asnan(self.GT)
        self.P = self.__noval_asnan(self.P)
    
    
    #RMSE: Root Mean Square Error
    def rmse(self):#, computed, groundtruth):
        computed = self.P.flatten()
        groundtruth = self.GT.flatten()
        #rmse
        k = (computed - groundtruth)**2
        sub = np.sum(np.isnan(k))
        Csum = np.nansum(k)
        Dim = len(k) - sub
        if Dim == 0:
            return np.nan
        else:
            return ( Csum / Dim )**0.5


    #MAE: Mean Absolute Error
    def mae(self):#, computed, groundtruth):
        computed = self.P.flatten()
        groundtruth = self.GT.flatten()
        #mae
        err = np.abs(computed - groundtruth)
        sub = np.sum(np.isnan(err))
        Csum = np.nansum(err)
        Dim = len(err) - sub
        if Dim == 0:
            return np.nan
        else:
            return Csum / Dim


#    #TODO: recheck
#    #SLOPE: Evaluation based on area
#    def hw_slope(self, widths, plotname=False):
#        import matplotlib.pyplot as plt
#        #take half-width:
#        k = int(self.GrBathy.shape[0]/4.0)
#        TDhw = self.GrBathy[k:-k,:] #testing half-width
#        GThw = self.GrDepth[k:-k,:] #ground-truth half-width
#        numlines = TDhw.shape[0]
#        
#        #TODO: find a better way to fill this 
#        #filling of ground-truth half-widths
#        if np.sum(np.isnan(GThw[:,0])) == numlines:
#            GThw[:,0] = TDhw[:,0]
#        for j in range(1,GThw.shape[1]):
#            if np.sum(np.isnan(GThw[:,j])) > numlines-2:
#                GThw[:,j] = GThw[:,j-1]
#
#        
#        #real width
#        w = ( float(numlines-1) / (self.GrBathy.shape[0]-1) ) * np.array(widths)
#        
#        slopeDiff = []
#        a = []
#        b = []
#        gt_a = []
#        gt_b = []
#        AREA = []
#        tdArea = []
#        gtArea = []
#        for j in range(TDhw.shape[1]): #for each cross-section
#            cs_x = np.linspace(0.0,1.0,numlines) * w[j]
#            #test
#            cs_y = TDhw[:,j]
#            excl = ~np.isnan(cs_y)
#            testslope = np.polyfit(cs_x[excl],cs_y[excl],1)  #f(x) = ax+b
#            #ground-truth            
#            cs_y = GThw[:,j]
#            excl = ~np.isnan(cs_y)
#            if np.sum(excl)<=1:
#                print "This should not happen!"
#                if j == 0:
#                    linslope = testslope
#                else:
#                    linslope = [gt_a[-1],gt_b[-1]]
#            else:
#                linslope = np.polyfit(cs_x[excl],cs_y[excl],1)   #f(x) = ax+b
#            if np.isnan(linslope[0]):
#                print linslope
#                print cs_x[excl],cs_y[excl]
#            
#            a.append(testslope[0])
#            b.append(testslope[1])
#            gt_a.append(linslope[0])
#            gt_b.append(linslope[1])
#            
#            #slope difference
#            slopeDiff.append( abs(a[-1] - gt_a[-1]) )
#            
#            #areas:
#            base = w[j]
#            sideA = abs(testslope[0]*cs_x[0]+testslope[1])
#            sideB = abs(testslope[0]*cs_x[-1]+testslope[1])
#            testarea = base * (sideA+sideB) / 2.0
#            tdArea.append(base*abs(sideA-sideB)/2.0)
#            sideA = abs(linslope[0]*cs_x[0]+linslope[1])
#            sideB = abs(linslope[0]*cs_x[-1]+linslope[1])
#            linarea = base * (sideA+sideB) / 2.0
#            gtArea.append(base*abs(sideA-sideB)/2.0)
#
#            AREA.append( abs(testarea-linarea)/linarea )
#        
#        #slope's mae:
#        slopeDiff = np.array(slopeDiff)
#        sub = np.sum(np.isnan(slopeDiff))
#        Csum = np.nansum(slopeDiff)
#        Dim = len(slopeDiff) - sub
#        slope_mae = Csum / Dim
#        print "Slope's mae=",slope_mae
#        
#        area_sum = sum(AREA)
#        
#        #PLOTS:
#        if plotname:
#             #area
#            plt.close(plotname+' area')
#            plt.figure(plotname+' area')
#            plt.hold(True)
#            plt.plot(tdArea,'r.-')
#            plt.plot(gtArea,'b.-')
#            plt.hold(False)            
#            
##             #alpha
##            plt.close(plotname+' a')
##            plt.figure(plotname+' a')
##            plt.hold(True)
##            plt.plot(a,'r.-')
##            plt.plot(gt_a,'b.-')
##            plt.hold(False)
##             #beta
##            plt.close(plotname+' b')
##            plt.figure(plotname+' b')
##            plt.hold(True)
##            plt.plot(b,'r.-')
##            plt.plot(gt_b,'b.-')
##            plt.hold(False)
##             #area difference
##            plt.close(plotname+' area')
##            plt.figure(plotname+' area')
##            plt.hold(True)
##            plt.plot(AREA,'r.-')
##            plt.hold(False)
##            
##             #error(slope differences)
##            plt.close(plotname+' errorgraph')
##            plt.figure(plotname+' errorgraph')
##            plt.plot(slopeDiff,'g.-')
#        
#        return area_sum



#    #TODO: recheck
#    #SLOPE: Evaluation based on area
#    def hw_slope(self, widths, lengths, WL, plotname=False):
#        import matplotlib.pyplot as plt
#        #take half-width:
#        k = int(self.GrBathy.shape[0]/4.0)
#        TDhw = self.GrBathy[k:-k,:] #testing half-width
#        GThw = self.GrDepth[k:-k,:] #ground-truth half-width
#        numlines = TDhw.shape[0]
#        
#        #TODO: find a better way to fill this 
#        #filling of ground-truth half-widths
#        if np.sum(np.isnan(GThw[:,0])) == numlines:
#            GThw[:,0] = TDhw[:,0]
#        for j in range(1,GThw.shape[1]):
#            if np.sum(np.isnan(GThw[:,j])) > numlines-2:
#                GThw[:,j] = GThw[:,j-1]
#
#        #real width
#        w = ( float(numlines-1) / (self.GrBathy.shape[0]-1) ) * np.array(widths)
#        
#        slopeDiff = []
#        a = []
#        gt_a = []
#        AREA = []
#        tdArea = []
#        gtArea = []
#        
#        
#        for j in range(TDhw.shape[1]): #for each cross-section
#            cs_x = np.linspace(0.0,1.0,numlines) * w[j]
#            #test
#            cs_y = TDhw[:,j]
#            excl = ~np.isnan(cs_y)
#            testslope = np.polyfit(cs_x[excl],cs_y[excl],1)  #f(x) = ax+b
#            #ground-truth            
#            cs_y = GThw[:,j]
#            excl = ~np.isnan(cs_y)
#            linslope = np.polyfit(cs_x[excl],cs_y[excl],1)   #f(x) = ax+b
#            if np.isnan(linslope[0]):
#                print linslope
#                print cs_x[excl],cs_y[excl]
#            
#            a.append(testslope[0])
#            gt_a.append(linslope[0])
#            #slope difference
#            slopeDiff.append( abs(a[-1] - gt_a[-1]) )
#            
#            #signed areas:
#            base = w[j]
#            sideA = WL - testslope[0]*cs_x[0]+testslope[1]
#            sideB = WL - testslope[0]*cs_x[-1]+testslope[1]
#            trapezoidArea = base*(sideA+sideB)/2.0
#            tdArea.append(trapezoidArea)
#            #tdArea.append(np.sign(testslope[0])*trapezoidArea)
#            
#            sideA = WL - linslope[0]*cs_x[0]+linslope[1]
#            sideB = WL - linslope[0]*cs_x[-1]+linslope[1]
#            trapezoidArea = base*(sideA+sideB)/2.0
#            gtArea.append(trapezoidArea)
#            #gtArea.append(np.sign(linslope[0])*trapezoidArea)
#            
#            AREA.append( abs((gtArea[-1]-tdArea[-1])/gtArea[-1]) )
#        
#        #slope's mae:
#        slopeDiff = np.array(slopeDiff)
#        sub = np.sum(np.isnan(slopeDiff))
#        Csum = np.nansum(slopeDiff)
#        Dim = len(slopeDiff) - sub
#        slope_mae = Csum / Dim
#        print "Slope's mae=",slope_mae
#        
#        area_sum = sum(AREA)/len(AREA)
#        print area_sum
#        
#        #PLOTS:
#        if plotname:
#             #area
#            plt.close(plotname+' area')
#            plt.figure(plotname+' area')
#            plt.xlabel('River span (m)', weight='bold', size='14')
#            plt.ylabel('Enclosed Area ($m^2$)', weight='bold', size='14')
#            plt.hold(True)
#            plt.plot(lengths,tdArea,'r.-')
#            plt.plot(lengths,gtArea,'b.-')
#            plt.hold(False)
#             #alpha
#            plt.close(plotname+' a')
#            plt.figure(plotname+' a')
#            plt.xlabel('River span (m)', weight='bold', size='14')
#            plt.ylabel('Half-Width Slope (radians)', weight='bold', size='14')
#            plt.hold(True)
#            plt.plot(lengths,a,'r.-')
#            plt.plot(lengths,gt_a,'b.-')
#            plt.hold(False)
#        
#        return tdArea, gtArea


    #TODO: recheck
    #SLOPE: Evaluation based on slope only
    def hw_slope(self, widths, lengths, plotname=False):
        import matplotlib.pyplot as plt
        #take half-width:
        k = int(self.GT.shape[0]/4.0)
        CDhw = self.P[k:-k,:] #computed dataset half-width
        GThw = self.GT[k:-k,:] #ground-truth dataset half-width
        numlines = CDhw.shape[0]
        
        #TODO: find a better way to fill this 
        #filling of ground-truth half-widths // by nearest grid neighbour
        z = GThw.flatten()
        m = range(CDhw.shape[1]) * CDhw.shape[0]
        n = sorted(m)
        gpnts = np.vstack((m,n)).T
        znear = griddata(gpnts[~np.isnan(z)], z[~np.isnan(z)], gpnts[np.isnan(z)], method='nearest')
        z[np.isnan(z)] = znear
        GThw = common.to_grid(z, CDhw.shape[0], CDhw.shape[1])


        #real half-width
        w = ( float(numlines-1) / (self.P.shape[0]-1) ) * np.array(widths)
        
        slopeDiff = []
        a = []
        gt_a = []
        testNHWS = []
        gtNHWS = []
        for j in range(CDhw.shape[1]): #for each cross-section
            excl = ~np.isnan(GThw[:,j])  *  ~np.isnan(CDhw[:,j])  #mutual exclusion
            cs_x = np.linspace(0.0,1.0,numlines) * w[j]
            #test
            cs_y = CDhw[:,j]
            testslope = np.polyfit(cs_x[excl],cs_y[excl],1)  #f(x) = ax+b
            #ground-truth            
            cs_y = GThw[:,j]
            gtslope = np.polyfit(cs_x[excl],cs_y[excl],1)   #f(x) = ax+b
            if np.isnan(gtslope[0]):
                #print gtslope
                #print cs_x[excl],cs_y[excl]
                print "CHECK"

            #slope difference
            a.append(testslope[0])
            gt_a.append(gtslope[0])
            slopeDiff.append( abs(a[-1] - gt_a[-1]) )
            
            #Normalized Half-Width Slope
            H = testslope[0]*cs_x[numlines/2] + testslope[1]
            Hgt = gtslope[0]*cs_x[numlines/2] + gtslope[1]
            testNHWS.append(testslope[0]*1./H)
            gtNHWS.append(gtslope[0]*1./Hgt)
        #//end for-loop
                       
        #slope's mae (not necessary):
        print "SLOPE metrics:"
        slopeDiff = np.array(slopeDiff)
        sub = np.sum(np.isnan(slopeDiff))
        Csum = np.nansum(slopeDiff)
        Dim = len(slopeDiff) - sub
        slope_mae = Csum / Dim
        print "mae:",slope_mae
        slopeDiff = np.array(slopeDiff)**2
        sub = np.sum(np.isnan(slopeDiff))
        Csum = np.nansum(slopeDiff)
        Dim = len(slopeDiff) - sub
        slope_rmse = np.sqrt(Csum / Dim)
        print "rmse:",slope_rmse
        nsd = np.sum( np.abs(np.array(gtNHWS) - np.array(testNHWS)) ) / len(gtNHWS)
        print "nsdiff:", nsd
        
        #PLOTS:
        if plotname:
             #alpha "normalized"
            plt.close(plotname+' NHWS')
            plt.figure(plotname+' NHWS')
            plt.xlabel('River span (m)', weight='bold', size='14')
            plt.ylabel('Normalized Half-Width Slope (1/m)', weight='bold', size='14')
            plt.hold(True)
            plt.plot(lengths,testNHWS,'r.-')
            plt.plot(lengths,gtNHWS,'b.-')
            plt.hold(False)

        #return calculations; in order of output:
        #[slope of dataset to be tested, slope of ground-truth data, MAE of slope for dataset to be tested,
        # RMSE of slope for dataset to be tested, absolute slope difference]
        return [testNHWS, gtNHWS, slope_mae, slope_rmse, nsd]
        
    
    #Absolute error map
    def error_map(self, x, y, clim, plotname):
        import matplotlib.pyplot as plt
        
        err = self.GT.flatten()-self.P.flatten()
        levels = np.linspace(clim[0],clim[1],11)
        n,m = self.GT.shape
        
        plt.figure("Error Map - "+str(plotname))
        err = common.to_grid(err,n,m)
        Gx = common.to_grid(x,n,m)
        Gy = common.to_grid(y,n,m)
        
        #If more than half the points are unknown, the dataset to plot logically must be a testing dataset
        if np.sum(np.isnan(err)) > err.shape[0]*err.shape[1]/2.0:
            print "-> Assuming plot of Testing Dataset."
            cs = plt.scatter(Gx, Gy, s=10, c=err, cmap=plt.cm.bwr, vmin=clim[0], vmax=clim[1], edgecolors='none')
        else:
            cs = plt.contourf(Gx, Gy, err, levels, cmap=plt.cm.bwr, extend='both')
        
        cs.cmap.set_over('deeppink')
        cs.cmap.set_under('k')
        cs.set_clim(np.floor(levels[0]), np.ceil(levels[-1]))
        cbar = plt.colorbar(cs)
        cbar.set_label('Bed Level Difference (m)', labelpad=15, weight='bold', size='14')
        plt.axis('equal')
        plt.ion()
        plt.show()


    #Change -999.0 values to NaN
    def __noval_asnan(self, data):
        data[data == -999.0] = np.nan
        return data
