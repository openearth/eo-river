# -*- coding: utf-8 -*-
"""
.. module:: bathymetry
    :platform: Windows
    :synopsis: This module provides the bathymetry based on the Width,
               Curvature, Discharge and Slope

.. module author:: Ymkje Huismans <ymkje.huismans@deltares.nl>
.. date:: 5th of May 2014

The optimized route module searches for the route which requiers least dredging
given the rules for navigation

"""

# to do:
# calculate volume
# calculate two cells if navwidth > cellwidth
# optimization of the route
# curvature
# write output to route

from __future__ import division
import numpy as np
import math as math
import matplotlib.pyplot as plt

class route():

    def __init__(self, grid, bath):
    
        self.gr = grid
        self.bath = bath
        self.route = {}

    def calc_route(self, branches='Bala', dep='h loaded', shipwidth = 5, draught = 8, jump = 20):

        # rules for navigation        
        navwidth = shipwidth + 15
        navdep = draught + 5 
        dist_to_bank = 20 + shipwidth/2
        
        # initialize and load values
        h = self.bath.bath[branches][dep]
        x = self.gr.grid[branches]['x']
        y = self.gr.grid[branches]['y']
        s = h.shape
        
        # exclude cells too close to the bank
        dis_n = np.zeros(x.shape[1])
        n = np.zeros(x.shape[1])
        for i in range(x.shape[1]):
            distemp = self._distance(x[0:2,i], y[0:2,i])
            n[i] = np.ceil(dist_to_bank/distemp[0])    
            dis_n[i] = distemp
            excl = n[i]+1    
            h[0:excl,i] = -10
            h[s[0]-excl:s[0],i] = -10

        print '\nSmallest cell size in n direction:'
        print min(dis_n)

        print '\nAverage cell size in n direction:'
        print np.mean(np.array(dis_n))
        
        # identify areas where no dredging is requiered
        depminlim = h-navdep
        navigable = depminlim.clip(0)
        for i in range(s[0]):
            for k in range(s[1]):
                if navigable[i,k]!=0:
                    navigable[i,k] = 1
        
        inav = np.max(navigable, 0)

        # first approx = deepest points, 1. boundary conditions: distance to bank
        # --> think of clever restriction mode!!!
        dredge = np.zeros([s[1]])
        deproute = np.zeros([s[1]])
        depcheck = np.zeros([s[1]])
        xroute = np.zeros([s[1]])
        yroute = np.zeros([s[1]])
        index = np.zeros([s[1]])
        for i in range(0,s[1]):
            if inav[i]==2:  # --> has been switched off!!!!
                index[i] = 'NaN'
                xroute[i] = 'NaN'
                yroute[i] = 'NaN'
                dredge[i] = 0
            else:
                index[i] = h[:,i].argmax()
                deproute = h[index[i],i]
                xroute[i] = x[index[i],i]
                yroute[i] = y[index[i],i]
                dredge[i] = navdep - deproute
        
        d = dredge.clip(0)
        print d.sum()
        
        # remove zigzags
        
        
        
        
        
        # all plots: (remove later)
        
#        plt.figure('no dredging')
#        cs = plt.contourf(x, y,navigable,levels = [0,1], cmap=plt.cm.Blues, extend = 'min')
#        cs.cmap.set_under('k')
#        plt.colorbar(cs)
#        plt.axis('equal')            
#        plt.ion()
#        plt.show()
        
        plt.figure('index')
        plt.plot(index, '-b.')
        plt.grid()
        
        self.bath.plot(fignum = 'bath', which = dep)
        plt.figure('bath')
        plt.plot(xroute, yroute, 'r', linewidth=2)
        plt.scatter(xroute, yroute, d*100, d, cmap=plt.cm.Reds)
        plt.colorbar()
#        
#        plt.figure('route no route')
#        plt.plot(xroute, yroute, 'k', linewidth=3)
#        xroute[20:120] = 'NaN'
#        plt.plot(xroute, yroute, 'r', linewidth=3)
        

    def _distance(self,xn,yn):
        L = len(xn)-1    
        dis = np.zeros(L)
        ind = 0
        for i in range(len(xn)-1):
            ind = ind + 1
            dx = xn[i]-xn[i+1]
            dy = yn[i]-yn[i+1]
            dis[i] = math.sqrt(dx**2+dy**2)
        
        return dis