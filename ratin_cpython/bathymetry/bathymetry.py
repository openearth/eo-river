 # -*- coding: utf-8 -*-

"""
.. module:: bathymetry
    :platform: Windows
    :synopsis: This module provides the bathymetry based on the Width,
               Curvature, Discharge and Slope

.. module author:: Ymkje Huismans <ymkje.huismans@deltares.nl>
.. modified by::   Dimitris Zervakis
.. date:: 9th of April 2014
.. mod. date:: 27th of January 2016

The bathymetry module provides the bed level in case no bed level measurements
are present. In a later stage this module will be extended by generating the
bathymetry based on single beam measurements or other data that have no full
area coverage.

The current module calculates the bathymetry in two steps:

1. 1D bed level, based on the formula of Chezy (zeroth order bed level)
2. Adding 2D features with the first order pertubation (axi-symmetric solution)

The theory described below is based on Crosato [1]_.

 .. [1] Analysis and modelling of river meandering, Thesis A. Crosato, TUDelft 2008


General theory
===============================

1. Zeroth order
-------------------------------

The zero-order solution to the momentum and continuity equations provides the
bed topography for a uniform flow in a straight and infinitely long channel:

.. math::
    h_{c} = (Q/(BC\sqrt i))^{2/3}

with \n

| h :sub:`c` - water depth at the centerline [m]
| Q - the river discharge [m :sup:`3` /s]
| B - the width [m]
| C - the Chezy roughness [m :sup:`1/2` /s]
| i - the slope [-]
|

2. Axi-Symmetric solution
-------------------------------

A first estimate for the effect of river bends on the bed level can be
estimated with the axi-symmetric solution of the equations, assuming mildly
curved channels:

.. math::
    h = h_{c}e^{Af(\\theta)n/R_{c}}

with \n

| h     -   water depth along n [m]
| hc    -   water depth at the centerline [m], calculated above
| A     -   coefficient (eqn 5.10), see expression below
| :math:`f(\\theta)`, weighing function (eqn 5.6), see expression below
| n     -   coordinate orthogonal to the streamline [m]
| R :sub:`c`    -   Radius of curvature [m]
|

In which:

.. math::
    A = 2\\alpha/\\kappa^2(1-\sqrt g/\\kappa C)

with

| :math:`\\alpha` - calibration coefficient
| :math:`\\kappa` - von Karman constant, usually 0.44
| g - gravitational acceleration [m/s2]
|

In which:

.. math::
    f(\\theta) = 0.85/E\sqrt(\\theta)

with

| :math:`\\theta` - shields parameter [-], see expression below
| E - calibration coefficient, see expression below
|

In which:

.. math::
    \\theta = (u^2+v^2)/(C^2\\Delta D_{50})

with

| u, v - velocity in resp. streamwise and transverse direction [m/s]
| :math:`\\Delta = (\\rho_{s} - \\rho_{w})/\\rho_{w}` [-], with :math:`\\rho_{s} = 2650 kg/m^{3}, \\rho_{w} = 1000 kg/m^{3}`
| :math:`D_{50}` - mean grainsize diameter
|

In which:

.. math::
    E = 0.0944(h/D_{50})^{0.3}
"""
#Internal notes: to be improved or tested
#-------------------------------
#
#* try different grid sizes (11, 13, 15 etc) to check robustness
#* different ways to specify branch & possibility to run through all braches
#* smooth function which causes a grid error
#* which discharge is representative for the bed level
#* finer array
#* side slopes from the banks, so that the deeper parts aren't always on the sides

import numpy as np
import os
import shapefile
import matplotlib.pyplot as plt
from ratin.common import common
from shapely.geometry import Polygon,LineString

class Bathymetry:

    def __init__(self, grid, banks = None):
        #Initialize
        self.gr = grid
        self.bath = {}

        # These values are set by default
        self.rhow = 1000.0        # density of water      [kg/m3]
        #TODO: Shouldn't this be provided as input?
        self.rhos = 2025.0        # density of sediment   [kg/m3]
        
        if type(banks) is list:
            #if banks are consistently defined as a certain width:
            try:
                self.banks = [abs(banks[0]),abs(banks[1])]
            except:
                print "Warning: bank profiles not described correctly. Assuming NO profile."
                self.banks = None
        else:
            self.banks = None
                
        

#TODO: #? Shouldn't Chezy be calculated in some way (chezy/manning)   ?
    def create_dep(self, branches='all', 
                   Q = 500, C = 45, i_slope = 9e-5, D50 = 300e-6, a = 10, WL = 0, 
                   deg = 20, system = 'm', banks_profile ='linear', curvature = 'global',
                   fix_CbyH = False, fix_CbyB = False, fix_abyB = False,
                   lonsmooth = False, transmooth = True, 
                   confined_width = True, adjust_byArea = False, remove_shallow_banks = False):
        
        if branches=='all':
            branches = self.gr.grid.keys()
        elif type(branches) != type([]):
            print "Please provide names of branches in a list.\n"
            return

        for branch in branches:

            print "Calculating bed levels for branch {0}".format(branch)
            
            try:
                #In order to extract the local curvature, it must be calculated on each point of the grid
                self.bath[branch] = {}
                Gx = self.gr.grid[branch]['x']
                Gy = self.gr.grid[branch]['y']
                k = self.__get_curvature(Gx,Gy,deg,which=curvature)
                widths = self.gr.widths[branch]
                
                #TODO: Experimenting with further smoothing: S-Golay smoothing on curvature (per longitudinal):
                if deg>=3:
                    if deg%2 != 1:
                        smoothfactor = int(deg - deg%2 + 1)
                    else:
                        smoothfactor = int(deg)
                    for q in range(k.shape[0]):
                        k[q,:] = common.savitzky_golay(k[q,:], smoothfactor, 1)
                
                
                #TODO: delete!
                s = np.array(self.gr._get_s_coordinates('1'))
                w = np.array(widths)
                if system=='ft':
                    s*=0.3048
                    w*=0.3048
                curv = k[k.shape[0]/2,:]
                """
                #TODO: LAG????
                curv_s = s + w/2.0 #ASSUMPTION!!!!
                curv_ls = LineString([(0,curv[0])]+zip(curv_s,curv))
                curv_new = []
                for i in s:
                    curv_new.append(curv_ls.interpolate(i).coords[:][0][1])
                curv_new = np.array(curv_new)
                for i in range(k.shape[0]):
                    k[i,:] = curv_new
                """
                plt.figure('WIDTH VS CURVATURE')
                plt.subplot(2,1,1)
                plt.plot(s,w)
                plt.ylabel('widths (m)')
                plt.subplot(2,1,2)
                plt.plot(s,curv)
                plt.ylabel('curvature [-]')
                plt.xlabel('Longitudinal distance (m)')
                plt.ion()
                plt.show()
                
                norow, nocon = self.gr.grid[branch]['x'].shape                
            except KeyError:
                print 'invalid key or parameter value'
                return
            
            #If the data are in feet, change widths into meters for equations to work
            if system=='ft':
                B = np.array(widths) * 0.3048
                if self.banks:
                    self.banks = [j*0.3048 for j in self.banks]
            else:
                B = np.array(widths)
            
            #TODO: DELETE!!
            if fix_CbyB:
                C *= B/np.average(B)
            #TODO: ???
            #k *= np.average(B)/B
            #a += (np.average(B)/B)
            if fix_abyB:
                a *= np.average(B)/B
            
            ### FIND BANKS AND ORIGINAL AREAS ###
            hc_ref = ( float(Q) / ( B * C * i_slope**0.5 ) ) ** (2.0/3.0)
            if fix_CbyH:
                href = np.average(hc_ref)
                p = 1./6
                hc_original = hc_ref**(3./(2*p+3.)) * href**(2*p/(2*p+3.))
                Cp = (hc_original/href)**p * C #fix Chezy
                #hc = hc_ref**(3./(3.-2*p)) * href**(-2*p/(3.-2*p))
                #Cp = (href/hc)**p * C #fix Chezy
                #TODO: delete
                plt.figure('WIDTH VS Chezy')
                plt.subplot(2,1,1)
                plt.plot(s,w)
                plt.ylabel('widths (m)')
                plt.subplot(2,1,2)
                plt.plot(s,Cp)
                plt.ylabel('Chezy coeff')
                plt.xlabel('Longitudinal distance (m)')
                plt.ion()
                plt.show()
            else:
                hc_original = hc_ref
                Cp = C * np.ones(hc_original.shape)
            E = 0.0944 * ((hc_original/D50)**0.3)
            delta = (self.rhos-self.rhow)/self.rhow
            theta = 1.0 / (Cp**2 * delta * D50)
            ft = (0.85/E) * (theta**0.5)
            Rcinv = -k
            karman = 0.44
            g = 9.80665
            A = (2.0*a)/(karman**2) * (1.0 - np.sqrt(g)/(karman*Cp) )
            hf = np.ones((norow, nocon)) * np.nan #initialize
            bankslope_width_r = []
            bankslope_width_l = []
            original_areas = []
            #for each cross-section:
            for j in range(nocon):
                n = np.linspace(-B[j]/2.0,B[j]/2.0,norow)
                power = A[j] * ft[j] * n * Rcinv[:,j]
                hf[:,j] = hc_original[j] * np.exp(power)
                
                #areas of original theory calculations
                xy = [(-B[j]/2.0,0.0)] + zip(n,hf[:,j]) + [(B[j]/2.0,0.0)] #close polygon
                p = Polygon(xy)
                original_areas.append(p.area)
                
                #TODO: DELETE!
                if j==300:
                    store1 = np.array(hf[:,j][::])
                #TODO: DELETE!
                if j==350:
                    store2 = np.array(hf[:,j][::])
                
                #compute banks widths:
                if self.banks:
                    if self.banks[0] > B[j]/3.0:
                        #TODO: delete next line's print!
                        print j, self.banks[0], B[j], self.banks[0]/B[j]
                        bankslope_width_r.append(B[j]/3.0)
                    else:
                        bankslope_width_r.append(self.banks[0])
                    if self.banks[1] > B[j]/3.0:
                        #TODO: delete next line's print!
                        print j, self.banks[1], B[j], self.banks[1]/B[j]
                        bankslope_width_l.append(B[j]/3.0)
                    else:
                        bankslope_width_l.append(self.banks[1])
                else:
                    add = max(hf[:,j]) #find maximum depth calculated
                    #check that banks do not exceed computable size
                    if add > B[j]/3.0:
                        #TODO: delete print!
                        print j, add, B[j], add/B[j]
                        add = B[j]/3.0
                    bankslope_width_r.append(add)
                    bankslope_width_l.append(add)
            bankslope_width_r = np.array(bankslope_width_r)
            bankslope_width_l = np.array(bankslope_width_l)
            
            #If chosen, the banks on the shallow side of the river bend can be excluded
            if remove_shallow_banks:
                sign = np.sign(Rcinv[int(norow/2),:])
                bankslope_width_r[sign>0] = 0.0
                bankslope_width_l[sign<0] = 0.0
                
            """
            #TODO: delete!
            shallowest_points = []
            for j in range(nocon):
                shallowest_points.append(min(hf[:,j]))
            plt.figure('Shallows')
            plt.plot(s,shallowest_points)
            plt.ion()
            plt.show()
            """
            
            ### ZERO-ORDER ####################################################
            if confined_width:
                B_model = B - (bankslope_width_r + bankslope_width_l)
            else:
                B_model = B
            hc_ref = ( float(Q) / ( B_model * C * i_slope**0.5 ) ) ** (2.0/3.0)
            if fix_CbyH:
                href = np.average(hc_ref)
                p = 1./6
                hc = hc_ref**(3./(2*p+3.)) * href**(2*p/(2*p+3.))
                Cp = (hc/href)**p * C #fix Chezy
                #hc = hc_ref**(3./(3.-2*p)) * href**(-2*p/(3.-2*p))
                #Cp = (href/hc)**p * C #fix Chezy
            else:
                hc = hc_ref
                Cp = C * np.ones(hc.shape)
            ###################################################################
            
            ### AXI-SYMMETRIC: ################################################
            E = 0.0944 * ((hc/D50)**0.3)
            delta = (self.rhos-self.rhow)/self.rhow
            theta = 1.0 / (Cp**2 * delta * D50)
            ft = (0.85/E) * (theta**0.5)
            Rcinv = -k
            karman = 0.44
            g = 9.80665
            A = 2.0*a/karman**2 * (1.0 - np.sqrt(g)/(karman*Cp) )
            hf = np.ones((norow, nocon)) * np.nan #initialize
            new_areas = []
            for j in range(nocon): #for each cross-section
                n = np.linspace(-B[j]/2.0,B[j]/2.0,norow) #calculate for whole area; excess will be "cut down" by the bank profiles
                power = A[j] * ft[j] * n * Rcinv[:,j]
                hf[:,j] = hc[j] * np.exp(power)
                #if there are banks described:
                discrete_chunks = B[j] / (norow-1.0)
                #right bank:
                if (bankslope_width_r[j]>0.0):
                    mr = int(np.round(bankslope_width_r[j]/discrete_chunks))
                    if mr == 0:
                        mr = 1 #at least 1
                    if banks_profile == 'linear':
                        sloppy_r = np.linspace( 0.0, hf[mr,j], mr+1 ) #right
                    elif banks_profile == 'hyper':
                        spread = np.linspace( -np.sinh(1.0), np.sinh(1.0), mr+1 )
                        sloppy_r = (np.arcsinh(spread)+1.0) * hf[mr,j]/2.0 #right
                    hf[:mr+1,j] = sloppy_r #attach right bank
                    bankslope_width_r[j] = mr*discrete_chunks
                #left bank:
                if (bankslope_width_l[j]>0.0):
                    ml = int(np.round(bankslope_width_l[j]/discrete_chunks))
                    if ml == 0:
                        ml = 1 #at least 1
                    if banks_profile == 'linear':
                        sloppy_l = np.linspace( 0.0, hf[norow-ml-1,j], ml+1 ) #left
                    elif banks_profile == 'hyper':
                        spread = np.linspace( -np.sinh(1.0), np.sinh(1.0), ml+1)
                        sloppy_l = (np.arcsinh(spread)+1.0) * hf[norow-ml-1,j]/2.0 #left
                    hf[norow-ml-1:,j] = sloppy_l[::-1] #attach left bank
                    bankslope_width_l[j] = ml*discrete_chunks
                    
                #areas
                xy = zip(n,hf[:,j])
                if hf[-1,j] != 0:
                    xy = xy + [(B[j]/2.0,0.0)]
                if hf[0,j] != 0:
                    xy = [(-B[j]/2.0,0.0)] + xy
                new_areas.append(Polygon(xy).area)
                
                #if an adjustment to initially computed areas is chosen:
                if adjust_byArea:
                    Adiff = original_areas[j] - new_areas[j]
                    h_plus = Adiff / ( B[j] - (bankslope_width_l[j]+bankslope_width_r[j])/2.0 )
                    hf[:,j] += h_plus
                    hc[j] += h_plus
                    #FIX banks:
                    if (bankslope_width_r[j]>0.0):
                        if banks_profile == 'linear':
                            sloppy_r = np.linspace( 0.0, hf[mr,j], mr+1 ) #right
                        elif banks_profile == 'hyper':
                            spread = np.linspace( -np.sinh(1.0), np.sinh(1.0), mr+1 )
                            sloppy_r = (np.arcsinh(spread)+1.0) * hf[mr,j]/2.0 #right
                        hf[:mr+1,j] = sloppy_r #attach right bank
                    #left bank:
                    if (bankslope_width_l[j]>0.0):
                        if banks_profile == 'linear':
                            sloppy_l = np.linspace( 0.0, hf[norow-ml-1,j], ml+1 ) #left
                        elif banks_profile == 'hyper':
                            spread = np.linspace( -np.sinh(1.0), np.sinh(1.0), ml+1)
                            sloppy_l = (np.arcsinh(spread)+1.0) * hf[norow-ml-1,j]/2.0 #left
                        hf[norow-ml-1:,j] = sloppy_l[::-1] #attach left bank
                    #AREAS
                    #recalculate final areas for comparison purposes
                    xy = zip(n,hf[:,j])
                    if hf[-1,j] != 0:
                        xy = xy + [(B[j]/2.0,0.0)]
                    if hf[0,j] != 0:
                        xy = [(-B[j]/2.0,0.0)] + xy
                    new_areas[j] = Polygon(xy).area
                        
                        
            #store finalized 0-order:
            if system == 'ft': #if specified system is in feet
                self.bath[branch]['hc'] = hc / 0.3048 #return notation back to feet
            else:
                self.bath[branch]['hc'] = hc                        
                        
            #TODO: delete graph
            plt.figure("MODEL_AREAS")
            plt.hold(True)
            plt.plot(s,original_areas)
            plt.plot(s,new_areas)
            plt.hold(False)
            plt.ion()
            plt.show()
            
            
                
            

##TODO: Inspect: 2D EXPERIMENTAL SMOOTHING:
#            hf = common.smooth2D(hf, 0.0, deg/2, 3.0, 1.0, 5.0) #OR:
#            hf[bankslope_width:norow-bankslope_width+1,:] = \
#                     common.smooth2D(hf[bankslope_width:norow-bankslope_width+1,:], 0.0, deg/2, 5.0, 1.0, 5.0)
#
            ### SMOOTHING ###
            #transverse
            if transmooth:
                for j in range(nocon): #smoothen each cross-section
                    hf[:,j] = common.smooth(hf[:,j], int(0.2*norow)) #Smoothing window: 20% of cross-section width
                #TODO: Should the side banks be returned to expected 0-water level after smoothing?
                #hf[0,:] = 0.0
                #hf[-1,:] = 0.0
            #longitudinal
            if lonsmooth and deg>0:
                for i in range(norow):  #smoothen each row on grid
                    hf[i,:] = common.smooth(hf[i,:], deg)
                
#TODO: delete - creates a few cross-sectional profiles for representation reasons ## all based on water level = 0
            plt.figure("a cross-section")#, figsize=(10,10))
            plt.hold(True)
            for j in range(hf.shape[1]): #smoothen each cross-section
                if j==300: #in range(0,1500,120):
                    bbb = np.linspace(-B[j]/2.0, B[j]/2.0, norow)
                    plt.plot(bbb,-hf[:,j],'.-',label=str(j))
                    plt.plot(bbb,-store1,'.-',label=str(j))
                if j==350: #in range(0,1500,120):
                    bbb = np.linspace(-B[j]/2.0, B[j]/2.0, norow)
                    plt.plot(bbb,-hf[:,j],'.-',label=str(j))
                    plt.plot(bbb,-store2,'.-',label=str(j))
            plt.hold(False)
            axes = plt.gca()
            axes.set_ylim([-60,40])
            plt.axis('equal')
            plt.show()

            if system == 'ft': #if specified system is in feet
                hf /= 0.3048 #return notation back to feet
                if self.banks:
                    self.banks = [j/0.3048 for j in self.banks]
            self.bath[branch]['h'] = hf #water levels
            
            Salong = self.gr._get_s_coordinates(branch)
            #Adjust to bed levels (assuming WL = Water Level at start of river section, following the water flow)
            WLdecrease = ( Salong * i_slope ) / ( (i_slope**2 + 1.0)**0.5 )
            self.bath[branch]['bedlvl'] = (WL-WLdecrease) - hf  #returns feet if feet, or meters if meters
            
            print "Model bed levels calculated for branch: {0}!".format(branch)



    def plot(self, fignum = 'dep', branch = '1', clim = 'auto', which = 'bedlvl', classes = 100, system = 'xy'):
        
        allowed_plots = ['h', 'hc', 'bedlvl', 'h_load']
        if which not in allowed_plots:
            print "Please specify plot ('h', 'hc', 'bedlvl' or 'h_load')."
            return
        else:
            hf = self.bath[branch][which]
        
        if which == 'hc':
            hf = np.repeat( np.array(hf,ndmin = 2), self.bath[branch]['h'].shape[0], axis=0 )
        
        #append nan values
        nanvalues = zip(*np.where(hf==-999.)) #assuming -999 is the "no value"
        for n in nanvalues:
            hf[n] = np.nan
            
        if system[1] == 'n':
            Salong = self.gr._get_s_coordinates(branch)
            ly, lx = self.gr.grid[branch]['y'].shape #rows,columns
        
        if system == 'xy':
            Gx = self.gr.grid[branch]['x']
            Gy = self.gr.grid[branch]['y']
        elif system == 'sn':
            Gx = np.zeros([ly,lx])
            Gy = np.zeros([ly,lx])
            for i in range(lx): #for each column
                Gx[:,i] = Salong[i]
                temp = 0.0
                for j in range( (ly-1)/2 - 1, -1, -1 ):
                    s1, n1 = self.gr.grid[branch]['x'][j+1,i], self.gr.grid[branch]['y'][j+1,i]
                    s2, n2 = self.gr.grid[branch]['x'][j,i], self.gr.grid[branch]['y'][j,i]
                    temp -= common.distance([s1,n1],[s2,n2])
                    Gy[j,i] = temp
                temp = 0.0
                for j in range( (ly-1)/2 + 1, ly ):
                    s1, n1 = self.gr.grid[branch]['x'][j-1,i], self.gr.grid[branch]['y'][j-1,i]
                    s2, n2 = self.gr.grid[branch]['x'][j,i], self.gr.grid[branch]['y'][j,i]
                    temp += common.distance([s1,n1],[s2,n2])
                    Gy[j,i] = temp
        elif system == 'mn':
            Gx, Gy = np.meshgrid( np.arange(lx), np.arange(ly) )
        else:
            print "Wrong coordinate system chosen. Please specify input system as \'xy\', \'sn\' or \'mn\'."        

        #plt.close(fignum)
        plt.figure(fignum)
        #plt.clf()
        
        # set colorlimits
        minn= np.nanmin(hf)
        maxx= np.nanmax(hf)
#       #initial 'auto': did not perform well with varying extents; changed
        if clim == 'auto':
            mean = np.mean(np.ma.masked_array(hf,np.isnan(hf)))
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
        
        cs = plt.contourf(Gx, Gy, hf, levels, cmap=plt.cm.jet, extend="both")
        cs.cmap.set_under('k')
        cs.set_clim(np.floor(levels[0]), np.ceil(levels[-1]))
        plt.colorbar(cs)
        plt.axis('equal')
        plt.ion()
        plt.show()
        
        #fix back:
        for n in nanvalues:
            hf[n] = -999.


    def load_dep(self, inputfile, branch = '1'):
        
        import warnings
        
        inputfile = inputfile.encode('string-escape')
        inputfile = inputfile.replace('\\','/')
        inputfile = inputfile.replace('//','/')
        print "Trying to load: "+str(inputfile)
        
        # Extract file extension
        extension = os.path.splitext(inputfile)[-1]
        
        if extension == '.dep':
            #input must be of same extents with the grid
            S = self.gr.grid[branch]['x'].shape
            cols = 12
            lines_per_chunk = int((S[1]+1)/cols)
            extra = (S[1]+1)%cols
            if extra > 0:
                lines_per_chunk += 1
            #open file
            with open(inputfile) as f:
                #1st chunk:
                data = []
                for l in range(lines_per_chunk):
                    line = f.readline()
                    temp = line.split()
                    data.extend([float(d) for d in temp])
                A = np.array(data)
                #for each other chunk:
                for i in range(1,S[0]): #excludes last chunk of -999.0's values
                    data = []
                    for l in range(lines_per_chunk):
                        line = f.readline()
                        temp = line.split()
                        data.extend([float(d) for d in temp])
                    A = np.vstack((A,data))
                A = A[:,:-1] #exclude last row and last column of -999.0's
            f.close()
            
        elif extension == '.xyz':
            warnings.warn("Error: Cannot load sample file. Must be a depth file ('.dep')")
        else:
            warnings.warn("Error: Input file is not of an implemented extension!")
        
        try:
            self.bath[branch]['h_load'] = A
        except:
            self.bath[branch] = {}
            self.bath[branch]['h_load'] = A
        
        print "Load successful.\n"
        
        

    def save_as_samp(self, outputfile, branch = '1'):

        hf = self.bath[branch]['bedlvl']
        x = self.gr.grid[branch]['x']
        y = self.gr.grid[branch]['y']
        
        S = hf.shape
        xyz = np.zeros((S[0]*S[1], 3))
        xyz[:,0] = np.reshape(x, S[0]*S[1])
        xyz[:,1] = np.reshape(y, S[0]*S[1])
        xyz[:,2] = np.reshape(hf, S[0]*S[1])
        #small check
        outputfile = str(outputfile)
        if not(outputfile.endswith('.xyz')):
            outputfile += '.xyz'
        np.savetxt(outputfile, xyz, fmt='%.17E', delimiter='\t')


    #TODO: More graceful code-writing?
    #Normal structure of 12 columns:  #?? very inefficient (O(n^3)), but correct
    def save_as_dep(self, outputfile, branches = '1'):

        hf = self.bath[branches]['bedlvl']
        S = hf.shape
        nanum = -999.0
        cols = 12
        extra = (S[1]+1)%cols
        fits = (S[1]-extra+1)
        slices = fits / cols
        newrow = np.tile(nanum, S[1])
        newcolumn = np.tile(nanum, S[0]+1)
        dep = np.vstack([hf, newrow])
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
        
        
    def save_as_shape(self, outputfile, branches = '1'):
        hf = self.bath[branches]['bedlvl'].flatten()
        keep = (hf!=-999.0)
        z = hf[keep]
        x = self.gr.grid[branches]['x'].flatten()[keep]
        y = self.gr.grid[branches]['y'].flatten()[keep]
        
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
        
        
    #Possibly best way... (instead of even averages takes odd ones, with central value the current id)
    def _runningMean(self, x, N):
        N = int(N)
        if N<=0:
            return x
        else:
            y = np.convolve(x, np.ones(2*N+1,)/(2*N+1), 'same')
            for i in range(N):
                y[i] = np.sum(x[:2*i+1])/(2*i+1)
                y[-i-1] = np.sum(x[-2*i-1:])/(2*i+1)
            return y
            
    #TODO: Fix. Ambiguous idea: have a forward running mean. N forward values and N/2 back are averaged
    def _fwdrunningMean(self, x, N):
        N = int(N)
        if N<=0:
            return x
        else:
            w = 1./(np.arange(1,N+2))
            b = int(np.ceil(N/2.0))
            w = np.hstack((w[1:b+1][::-1],w))
            w = w/np.sum(w)
            start = x[:b]
            y = []
            end = x[-N:]
            for i in range(b,len(x)-N):
                y.append(np.sum(x[i-b:i+N+1]*w))
            y = np.array(y)
            y = np.hstack((start,y))
            y = np.hstack((y,end))
            return y


    #Local cell edge curvature (ds/dtheta).
    #dtheta can be equal on the same transverse, but ds varies
    def __get_curvature(self, Gx, Gy, deg, which='local'):
        
        k = np.ones(Gx.shape)
        
        if which == 'local':
            for i in range(Gx.shape[0]): #for each row of the grid
                x = Gx[i,:]
                y = Gy[i,:]
                #calculate local chainage (for the currect grid row)
                s = common.get_chainage(x,y)
                #calculate local direction
                theta = common.get_direction(x,y, smoothdegree=deg, units='radians')
                #curvature calculation
                curvature = np.diff(theta)/np.diff(s)
                #starting curvature set to 0.0
                k[i,:]  = np.append(curvature[0],curvature)
        elif which == 'global':
            mid = int(Gx.shape[0]/2)
            x = Gx[mid,:]
            y = Gy[mid,:]
            #calculate global chainage
            s = common.get_chainage(x,y)
            #calculate local direction
            theta = common.get_direction(x,y, smoothdegree=deg, units='radians')
            #curvature calculation
            curvature = np.diff(theta)/np.diff(s)
            #starting curvature set to 0.0
            k  = k * np.append(curvature[0],curvature)
            
        return k
        