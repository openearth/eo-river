# -*- coding: utf-8 -*-
"""
.. module:: grid
    :platform: Windows
    :synopsis: 
    
.. moduleauthor:: Koen Berends <koen.berends@deltares.nl>
.. modified by::  Dimitris Zervakis


"""
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
from ratin.common import common


class Grid:
    '''
    Creates a Grid class instance to hold the positions of the node intersection
    locations, where data are.
    
    >>> grid_instant = Grid()
    '''
    
    def __init__(self):
        self.grid = {}
        self.widths = {}

        
    def build(self, network, branches=['all'], smooth=True, nlines=5, degree=5):
        '''
        Build a grid for the given network branches. 
        
        Kwargs:
            branches (list): list of branches to make a grid for. By default all
            branches
        '''
        
        try:
            nlines = int(abs(nlines))
            degree = int(abs(degree))
        except:
            print "Wrong input for Grid:build() function!"
            raise
        
        if branches == ['all']:
            branches = network.library.keys()

        for branch in branches:
            try:
                x = np.array( network.library[branch]['data']['x'] )
                y = np.array( network.library[branch]['data']['y'] )
                width = network.library[branch]['data']['width'] #"cross-sectional" width - input might be smoothened prior
                theta = network.get_direction(branch=branch, units='radians', smoothdegree=degree)  #?? here the direction is smoothened!
            except KeyError:
                print 'GridError: invalid key or network'
                return
            
            print "\nBuilding Grid for branch:",branch

            ##?? Why resmoothen a smoothened direction?? (see above ~ln.46)
#            if smooth:
#                theta = np.array(common.smooth(theta, degree))
#            else:
#                theta = np.array(theta)

            nocon = 2*nlines+1 #equal number of lines from 'top' and 'bottom' of centerline (the +1 factor)
            L = len(x)            
            xn = np.zeros((nocon, L))
            yn = np.zeros((nocon, L))

            perpendicular_theta = theta+0.5*np.pi #the new thetas to compute the new points by trigonometric functions
            perpendicular_theta[perpendicular_theta>np.pi] -= 2*np.pi  ##?can delete, because it doesn't matter (cos/sin functions follow)
            
            self.grid[branch] = {'x': [], 'y': []}
            #Construct grid
            for i in range(-1*nlines, nlines+1):
                n_widths = 0.5*i*(np.array(width)/float(nlines)) #'hypotenuses' (=width segment on n direction)
                index = i + nlines
                xn[index, :] = x + n_widths*np.array(np.cos(perpendicular_theta))
                yn[index, :] = y + n_widths*np.array(np.sin(perpendicular_theta))

            self.grid[branch]['x'] = xn
            self.grid[branch]['y'] = yn
            self.widths[branch] = np.array(width)
            
        print "Finished!"
              

        
    def plot(self, fignum=1, gridcolor='k', system='xy', block=False):
        '''
        Plot grid 
        '''
    
        fig = plt.figure(fignum)
        ax = fig.add_subplot(111)    
        
        for branch in self.grid:
    
            if system[1] == 'n':
                Salong = self._get_s_coordinates(branch)
                ly, lx = self.grid[branch]['y'].shape #rows,columns
            
            if system == 'xy':
                Gx = self.grid[branch]['x']
                Gy = self.grid[branch]['y']
            elif system == 'sn':
                Gx = np.zeros([ly,lx])
                Gy = np.zeros([ly,lx])
                for i in range(lx): #for each column
                    Gx[:,i] = Salong[i]
                    temp = 0.0
                    for j in range( (ly-1)/2 - 1, -1, -1 ):
                        s1, n1 = self.grid[branch]['x'][j+1,i], self.grid[branch]['y'][j+1,i]
                        s2, n2 = self.grid[branch]['x'][j,i], self.grid[branch]['y'][j,i]
                        temp -= common.distance([s1,n1],[s2,n2])
                        Gy[j,i] = temp
                    temp = 0.0
                    for j in range( (ly-1)/2 + 1, ly ):
                        s1, n1 = self.grid[branch]['x'][j-1,i], self.grid[branch]['y'][j-1,i]
                        s2, n2 = self.grid[branch]['x'][j,i], self.grid[branch]['y'][j,i]
                        temp += common.distance([s1,n1],[s2,n2])
                        Gy[j,i] = temp
            elif system == 'mn':
                Gx, Gy = np.meshgrid( np.arange(lx), np.arange(ly) )
            else:
                print "Wrong coordinate system chosen. Please specify input system as \'xy\', \'sn\' or \'mn\'."      
    
        
            S = self.grid[branch]['x'].shape
            for i in range(0,S[0]):
                xtemp = Gx[i,:]
                ytemp = Gy[i,:]
                ax.plot(xtemp,ytemp,color=gridcolor)
            for i in range(0,S[1]):
                xtemp = Gx[:,i]
                ytemp = Gy[:,i]
                ax.plot(xtemp,ytemp,color=gridcolor)
            
        plt.axis('equal')
        if block:
            plt.ioff()
        else:
            plt.ion()
        fig.show()
        
        
    def save(self, outputfile, branch='1'):
        '''
        Export current grid object to *.grd. This file format is used by e.g. 
        Delft3D and can be opened and viewed using `Delft3D-QUICKPLOT <https://publicwiki.deltares.nl/display/OET/OPeNDAP+access+with+Delft3D-Quickplot>`_.
        
        
        Args:
            outputfile (str): name of the output grid file
        '''
        #[branch] = self.grid.keys()
        try:
            self.grid[branch] #? unused?? why could it raise the error when the key is drawn 2 lines above? changed it and set branch as input for function.
        except KeyError:
            print 'Invalid branch'
            return

        print "\nSaving grid..."
        with open(outputfile,'w') as f:
            x = self.grid[branch]['x']
            y = self.grid[branch]['y']
            mn = x.shape
            
            f.write('* Created at '+datetime.now().strftime('%Y-%b-%d %H:%M:%S \n'))
            f.write('* by RAT-IN \n')
            f.write('Coordinate System= Cartesian \n')
            f.write(str(mn[1])+ ' ' +str(mn[0])+'\n')
            f.write('0 0 0 \n') #? what is this?
            
            # Write in 5 columns - .grd file format
            rows = int(np.ceil(mn[1]/float(5)))
            
            # X-coordinates
            for n in range(mn[0]):
                start_index = 0
                end_index   = min(5,mn[1])
                
                n_string = "%5g" % (n+1)
                
                f.write(' ETA='+n_string+'  '+' '.join(("%.17E" % i) for i in x[n,start_index:end_index])+ '\n')
                
                if not(rows < 2):
                    for row in range(rows):
                        if not(row==0):
                            start_index = row*5
                            end_index = min((row+1)*5,mn[1])
                            f.write('            '+ ' '.join(("%.17E" % i) for i in x[n,start_index:end_index])+ '\n')
              
             # Y-coordinates
            for n in range(mn[0]):
                start_index = 0
                end_index   = min(5,mn[1])
                
                n_string = "%5g" % (n+1)
                
                f.write(' ETA=' + n_string+'  '+ ' '.join(("%.17E" % i) for i in y[n,start_index:end_index])+ '\n')
                
                if not(rows < 2):
                    for row in range(rows):
                        if not(row==0):
                            start_index = row*5
                            end_index = min((row+1)*5,mn[1])
                            f.write('            '+ ' '.join(("%.17E" % i) for i in y[n,start_index:end_index])+ '\n')
        print "Success!"
    
        
        
    def save_as_shp(self, outputfile, branch='1'):
        import shapefile as sf
        try:
            self.grid[branch] #? unused?? why could it raise the error when the key is drawn 2 lines above? changed it and set branch as input for function.
        except KeyError:
            print 'Invalid branch'
            return
        print "\nSaving grid as shapefile lines..."
        x = self.grid[branch]['x']
        y = self.grid[branch]['y']
        mn = x.shape
        mnum = len(str(mn[0]))
        nnum = len(str(mn[1]))
        w = sf.Writer(sf.POLYLINE)
        w.field('ID','N',max([mnum,nnum]),0)
        w.field('DIRECTION','C','15')
        #write to shapefile
        for i in range(mn[0]):
            line = [list(l) for l in zip(x[i,:],y[i,:])]
            w.line(parts=[line])
            w.record(i,'longitudinal')
        for i in range(mn[1]):
            line = [list(l) for l in zip(x[:,i],y[:,i])]
            w.line(parts=[line])
            w.record(i,'transverse')
        w.save(outputfile)


    def load(self, gridfile):
        #Assumption of 1 single grid.
        self.grid['1'] = {}
        self.grid['1']['x'] = []
        self.grid['1']['y'] = []
        
        print "\nLoading grid..."
        f = open(gridfile, 'r')
        ETAnum = 0
        #start
        for line in f:
            words = line.split()
            if words[0] != "ETA=":
                continue
            else:
                xs = [float(x) for x in words[2:]]
                self.grid['1']['x'].append(xs)
                break
        #x-coordinates
        for line in f:
            words = line.split()
            if words[0] != "ETA=":
                xs = [float(x) for x in words]
                self.grid['1']['x'][ETAnum].extend(xs)
            elif words[1] == "1":
                ys = [float(y) for y in words[2:]]
                self.grid['1']['y'].append(ys)
                ETAnum = 0
                break
            else:
                ETAnum += 1
                xs = [float(x) for x in words[2:]]
                self.grid['1']['x'].append(xs)
        #y-coordinates
        for line in f:
            words = line.split()
            if words[0] != "ETA=":
                ys = [float(y) for y in words]
                self.grid['1']['y'][ETAnum].extend(ys)
            else:
                ETAnum += 1
                ys = [float(y) for y in words[2:]]
                self.grid['1']['y'].append(ys)
        f.close()
        
        #structure grid points as numpy arrays
        self.grid['1']['x'] = np.array(self.grid['1']['x'])
        self.grid['1']['y'] = np.array(self.grid['1']['y'])
        w = []
        for col in range(self.grid['1']['x'].shape[1]):
            chain = common.get_chainage( self.grid['1']['x'][:,col], self.grid['1']['y'][:,col] )
            w.append(chain[-1])
        self.widths['1'] = np.array(w)
                
        print "Grid loaded successfully!"

            
    # ----------------------------------------------------------------
    # Private Methods
    # ----------------------------------------------------------------
    # These methods are for internal use only       
    
    #allowed to be used outside:
    def _get_s_coordinates(self,branch):
        '''
        Returns chainage (of vertices on branch)
        '''
        
        if not(branch):
            pass
        else:
            assert type(branch) == str
            
            center = int(self.grid[branch]['x'].shape[0]/2)
            
            x = self.grid[branch]['x'][center,:]
            y = self.grid[branch]['y'][center,:]
            
            return common.get_chainage(x,y)
    
    #allowed to be used outside:
    def _smooth_widths(self,branch,degree):
        '''
        1D Gaussian smoothing of Grid width values.
        '''
        try:
            self.widths[branch] = common.smooth(self.widths[branch],degree)
        except:
            print "ERROR: Cannot smooth grid widths! [Key Error]."

    def __smooth_theta(self,alpha):    #TODO: Used anywhere? Delete??
        original_theta = np.radians(alpha)
        theta = np.copy(original_theta)
        for i in range(1,len(alpha)):
            theta1=original_theta[i-1]
            theta2 = original_theta[i]
            if theta1 < 0:
                if theta2 > 0:
                    theta1 = 2*np.pi+theta1
            else:
                if theta2 < 0:
                    theta2 = 2*np.pi+theta2
        
            theta[i] = (theta1+theta2)/2.
        return theta
