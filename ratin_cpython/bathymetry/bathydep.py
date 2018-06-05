# -*- coding: utf-8 -*-
"""
.. module:: bathymetry
    :platform: Windows
    :synopsis: This module provides the bathymetry based on the Width,
               Curvature, Discharge and Slope

.. module author:: Dimitris Zervakis

The bathymetry depths module aggregates a set of data points onto a specified
grid. The function at the moment automatically defines a threshold radius per
grid cell and snaps the closest data point onto the grid's edge from which the
radius is drawn. If no point lies within the radius, a NaN value is set.

"""

import os, copy, shapefile, itertools
import numpy as np
from scipy.spatial import KDTree
import matplotlib.pyplot as plt
from ratin.common import common


class BathyDep:

    def __init__(self, grid, samples, zfield=2):

        self.gr = grid
        self.griddeddata = {} #dictionary of data on grid (after averaging or loading)
        self.test = {}  #optional dictionary to hold testing dataset values
        self.samples = [] #sampled point list of (X,Y,Z) data
        
        if type(samples) is str: #expect filename
            self.__read_file(samples,zfield)
        elif type(samples) is list: #expect list of triplets of x,y,z
            temp = np.array(samples)
            if temp.shape[1]!=3:
                print "ERROR: No samples loaded - expected triplets of data samples! (X,Y,Z)"
            else:
                self.samples = samples
        else:
            print "ERROR: No samples loaded - incorrect samples type!"
    
    
    def avgdep(self, branches='all', intensity = 1.0):
        """
            Aggregates the data on the specified grid cell edges locations.
        """
        
        intensity = abs(intensity)
        T = KDTree( [(self.samples[i][0],self.samples[i][1]) for i in range(len(self.samples))] ) #build tree
        depth = [self.samples[i][2] for i in range(len(self.samples))]
        
        if branches=='all':
            print "Going through all branches."
            branches = self.gr.grid.keys()
        elif type(branches) != type([]):
            print "Please provide names of branches in a list.\n"
            return

        for branch in branches:
            print "Branch:", branch
            #Grid Polygons
            x = self.gr.grid[branch]['x']
            y = self.gr.grid[branch]['y']
            #dimensions
            gridheight = x.shape[0]
            gridlength = x.shape[1]
            
            hf = np.ones((gridheight, gridlength)) * (-999.0)
            
#TODO:      #AVERAGING: It can be performed as averaging the values within a circular/polygon/other neighbourhood
            #           Another way is to find the closest point in the neighbourhood (or the X closest points within a threshold)
                        #HERE: Closest one within cell circular neighbourhood is chosen, since data are usually much denser than the grid cells.
            print "Aggregating values..."
            for gh in range(gridheight):
                for gl in range(gridlength):
                    max_radius = []
                    try:
                        max_radius.append(self.__distance( [x[gh][gl],y[gh][gl]] , [x[gh-1][gl-1],y[gh-1][gl-1]] ))
                    except:
                        pass
                    try:
                        max_radius.append(self.__distance( [x[gh][gl],y[gh][gl]] , [x[gh+1][gl-1],y[gh+1][gl-1]] ))
                    except:
                        pass
                    try:
                        max_radius.append(self.__distance( [x[gh][gl],y[gh][gl]] , [x[gh+1][gl+1],y[gh+1][gl+1]] ))
                    except:
                        pass
                    try:
                        max_radius.append(self.__distance( [x[gh][gl],y[gh][gl]] , [x[gh-1][gl+1],y[gh-1][gl+1]] ))
                    except:
                        pass
                    #find max radius
                    if max_radius:
                        max_radius = max(max_radius)                    
                    else:
                        continue
                    #take half of the radius for circular neighbourhood (intesify by optional value of threshold)
                    r = 0.5 * max_radius * intensity
                    close_pnts = T.query_ball_point( (x[gh][gl],y[gh][gl]), r ) #find closest neighbouring samples
                    tempD = []
                    for p in close_pnts:
                        d = self.__distance( T.data[p], [x[gh][gl],y[gh][gl]] )
                        tempD.append([d,float(depth[p])])                        
                    #if at least 1 sample is found:
                    if tempD:
                        tempD.sort() #find closest one
                        if np.isnan(tempD[0][1]):
                            print "\nWARNING: Unexpected behaviour\n"
                        else:
                            hf[gh,gl] = tempD[0][1]
                print "-> longitudinal #{0} processed".format(gh)
            print "Finished aggregating."
            self.griddeddata[branch] = hf #store results
            
            
    def create_test(self, branch = '1', mode = 'cross', rowSkip = 'auto', columnSkip = 'auto', bankskip = [0,0]):
        '''
        Run after avgdep.
        '''
        
        print "\nCreating testing dataset."
        allowedmodes = ['cross','tracks']
        if mode not in allowedmodes:
            print "Invalid mode: Please specify either \'cross\' for cross-sections, or \'tracks\' for ship tracklines."
        
        try:
            x = self.gr.grid[branch]['x']
            y = self.gr.grid[branch]['y']
            z = copy.deepcopy(self.griddeddata[branch])
            self.test[branch] = {'x': x, 'y': y, 'depth': np.zeros(x.shape)-999.0}
            nocon = x.shape[0]
        except KeyError:
            print "Invalid branch name."
        
        #rows
        if rowSkip == 'auto':
            l = x.shape[0]
            rowSkip = int(np.floor(l/10.0))
            if rowSkip > l-2:
                rowSkip = l-2
        if rowSkip < 1:
            rowSkip = 0
        #columns
        if columnSkip == 'auto':
            l = x.shape[1]
            columnSkip = int(np.floor(l/20.0))
            if columnSkip > l-2:
                columnSkip = l-2
        if columnSkip < 1:
            columnSkip = 1  #at least 1 column distance spacing for testing dataset
        print "Row spacing: {0}\nColumn spacing: {1}\nBanklines skip: start({2}),end({3})".format(rowSkip,columnSkip,bankskip[0],bankskip[1])
    
        #cross-sections
        if mode == 'cross':
            print "->Mode: Cross-sections"
            print "->Row spacing denotes how many grid edge points are skipped vertically."
            print "->Column spacing denotes how many grid columns are skipped horizontally."
            #shifts to set the chosen points in the middle of the area
            Rshift = int( np.ceil((x.shape[0]-1)%(rowSkip+1) / 2.0) )     #row shift
            Cshift = int( np.ceil((x.shape[1]-1)%(columnSkip+1) / 2.0) )  #column shift  
            if rowSkip:
                temp = copy.copy(z[Rshift:z.shape[0]-Rshift+1:rowSkip+1,Cshift:z.shape[1]-Cshift+1:columnSkip+1])
                self.test[branch]['depth'][Rshift:z.shape[0]-Rshift+1:rowSkip+1,Cshift:z.shape[1]-Cshift+1:columnSkip+1] = temp
            else:
                temp = copy.copy(z[:,Cshift:z.shape[1]-Cshift+1:columnSkip+1])
                self.test[branch]['depth'][:,Cshift:z.shape[1]-Cshift+1:columnSkip+1] = temp
            
            if sum(bankskip)!=0:
                self.test[branch]['depth'][:bankskip[0],:] = -999.0
                self.test[branch]['depth'][nocon-bankskip[1]:,:] = -999.0
        #tracklines
        else:
            print "->Mode: Tracklines"
            print "->Row spacing denotes how many grid lines are skipped vertically."
            print "->Column spacing denotes how many grid columns are skipped horizontally."
            
            allowedrange = range(bankskip[0],nocon-bankskip[1])
            allowedrange = allowedrange + allowedrange[::-1]
            iterator = itertools.cycle(allowedrange)
            j = 0
            while j < x.shape[1]:
                i = iterator.next()
                self.test[branch]['depth'][i,j] = z[i,j]
                for k in range(rowSkip):
                    i = iterator.next()
                j += columnSkip  #(at least adds 1, because of ln.135-136)
                    
        print "*Testing dataset created successfully.*"
    
    
    def gridify_samples(self,branch):
        """
        Use to set the sample data on the grid, if it is known that the samples
        respond to the specified grid.
        """
        if self.samples:
            try:
                data = zip(*self.samples)[2] #retrive only z values
                b = str(branch)
                s = self.gr.grid[b]['x'].shape
                self.griddeddata[b] = common.to_grid(data, s[0], s[1])
                print "Success: Sample data set to grid notation (matrix).\n"
            except:
                print "ERROR: Could not transform samples onto gridded data."
        else:
            print "No data samples loaded."
               

    def plot(self, fignum = 'grid_dep', branch = '1', clim = 'auto', which = 'grid', classes = 100, system = 'xy', contour=True):
    
            
        if which == 'grid':
            hf = self.griddeddata[branch]
        elif which == 'test':
            hf =  self.test[branch]['depth']
        else:
            print "Please specify plot ('grid' or 'test')."
        
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
                    temp -= self.__distance([s1,n1],[s2,n2])
                    Gy[j,i] = temp
                temp = 0.0
                for j in range( (ly-1)/2 + 1, ly ):
                    s1, n1 = self.gr.grid[branch]['x'][j-1,i], self.gr.grid[branch]['y'][j-1,i]
                    s2, n2 = self.gr.grid[branch]['x'][j,i], self.gr.grid[branch]['y'][j,i]
                    temp += self.__distance([s1,n1],[s2,n2])
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
        
        if contour: #in case the data are almost fully covering the grid
            cs = plt.contourf(Gx, Gy, hf, levels, cmap=plt.cm.jet, extend="both")
        else: #in case the data are much scarcer than the grid cells
            cs = plt.scatter(Gx, Gy, s=5, c=hf, cmap=plt.cm.jet, vmin=clim[0], vmax=clim[1], edgecolors='none')
        cs.cmap.set_under('k')
        cs.set_clim(np.floor(levels[0]), np.ceil(levels[-1]))
        cbar = plt.colorbar(cs)
        cbar.set_label('Bed Level (m)'.format(), labelpad=15, weight='bold', size='14')
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
            print "Load successful.\n"
            
        elif extension == '.xyz':
            warnings.warn("Error: Cannot load sample file. Must be a depth file ('.dep')")
        else:
            warnings.warn("Error: Input file is not of an implemented extension!")
        
        self.griddeddata[branch] = A
        
        
        
    #Normal structure of 12 columns:  #?? very inefficient (O(n^3)), but correct
    def save_as_dep(self, outputfile, branch = '1', which = 'grid'):
        
        if which == 'grid':
            hf = self.griddeddata[branch]
        elif which == 'test':
            hf = self.test[branch]['depth']
        else:
            print "Wrong specified output."
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
        
    
    
    def save_as_samp(self, outputfile, branch = '1', which = 'grid'):        
    
        if which == 'grid':
            hf = self.griddeddata[branch]
        elif which == 'test':
            hf = self.test[branch]['depth']
        else:
            print "Wrong specified output."
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
        
        
        
    def save_as_shape(self, outputfile, branch = '1', which = 'grid'):
        if which == 'grid':
            hf = self.griddeddata[branch].flatten()
        elif which == 'test':
            hf = self.test[branch]['depth'].flatten()
        else:
            print "Wrong specified output."
        keep = (hf!=-999.0)
        z = hf[keep]
        x = self.gr.grid['1']['x'].flatten()[keep]
        y = self.gr.grid['1']['y'].flatten()[keep]
        
        w = shapefile.Writer(shapefile.POINT)
        
        sx = len(str(int(max(abs(x)))))+12
        sy = len(str(int(max(abs(y)))))+12
        sz = len(str(int(max(abs(z)))))+4
        
        w.field("X","N",sx,12)
        w.field("Y","N",sy,12)
        w.field("Z","N",sz,3)
        for i in range(len(z)):
            w.point(x[i],y[i])
            w.record(x[i],y[i],z[i])
        #small check
        outputfile = str(outputfile)
        if not(outputfile.endswith('.shp')):
            outputfile += '.shp'
        w.save(outputfile)

    
    
    # ----------------------------------------------------------------
    # Private Methods
    # ----------------------------------------------------------------

    def __read_file(self, inputfile, zfield):
        
        inputfile = inputfile.encode('string-escape')
        inputfile = inputfile.replace('\\','/')
        inputfile = inputfile.replace('//','/')
        print "Trying to load: "+str(inputfile)
        
        # Extract file extension
        extension = os.path.splitext(inputfile)[-1]
        allowed_extensions = ['.shp','.xyz','.dep']
        
        # Test if extension is allowed
        if extension not in allowed_extensions:
            raise KeyError('Extension not recognized. Invalid input file.')
            return
        
        #Load data:
        if extension == '.shp':
            #Read point data shapefile: The points must have attributes as X,Y,Z
            sf = shapefile.Reader(inputfile)
            if (sf.shapeType != 1):
                raise Exception('Input is not a POINT file!')
            data = sf.records()
            pnts = sf.shapes()
            #Close file:
            sf.shp.close()
            sf.dbf.close()
            sf.shx.close()
            self.samples = [ (pnts[i].points[0][0],pnts[i].points[0][1],data[i][zfield]) for i in range(len(pnts))  ]
        elif extension == '.xyz':
            self.samples = []
            with open(inputfile, 'r') as f:
                #START: Parsing
                for line in f:
                    words = line.split() #(x,y,z) space-delimited format
                    if (not words) or (words[0][0].isdigit() == False):
                        continue #skip empty lines or comment lines
                    self.samples.append( (float(words[0]),float(words[1]),float(words[zfield])) )
            f.close()
        elif extension == '.dep':
            print "Warning: Assuming match between GRID {0} and SAMPLES.".format(self.gr.grid.keys()[0])
            shape = self.gr.grid[self.gr.grid.keys()[0]]['x'].shape
            #input must be of same extents with the grid
            cols = 12 #standard Delft3D .dep file format of 12 columns
            lines_per_chunk = int((shape[1]+1)/cols)
            extra = (shape[1]+1)%cols
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
                h = np.array(data)
                #for each other chunk:
                for i in range(1,shape[0]): #excludes last chunk of -999.0's values
                    data = []
                    for l in range(lines_per_chunk):
                        line = f.readline()
                        temp = line.split()
                        data.extend([float(d) for d in temp])
                    h = np.vstack((h,data))
                h = h[:,:-1] #exclude last row and last column of -999.0's
            f.close()
            self.griddeddata[self.gr.grid.keys()[0]] = h #assume data on grid (if not, run avgdep())
            x = [i for i in self.gr.grid[self.gr.grid.keys()[0]]['x'].flatten()]
            y = [i for i in self.gr.grid[self.gr.grid.keys()[0]]['y'].flatten()]
            z = [i for i in h.flatten()]
            self.samples = zip(x,y,z)            
            
        print "Load successful.\n"

        
    def __distance(self, p1, p2):
        return np.sqrt( (p2[0]-p1[0])**2 + (p2[1]-p1[1])**2 )
