#!/usr/bin/env python
# coding: utf-8
"""
.. module:: virtualearth
    :platform: Windows
    :synopsis: Provides the medial_axis_transform_class
    
.. moduleauthor:: Koen Berends <koen.berends@deltares.nl>

This module provides the VirtualEarth class that allows the user to plot 
images from Microsoft's Bing Map servers using the Bing Tile System. This system
is excellenty explained `here <http://msdn.microsoft.com/en-us/library/bb259689.aspx>`_. 

Dependencies
    * numpy (1.7+)
    * matplotlib    
"""
import os
import matplotlib.image as mpimg
import numpy as np
import random
import urllib

class VirtualEarth():
    '''
    Load images from virtual earth. Save to disk and plot if desired. Images
    are only downloaded if there is no local copy available. 
    
    
    '''
    def __init__(self):

        # Options        
        self.layers = 'hybrid'
        self.outputdir = './rat_output'
        self.urltemplate = 'http://ecn.t{4}.tiles.virtualearth.net/tiles/{3}{5}?g=0'     
        self.quality = 'medium'
        
        # no options
        self._get_quality()
        self.layerdict = {'satellite': 'a', 'hybrid': 'h', 'roads': 'r'}
        self.imdict = {}
        self.imglatlon = {}
        self.surfdict = {}
        
    
    def plot(self,lat,lon):
        '''
        Plot a selected location
        
        Args:
            lat (list) : latitude [min,max]
            lon (list) : longitude [min,max]
            
        Downloaded temporary files are stored in VirtualEarth.outputdir
        '''
        import matplotlib.pyplot as plt
        coords = {'lat':lat,'lon':lon}
        self._get_zoom_level(coords)
        self._get_tile_xy(coords)
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for img in self.imdict:
            ax.imshow(self.imdict[img],extent=self.imglatlon[img]['lon']+self.imglatlon[img]['lat'],origin = 'upper',aspect =1.7)
            
        ax.set_xlim(coords['lon'])
        ax.set_ylim(coords['lat'])
        
        ax.grid(color = 'w')
        ax.patch.set_alpha(0)
        self.fig = fig
        self.ax = ax
        plt.ion()
        plt.show()

    def setquality(self,quality):
        '''
        Set quality of images. The quality determines the zoom level. 
        
        Args:
            quality (str): quality of the image. ['low','medium','high']
            
        
        >>> from rattool.topography import virtualearth
        >>> ts = virtualearth.VirtualEarth()
        >>> ts.outputdir = './rat_output/'
        >>> ts.setquality('high')
        >>> ts.plot([51.90,51.98],[5.90,6])
        
        This produces following image:
        
        .. figure::  _static/ve_arnhem.png
           :width: 500px
           :align: center
           :height: 400px
           
        `Hybrid Bing Maps result for the bifurcation of the Pannerden Canal into
        the River Nederrijn and River IJssel at Arnhem`
           
        '''
        self.quality = quality
        self._get_quality()
        

    # ----------------------------------------------------------------
    # Public Methods
    # ----------------------------------------------------------------
    # These methods should be relatively easy accessable from outside    
    def _get_quality(self):
        qualitydict = {'low': 256, 'medium':512, 'high':1024}
        self.target_resolution = qualitydict[self.quality.lower()]
    # ----------------------------------------------------------------
    # Private Methods
    # ----------------------------------------------------------------
    # These methods are for internal use
       
    def _get_zoom_level(self,coords):
        dlat = coords['lat'][1]-coords['lat'][0]
        dlat_virtualearth = 85.05112878*2 
        map_width = np.ceil(dlat_virtualearth/dlat)*self.target_resolution
        level_of_detail = np.floor(np.log(map_width/256)/np.log(2))
        self.level_of_detail = level_of_detail

    def _get_tile_xy(self,coords):
        x_tileminmax = []
        y_tileminmax = []
        for Lon in coords['lon']:
            for Lat in coords['lat']: 
                pixelX,pixelY = self._latlon_to_xy(Lat,Lon,self.level_of_detail)

                tileX = np.floor(pixelX/256)
                tileY = np.floor(pixelY/256)
                
                x_tileminmax.append(int(tileX))             
                y_tileminmax.append(int(tileY))

        x_tilerange = range(np.min(x_tileminmax),np.max(x_tileminmax)+1)
        y_tilerange = range(np.min(y_tileminmax),np.max(y_tileminmax)+1)


        for tileX in x_tilerange:
            for tileY in y_tilerange:             
                tilekey = (int(tileX), int(tileY), int(self.level_of_detail))
                self.imdict[tilekey] = self._tile_as_image(tilekey)
                self.imglatlon[tilekey] = {}
                lonmin,latmax = self._xy_to_latlon(tileX*256,tileY*256,self.level_of_detail)
                lonmax,latmin = self._xy_to_latlon((tileX+1)*256,(tileY+1)*256,self.level_of_detail)
                self.imglatlon[tilekey]['lat'] = [latmin,latmax]
                self.imglatlon[tilekey]['lon'] = [lonmin,lonmax]
                self.imglatlon[tilekey]['x'] = [tileX*256,(tileX+1)*256]
                self.imglatlon[tilekey]['y'] = [(tileY+1)*256,(tileY)*256]
     
    def _latlon_to_xy(self,Lat,Lon,n):
        # see http://msdn.microsoft.com/en-us/library/bb259689.aspx for source if these formulas
        sinLatitude = np.sin(float(Lat)*np.pi/180)
        pixelX = ((float(Lon)+180)/360)*256*2**n
        pixelY = (0.5-np.log((1+sinLatitude)/(1-sinLatitude))/(4*np.pi))*256*2**n
        
        return pixelX,pixelY
        
    def _xy_to_latlon(self,x,y,n):
        
        lon = x*360/(256*2**n)-180
        alpha = np.exp(4*np.pi*(y*256**-1*2**-n-0.5))
        slat = (1-alpha)/(alpha+1)
        lat = np.arcsin(slat)*180/np.pi
        
        return lon,lat
        
        
    def _tiletoquadkey(self, xi, yi, z):
        quadKey = ''
        for i in range(z, 0, -1):
            digit = 0
            mask = 1 << (i - 1)
            if(xi & mask) != 0:
                digit += 1
            if(yi & mask) != 0:
                digit += 2
            quadKey += str(digit)
            #print quadKey
        return quadKey

    def _loadimage(self, fullname, tilekey):
        im = mpimg.imread(fullname)
        self.imdict[tilekey] = im
        return self.imdict[tilekey]

    def _tile_as_image(self,tilekey):
        xi,yi,zoom = tilekey
        result = None
        try:
            result = self.imdict[tilekey]
        except:
            filename = '{}_{}_{}_{}.jpg'.format(zoom, xi, yi, self.layerdict[self.layers])
            fullname = os.path.join(self.outputdir,filename)
            try:
                result = self._loadimage(fullname, tilekey)
            except:
                server = random.choice(range(1,4))
                quadkey = self._tiletoquadkey(*tilekey)
                #print quadkey
                url = self.urltemplate.format(xi, yi, zoom, self.layerdict[self.layers], server, quadkey)
                print "Downloading tile %s to local cache." % filename
                urllib.urlretrieve(url, fullname)
                result = self._loadimage(fullname, tilekey)
        return result

