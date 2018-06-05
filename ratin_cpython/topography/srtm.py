# -*- coding: utf-8 -*-
"""
.. module:: srtm
    :platform: Windows
    :synopsis: 
    
.. moduleauthor:: Koen Berends <koen.berends@deltares.nl>

    



"""
# Standard library imports
import os
import urllib
import requests
import zipfile
import subprocess as sbp
import glob
from matplotlib import cm
import matplotlib.pyplot as plt

# External dependencies
import numpy as np
from osgeo import gdal
import gdal_merge

class srtm(object):
    '''
    Load data from SRTM (Shuttle Radar Topography Mission)()
    
    Dependencies:
        * numpy (1.7+)
        * GDAL 
        
        `GDAL <www.gdal.org>`_ is an Open-Source library for geospatial raster 
        format translation. It comes bundled with the OGR library for geospatial 
        vector data. Look for installation instruction at `pypi <https://pypi.python.org/pypi/GDAL/>`_. 
        
        Tested with version 19 (MSVC2010)
    
    The code in this class is partly based on the OSM2Hydro source.      
    '''
    def __init__(self):
        self._set_defaults()
        
        # SRTM uses the WGS84 (EPSG:4326) coordinate system
        self._CoordinateSystem = 'epsg:4326'
        
    # ----------------------------------------------------------------
    # Public Methods
    # ----------------------------------------------------------------
    # These methods should be relatively easy accessable from outside    
        
    def set_srtm_url(self,mirror):
        '''
        Choose a srtm 4.1 mirror out of available mirrors:
        
        Args:
            Mirror (int) : 
                * #0 - droppr.org (DEFAULT)
                * #1 - srtm.csi.cgiar.org
        
        '''
        self.available_mirrors = {0: 'http://droppr.org/srtm/v4.1/6_5x5_TIFs/', \
                                  1:'ftp://srtm.csi.cgiar.org/SRTM_v41/SRTM_Data_GeoTIFF/'}                                  
        try:
            self.srtm_url = self.available_mirrors[mirror]
        except KeyError:
            print 'Mirror number is not availabe. Available mirrors:'
            for mirrornumber in self.available_mirrors:
                print str(mirrornumber)+': '+self.available_mirrors[mirrornumber]
            
    def build_srtm(self):
        '''
        Call this method to start building a SRTM tiff for the selected
        coordinates. Set coordinates:
        
        >>> # Set coordinates for Liege, Belgium
        >>> srtm.coords['longitude'] = [5.25,5.75]
        >>> srtm.coords['latitude'] = [50.53,50.75]
        >>> # Call build
        >>> srtm.build()
        
        Output is being written to the output folder. Set outputfolder using:
        
        >>> srtm.outputdir = r'C:\tmp'
        
        By default, the output directory is located in the workdirectory/rat_output
        '''
        self._retrieve_srtm()
        self._merge_srtm()  
    
    def transform_coordinates(self,new_coordinatesystem,inputfile = 'none',outputfile='none'):
        '''
        Transform coordinates of geotif data using GDAL. SRTM data is by default
        encoded in the WGS84 coordinate system (EPSG:4326). However, for most
        appliances a cartesian (WGS84 is geodetic) is needed. 
        
        In the following example a geotif is build (an active internet 
        connection is generally required) and converted to a local coordinate
        system. In this case, the SRTM covers the city of Arnhem in the
        Netherlands. The local coordinate system is 'Rijksdriehoek' (EPSG:28992)
        
        >>> from rattool.topography import srtm
        >>> dem = srtm.srtm()
        >>> # Set bounds of model area
        >>> dem.coords['latitude'] = [51.965,52.03]
        >>> dem.coords['longitude'] = [5.9,6.1]
        >>> # Build SRTM DEM 
        >>> dem.build_srtm()
        >>> dem.contour(fignum = 1)
        >>> dem.transform_coordinates('epsg:28992')
        >>> dem.contour(fignum=2)
        
        This will produce two figures. In the output folder, a projected geotif
        is written. By default, the projected (transformed) tif has the 
        suffix '_proj'. The .tif file is accompanied by a .prj file containing
        information about the coordinate system of the tif. 
        
        The tif file can be loaded by most GIS application, e.g. ARC-GIS. The 
        following figure shows the projected map in GIS, with supporting
        additional information about the river geometry. 
        
        .. figure::  _static/arnhem_srtm.png
           :width: 500px
           :align: center
           :height: 400px
           
           `SRTM for the Rhine at Arnhem, The Netherlands, in Rijksdriehoek coordinate system`
        
        '''
        if inputfile == 'none':
            inputfile = self.merged_srtm_file
        if outputfile == 'none':
            inputfile_part,ext = os.path.splitext(inputfile)
            outputfile = os.path.join(inputfile_part+'_proj'+ext)
        
        # Transform .tif file        
        requestsObject = self._getCSInfo(self._CoordinateSystem,'proj4')
        oldstring = requestsObject.text
        requestsObject = self._getCSInfo(new_coordinatesystem,'proj4')
        newstring = requestsObject.text
        command= 'gdalwarp -s_srs "%s" -t_srs "%s" %s %s' % (oldstring, newstring, inputfile, outputfile)
        print command
        sbp.call(command)
        
        # Write accompanying .prj file containing the coordinate system string
        outputfilepart,ext = os.path.splitext(outputfile)
        output_prj = outputfilepart+'.prj'
        requestsObject = self._getCSInfo(new_coordinatesystem,'proj4')
        with open(output_prj,'w') as f:
            f.write(requestsObject.text)
        
        self.merged_srtm_file = outputfile
        return
        
    @staticmethod
    def _getCSInfo(cs,info):
        csp = cs.split(':')
        string= r'http://spatialreference.org/ref/'+csp[0]+\
                                          r'/'+csp[1]+r'/'+info+'/'
        print string                               
        requestsObject = requests.get(string) 
        return requestsObject      

           
    def convert_tif(self,inputfile,outputfile,outputfiletype):    
        
        filetype_options = {'xyz':
        'gdal_translate -of XYZ '+inputfile+' '+ outputfile,\
        'ascii':\
        'gdal_translate -of AAIGrid '+inputfile+' '+outputfile}
        
        try:
            print 'Converting to '+outputfiletype+'...'
            sbp.call(filetype_options[outputfiletype])
            print 'File saved to '+ os.getcwd()
        except KeyError:
            'Illegal filetype. Available filetypes:'
            for key in filetype_options:
                print key  
                
    def getValueAtCoords(self,xy,inputfile='none'):
        '''
        Extract values from GeoTiff at given x,y coordinates.
        Args:
            xy (list): [[x],[y]] coordinates (or lon/lat depending on coordinate system)
        Kwargs:
            inputfile (str): GeoTiff file. 
        '''
        if inputfile == 'none':
            inputfile = self.merged_srtm_file
        if type(xy[0]) == type([]):
            if not(len(xy[0]) == len(xy[1])):
                print 'Error: x and y not of equal length'
                return
        
        ds = gdal.Open(inputfile)
        dem = ds.ReadAsArray()
        gt = ds.GetGeoTransform()

        xres = gt[1]
        yres = gt[5]
        
        X = np.linspace(gt[0], gt[0] + dem.shape[1]*xres, dem.shape[1])
        Y = np.linspace(gt[3], gt[3] + dem.shape[0]*yres, dem.shape[0])
        
        output_z = []
        if type(xy[0]) == type([]):
            for i in range(len(xy[0])):
                nearest_x = np.argmin(np.abs(X-xy[0][i]))
                nearest_y = np.argmin(np.abs(Y-xy[1][i]))
                output_z.append(dem[nearest_y,nearest_x])
        else:
            nearest_x = np.argmin(np.abs(X-xy[0]))
            nearest_y = np.argmin(np.abs(Y-xy[1]))
            output_z.append(dem[nearest_y,nearest_x])
            
        return output_z
            
        
          
    def plot3d(self,inputfile = 'none'):
        '''
        Plot a geospatial tiff. Supply file name. If  build_srtm has been used, 
        and no kwarg is supplied, plots the output from build_srtm. 
        
        Note that this method is potentially very slow. An example output
        result is shown below
        
        
        .. figure::  _static/srtm_plot3d.png
           :width: 500px
           :align: center
           :height: 400px
           
           `SRTM for Meuse at Liege, Belgium (longitude 5.25-5.7, latitude 50.53-50.75)`

        
        Kwargs:
            FileIn (str): filename (include absolute path)
            
        Dependencies:
            * matplotlib
            
        '''
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        if inputfile == 'none':
            inputfile = self.merged_srtm_file
                
        ds = gdal.Open(inputfile)
        dem = ds.ReadAsArray()
        gt = ds.GetGeoTransform()
        ds = None
        fig,ax = plt.subplots(figsize=(16,8), subplot_kw={'projection': '3d'})
        xres = gt[1]
        yres = gt[5]
        
        X = np.linspace(gt[0], gt[0] + dem.shape[1]*xres, dem.shape[1])
        Y = np.linspace(gt[3], gt[3] + dem.shape[0]*yres, dem.shape[0])
        
        X, Y = np.meshgrid(X, Y)
        
        surf = ax.plot_surface(X,Y,dem, rstride=1, cstride=1, cmap=plt.cm.RdYlBu_r, vmin=0, vmax=250, linewidth=0, antialiased=True)
        ax.set_zlim(0, 600) # to make it stand out less
        ax.view_init(60,-105)
        
        fig.colorbar(surf, shrink=0.4, aspect=20)
        
        plt.show() 
    
    def contour(self,inputfile = 'none',fignum=1):
        '''
        Plots a countour plot of geotiff data
        
        Kwargs:
            FileIn (str): filename (include absolute path)
            
        Dependencies:
            * matplotlib
        '''
        import matplotlib.pyplot as plt
        from matplotlib import cm
        if inputfile == 'none':
            inputfile = self.merged_srtm_file
                
        ds = gdal.Open(inputfile)
        dem = ds.ReadAsArray()
        gt = ds.GetGeoTransform()
        ds = None
        fig = plt.figure(fignum)
        ax =  fig.add_subplot(111)
        xres = gt[1]
        yres = gt[5]
        
        X = np.linspace(gt[0], gt[0] + dem.shape[1]*xres, dem.shape[1])
        Y = np.linspace(gt[3], gt[3] + dem.shape[0]*yres, dem.shape[0])
        
        X, Y = np.meshgrid(X, Y)
        
        cs = ax.contourf(X,Y,dem,30, antialiased=True,cmap = cm.gist_earth)
        ax.contour(X,Y,dem,10, antialiased=True,colors='k')
        ax.clabel(cs,inline=1,fontsize=10)

        plt.ion()
        plt.show() 
        
    # ----------------------------------------------------------------
    # Private Methods
    # ----------------------------------------------------------------
    # These methods are only for internal use

    def _retrieve_srtm(self):
        '''
        This method download and unzips SRTM tiles        
        '''
        
        # Get area for which to extract SRTM
        xmin = self.coords['longitude'][0]
        xmax = self.coords['longitude'][1]
        ymin = self.coords['latitude'][0]
        ymax = self.coords['latitude'][1]
            
        # Get SRTM tiles from latitude/longitude 
        tileMinX=np.int((np.round(xmin) + 180 ) / 5 + 1)
        tileMaxX=np.int((np.round(xmax) + 180 ) / 5 + 1)
        tileMinY=np.int((60 - np.round(ymax) ) / 5 + 1)
        tileMaxY=np.int((60 - np.round(ymin) ) / 5 + 1)        


        tileLat=tileMinY-1
        tileLon=tileMinX-1
        
        # Loop through the required SRTM tiles, download and unzip
        for tileLon in range(tileMinX, tileMaxX+1):
            for tileLat in range(tileMinY, tileMaxY+1):
                
                fileName = str(self.file_prefix + '%02.f_%02.f' + self.url_suffix) % (tileLon, tileLat)
                url = self.srtm_url + fileName
                fileTarget = os.path.join(self.outputdir, fileName)

 
                # If the file is not locally already available, download
                if not(os.path.exists(fileTarget)):
                    print 'Retrieving from ' + url
                    urllib.urlretrieve(url,fileTarget,reporthook = self._dlprogress)
                else:
                    print('Target file locally found. Not downloading from external source')
                    
                print('Unzipping %s' %(fileTarget))
                zf = zipfile.ZipFile(fileTarget , 'r') 
                nameList = zf.namelist()
                for n in nameList:
                    outFile = open(os.path.join(self.outputdir, n), 'wb')
                    outFile.write(zf.read(n))
                    outFile.close()
                zf.close()


 
    def _set_defaults(self):
        self.set_srtm_url(0)
        self.coords = {}
        self.coords['longitude'] = [5.25,5.7]
        self.coords['latitude'] = [50.533,50.75]
        self.file_prefix='srtm_'
        self.url_suffix='.zip'
        self.__oldprcnt=-1 # for progress bar
        # Output directory
        self.outputdir = os.path.join(os.getcwd(),'rat_output')
        try:
            os.mkdir(self.outputdir)
        except WindowsError:
            print 'outputfolder already exists'
            
        print 'Printing output to ' + self.outputdir


    def _cutMap(self,InMap,OutMap):
        """
        asdf
        """
        command= 'gdal_translate -ot Float32 -projwin %f %f %f %f %s %s' \
             % (self.coords['longitude'][0], self.coords['latitude'][1],\
                self.coords['longitude'][1], self.coords['latitude'][0], InMap, OutMap)
        print command  
        try:
            sbp.call(command,shell=True)
        except OSError as msg:
            print 'RATTOOLS | Execution failed: ' + msg
            print 'Is GDAL installed?'

    def _merge_srtm(self):
        temporary_dem = os.path.join(self.outputdir, 'temp.tif')
        cut_dem       = os.path.join(self.outputdir, '_latlon.tif')
        source_dems   = os.path.join(self.outputdir, 'srtm*.tif')
        
        if not os.path.isdir(os.path.join(self.outputdir)):
            raise Exception("Output folder does not exist! Tip: run retrieve_strm()")
        
        gdal_merge.main(argv=['dummy','-o', temporary_dem, source_dems])
        
        # this is the final lat lon map          
        self._cutMap(temporary_dem,cut_dem) 

        # remove all redundant files
        self._removeFiles(os.path.join(self.outputdir, 'srtm*.tif'))
        self._removeFiles(os.path.join(self.outputdir, 'srtm*.hdr'))
        self._removeFiles(os.path.join(self.outputdir, 'srtm*.tfw'))
        os.unlink(os.path.join(self.outputdir, 'readme.txt'))
        os.unlink(os.path.join(self.outputdir, 'temp.tif'))
        
        self.merged_srtm_file = cut_dem
        
    
    def _removeFiles(self,wildCard):
        filelist = glob.glob(wildCard)
        for filename in filelist:
            os.unlink(filename)
    


            
    def _dlprogress(self,count, blockSize, totalSize):
        percent = int(count*blockSize*100/totalSize)
         
        if np.mod(percent,10)==0 and not(self.__oldprcnt == percent):
            print "..."+str(percent)+"% \r",
        self.__oldprcnt = percent   
