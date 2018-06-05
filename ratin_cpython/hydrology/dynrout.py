# -*- coding: utf-8 -*-
"""
.. module:: dynrout
    :platform: Windows
    :synopsis: Retrieve global discharge data


    
"""

import matplotlib as mpl
import matplotlib.dates as mdates
from datetime import datetime
modules = {}
try:
    from pydap.client import open_url
except ImportError:
    modules['pydap'] = True
    print 'pydap not installed!'
    

try:
    import numpy as np
except ImportError:
    print 'Numpy not installed!!'
  

class dynrout(object):
    '''
    This class provides an interface to the PCR-GLOBWB model results from the Deltares 
    OpenDAP server. Results can be extracted for a certain point defined by a
    latitude and longitude. 
    '''
    def __init__(self):    
        self.__set_default_parameters()
        #  GLOFRIS/PCR-GLOBWB dynrout results use the WGS84 (EPSG:4326) coordinate system
        self._CoordinateSystem = 'epsg:4326'
     
    
    def read_from_opendap(self,years):
        '''
        Read from the `OpenDAP <http://opendap.deltares.nl/thredds/dodsC/opendap/deltares/global_floods/EU-WATCH/dynrout/>`_ 
        server. 
        
        Args:
            years (list): list with year(s) for which to extract global discharge data
            
        Example:
        
        >>> dynrout.read_from_opendap([1999,2000])
        
        
        '''
       
        self.datalink = list()
        # Get data from OpenDAP server
        self.__get_url(years)
        for filename in self.file_url:
            print 'resolving '+filename
            data = open_url(filename)
            self.datalink.append(data)
        
            
    def extract_discharge_at_latlon(self,target_lat,target_lon):
        '''
        Extract data for a certain lat,long point. The method finds a point nearby
        from the model result grid. 'read_from_opendap' should be used before
        calling this method

        Args:
            target_lat (float): latitude
            target_lon (float): longitude
            
        Example:
        
        >>> # Extract discharge for the River Rhine 
        >>> dynrout.extract_discharge_at_latlon(51.74,6)
        
        '''
        if not(self.datalink):
            'First read data using read_from_opendap'
        else:
            print 'extracting at '+str(target_lat)+', '+str(target_lon)
            self.data = dict()
            self.data['time'] = []
            self.data['q']  = []
            for source in self.datalink:
                # Read latitude and longitude data
                lon = source['lon'][:]
                lat = source['lat'][:]
                
                

                # First item where lon > targetlon and lat is chosen as the closest to the
                # user defined approximate lat,lon coordinates
                approximate_lon = np.argmax(lon>target_lon)
                approximate_lat = np.argmax(lat>target_lat)
                
                # Extract data at the lat/lon combination
                self.data['time'].extend(source['time'][:])
                
                # Discharge comes out as nested lists and must be concatenated twice first
                q = source['qc']['qc'][:,approximate_lat,approximate_lon]            
                q = np.concatenate(q)
                q = np.concatenate(q)
                self.data['q'].extend(q)

    def plot(self):
        '''
        Simple plot.
        
        Be sure to call this method *after* calling extract_discharge_at_latlon,
        otherwise there is no data to plot. 
        
        This results in such a figure:
        
        .. figure::  _static/dynrout_rhine.png
           :width: 900px
           :align: center
           :height: 400px
        
        '''
 
        t = self.data['time']
        q = self.data['q']
        mpl.rc('font',**{'family':'serif','sans-serif':['Palatino']})
        time = []
        # Convert to datetime objects
        for i in t:
            time.append(datetime.fromordinal(int(i)+datetime.toordinal(datetime(1900,1,1))))
            
        fig,ax = mpl.pyplot.subplots(1)
        ax.plot(time,q,color='#551A8B',label='Discharge')
        
        ax.grid(True)
        mpl.pyplot.ylabel('Discharge [$m^3s^{-1}$]')
        mpl.pyplot.title('Rhine River at the border of The Netherlands and Germany')
        fig.autofmt_xdate()
        ax.fmt_xdata = mdates.DateFormatter('%y-%m-%d')
        mpl.pyplot.ion()
        mpl.pyplot.show()
        
        
    def __get_url(self,years):
        self.file_url = list()
        if self.parameters['default_fileformat'] == True:
            for year in years:
                self.file_url.append(self.opendap_url+str(year)+'1231_dynrout.nc')
        else:
            self.file_url.append(self.opendap_url+str(year))
                   
    def __set_default_parameters(self):
        self.opendap_url = r'http://opendap.deltares.nl/thredds/dodsC/opendap/deltares/global_floods/EU-WATCH/dynrout/'
        self.parameters = {}
        self.parameters['default_fileformat'] = True
        
    
        

        
        