# -*- coding: utf-8 -*-
"""
.. module:: mat
    :platform: Windows
    :synopsis: Provides the medial_axis_transform_class
    
.. moduleauthor:: Koen Berends <koen.berends@deltares.nl>

Algorithm background
-------------------------------
The medial axis of an object is the theoretical center-line. It
is also known as the *morphological skeleton* or *topological skeleton*. 

There are various ways of estimating the medial axis. An often-used method
is achieved by 'morphological shrinking' of the bitmap representation of the
object. Another method links the midpoints of lines resulting from delauny
triangulation.

Both methods are problematic in some way when applied to rivers. Morphological
shrinking requires very high resolution bitmap representation. Delauny
Triangulation needs sophisticated masking.

The aim is to present a robust, easy-to-use method to determine not only
the medial axis (center-line), but the river width as well. Therefore the
implementation makes use of the concept algorithm proposed by Vilaplana [1]_.

The algorithm makes use of 'testing discs' to find a way through the river system.

.. [1] Josep Vilaplana (2013), Computing the Medial Axis Transform of Polygonal Objects by Testing Discs 


"""

import numpy as np
import os
import matplotlib as mpl
from ratin.common import common


class MedialAxisTransformation:
    ''' 
    This class supplies the following attributes for configuration:
    
    * search_tolerance: Controls the width of the ring in meters
    * search_set_size: Controls the active domain in meters
    * step_size: Relative size of the radius of the active disc. Used when computing the next step
    * outlier_tolerance: ???
    '''

    def __init__(self):

        # These values are set by default
        self.search_tolerance   = 50.
        self.search_set_size    = 1500. #half-edge of square search area #? why not make this parametric to the active disc size?
        self.step_size          = 0.3   #how much the bubble "shifts" in the chosen direction
        self.current_direction  = 160.
        self._proposed_direction = 90.
        self.outlier_tolerance  = 110.
        self.minimal_width      = 200.  #below this, the bubble is too small and the algorithm ends
        self.maximum_transient_collection_size = 500
        self.maximal_width      = 2000. #over this, the bubble is too large and the algorithm ends
        self.shrink_parameter = 5.0 #how much the bubble shrinks
        #self.flush_memory_limit = 1000  #?unused
        self.output_folder      = os.getcwd()
        self.usepolygon         = True

        # Initial values. Should not be tampered with
        self.jump_start       = False #?? if not manually tampered, it never changes
        self.firstdisc        = True
        self.breakflag        = False
        self.flag_outlier     = False #?? unused/unchanged
        self.flag_grow        = False
        self._figholdon       = False
        self.advice           = 'none'
        self.breakmessage     = 'no breakmessage'


        # Initiate the warm and cold libraries
        self.active_disc = {'x':[],'y':[],'r':[],'d':[],'edges':[]}
        
        self.proposed_disc = self.active_disc.copy()
        self.disc_library = {}
        self.disc_library['discs'] = []
        self.disc_library['medial_axis'] = {}
        self.disc_library['medial_axis']['x'] = []
        self.disc_library['medial_axis']['y'] = []
        self.disc_library['medial_axis']['width'] = []
        self.disc_library['medial_axis']['theta'] = []

    # ----------------------------------------------------------------
    # Public Methods
    # ----------------------------------------------------------------
    # These methods should be relatively easily accessible from outside
    def load(self, first, second=False):
        '''
        Load a pointcloud or polygon 
        :param first:
        :param second:
        Args:
            first: a list of x values or polygon
            second: a list of y values (if first contains x values)
           
        '''
        if first and second:
            raise Exception('x,y input not supported anymore!')            
            self.x_collection = np.array(first)
            self.y_collection = np.array(second)
        elif first:
            self.usepolygon = True
            self.polygon = first

            self.x_collection = list()
            self.y_collection = list()
            self.p_collection = list()
            self.pi_collection= list()

            #if self.usepolygon: #? why is this 'if' here? (5 lines above it is set to True)
            # No support (yet) for more polygons!
            if len(self.polygon.data['data']) > 1:
                print 'More than 1 polygon Detected. Functionality not implemented!'+\
                      '\n Using first polygon in dict'
            poly = self.polygon.data['data'].keys()[0] #first polygon's name

            partnumber = 0
            for polypart in self.polygon.data['data'][poly]['parts']:
                self.x_collection.extend(polypart['x'])
                self.y_collection.extend(polypart['y'])
                self.p_collection.extend(np.array(polypart['x']).astype(int)**0*partnumber)
                self.pi_collection.extend(range(0,len(polypart['x'])))
                polypart['pi'] = range(0,len(polypart['x']))
                partnumber += 1

            self.x_collection = np.array(self.x_collection)
            self.y_collection = np.array(self.y_collection)
            self.p_collection = np.array(self.p_collection)
            self.pi_collection = np.array(self.pi_collection)

#            #?An other way (to avoid the 'for' loop and reassigning to arrays) [lines 125-137]:
#            self.x_collection = np.array(self.polygon.data['data'][poly]['x'])[~(np.isnan(self.polygon.data['data'][poly]['x']))]
#            self.y_collection = np.array(self.polygon.data['data'][poly]['y'])[~(np.isnan(self.polygon.data['data'][poly]['y']))]
#            temp = [self.polygon.data['data'][poly]['parts'][i]['pi'] for i in range(len(self.polygon.data['data'][poly]['parts']))]
#            from itertools import chain #place on top
#            self.p_collection = np.array(list(chain(*[len(temp[i])*[i] for i in range(len(temp))])))
#            self.pi_collection = np.array(list(chain(*temp)))

            #?? What is transient collection?
            self.transient_collection = {}
            self.transient_collection['x'] = np.array([])
            self.transient_collection['y'] = np.array([])
            self.transient_collection['p'] = np.array([])
            self.transient_collection['pi'] = np.array([])
            
            #?unused ??
#            self.original_x_collection = np.array(self.x_collection)
#            self.original_y_collection = np.array(self.y_collection)
#            self.original_p_collection = np.array(self.p_collection)
#            self.original_pi_collection = np.array(self.pi_collection)
            
    def reset(self):
        '''
        Resets the object to a pristine state. To be used after 'build'
        '''
        self.active_disc = {'x':[],'y':[],'r':[],'d':[]}
#        self.x_collection = self.original_x_collection
#        self.y_collection = self.original_y_collection
#        self.p_collection = self.original_p_collection
#        self.pi_collection = self.original_pi_collection
        self.disc_library = {}
        self.disc_library['discs'] = []
        self.disc_library['medial_axis'] = {}
        self.disc_library['medial_axis']['x'] = []
        self.disc_library['medial_axis']['y'] = []
        self.disc_library['medial_axis']['width'] = []
        self.disc_library['medial_axis']['theta'] = []
        #self.current_direction = [] #??why assign here? at ln ?? is another
        self._proposed_disc = self.active_disc.copy()
        self._proposed_direction = [0]
        self._blobtangent = [0]
        self.advice         =  'none'
        self.firstdisc        = True
        self.breakflag        = False
        self._breakmessage    = 'no breakmessage'
        self.flag_outlier     = False
        self.flag_grow        = False
        self.current_direction = 160.
        self.polylines = []  #??this is the first instance, does not exist in init
        #??Why reasigning so many parameters? Why not to simply recall init?

    def initialise(self, x, y, r, d):
        '''
        Start from a specific disc. If the attribute 'jump_start' is True, 
        it will immediatly propose a new disc. This is necessary if the 
        initial disc (specified by x,y,r and d) is also a valid disc. 
        
        Args:
            x (float): x-coordinate of the centre of the disc 
            y (float): y-coordinate of the centre of the disc 
            r (float): radius of the disc 
            d (float): path direction of the disc [in degrees!] 
            
        '''
        x = np.float(x)
        y = np.float(y)
        r = np.float(r)
        d = np.float(d)

        self.active_disc['x'] = x
        self.active_disc['y'] = y
        self.active_disc['r'] = [r,r+self.search_tolerance,self.search_set_size] #? I think the third value can be safely omitted (the global search_set_size can be used always)
        self.active_disc['d'] = d
        self._proposed_direction = d
        self.current_direction = d
        self._update_disc()
        #? jump_start - unless specified by the user - is never changing to True!
        if self.jump_start:
            # commit and directly propose new bubble
            self._commit_to_library()
            self.active_disc['x'] += np.cos(np.radians(d))*self.active_disc['r'][0]*self.step_size
            self.active_disc['y'] += np.sin(np.radians(d))*self.active_disc['r'][0]*self.step_size
            self._update_disc()
            self._test_disc()

    def build(self, maximum_iterations, edge_determined = False):
        '''
        Start building the medial axis. 

        '''
        if not(edge_determined):
            print '\n-----------(((MAT ALGORITHM)))-----------'
            for iteration in range(maximum_iterations):
                #self.outeriteration = iteration #? commented below, unused??
                # Print every 50 iteration
                if np.mod(iteration,50)==0:
                    print str(iteration)+'/'+str(maximum_iterations)
                #If everything is ok
                if not(self.breakflag):
                    if self.advice != 'none':
                        self._propose_disc()
                        self._check_progress()
                    self._edge_detection()
                    self._test_disc()
    
                elif self.breakflag:
                    print ('Algorithm ended with message: '+self.breakmessage)
                    print '--------------(((END MAT)))--------------\n'
                    break
    
    
            if iteration == maximum_iterations-1:
                print('Algorithm ended with message: maximum iterations achieved')
                print '--------------(((END MAT)))--------------\n'
                self._breakcode = 4
                self.breakdirections = [self._proposed_direction]
                self.breakwidth = 2*self.active_disc['r'][0]
                
        else:
            raise Exception('x,y input not supported anymore!')     


    def plot(self, fignum=1, block=False):
        '''
        Plots the medial axis. Call after build. Dependent on matplotlib
        
        Kwargs:
            fignum (int): figure number (default = 1)
            block (boolean): blocks the figure after plotting. Use this if the plot closes immediately after plotting.
        '''
        #region Create Figure
        fig = mpl.pyplot.figure(fignum)
        ax = fig.add_subplot(111)
        #endregion

        #region Plot medial axis seeker circle
        xl = mpl.pyplot.xlim()
        yl = mpl.pyplot.ylim()

        x_inner = self.active_disc['inner']['x']
        y_inner = self.active_disc['inner']['y']
        x_ring = self.active_disc['ring']['x']
        y_ring = self.active_disc['ring']['y']
        x = self.active_disc['x']
        y = self.active_disc['y']
        d = self.active_disc['d']
        r = self.active_disc['r']

        plot_circle1 =mpl.patches.Circle((self.active_disc['x'] ,self.active_disc['y']),self.active_disc['r'][1],facecolor='none',edgecolor='k' )
        ax.add_patch(plot_circle1)
        #endregion

        #region Plot vertices of the polygon
        ax.plot(self.x_collection, self.y_collection, '.', color=[0.6, 0.6, 0.6])
        ax.plot(x_inner,y_inner,'.',color=[0.0,0.2,0.2])
        ax.plot(x_ring,y_ring,'.',color=[1.0,0,0])
        ax.plot(self.disc_library['medial_axis']['x'],self.disc_library['medial_axis']['y'],'.-',color='r')
        #endregion\\

        #region Plot additional circle parameters
        minmax = [self._proposed_direction+self.outlier_tolerance, self._proposed_direction+2*self.outlier_tolerance]
        minmax[minmax > 360] -= 360
        try:
            ax.plot([x, x+np.cos(np.radians(d))*r[0]], [y, y+np.sin(np.radians(d))*r[0]], 'k')
        except TypeError:
            pass
        except:
            for i in d:
                ax.plot([x,x+np.cos(np.radians(i))*r[0]],[y,y+np.sin(np.radians(i))*r[0]],'k')
        for i in range(0,360,5):
            ax.plot([x+np.cos(np.radians(i))*r[0]*0.9,x+np.cos(np.radians(i))*r[0]],
                     [y+np.sin(np.radians(i))*r[0]*0.9,y+np.sin(np.radians(i))*r[0]],'k')
        for i in range(0,360,45):
            ax.plot([x+np.cos(np.radians(i))*r[0]*0.9,x+np.cos(np.radians(i))*r[0]],
                     [y+np.sin(np.radians(i))*r[0]*0.9,y+np.sin(np.radians(i))*r[0]],'r')
        for i in [minmax[0],minmax[1]]:
            try:
                ax.plot([x+np.cos(np.radians(i))*r[0]*0.9,x+np.cos(np.radians(i))*r[0]],
                         [y+np.sin(np.radians(i))*r[0]*0.9,y+np.sin(np.radians(i))*r[0]],'m',linewidth =4)
            except:
                pass
        try:
            for i in self._blobtangent:
                ax.plot([x+np.cos(np.radians(i))*r[0]*0.9,x+np.cos(np.radians(i))*r[0]],
                         [y+np.sin(np.radians(i))*r[0]*0.9,y+np.sin(np.radians(i))*r[0]],'g',linewidth = 4)
        except:
            pass

        if self._figholdon:
            mpl.pyplot.xlim(xl)
            mpl.pyplot.ylim(yl)

        plot_circle2 =mpl.patches.Circle((self.active_disc['x'] ,self.active_disc['y']),self.active_disc['r'][0],facecolor='none',edgecolor='b',linestyle='dashed' )
        ax.add_patch(plot_circle2)
        #endregion

        #region Plot width
        direction = common.get_direction(self.disc_library['medial_axis']['x'],
                                         self.disc_library['medial_axis']['y'],
                                         smoothdegree=5.0)   ##?? why 5 here???

        leftline = common.get_parallel_line(self.disc_library['medial_axis']['x'],
                                            self.disc_library['medial_axis']['y'],
                                            direction,
                                            np.array(self.disc_library['medial_axis']['width'])/2.)

        mpl.pyplot.plot(leftline[0], leftline[1],'-g')

        rightline = common.get_parallel_line(self.disc_library['medial_axis']['x'],
                                            self.disc_library['medial_axis']['y'],
                                            direction,
                                            -1*np.array(self.disc_library['medial_axis']['width'])/2.)

        mpl.pyplot.plot(rightline[0], rightline[1],'-g')
        #endregion
        #region Show figure
        mpl.pyplot.ion()
        mpl.pyplot.show(block=block)
        #endregion


    # ----------------------------------------------------------------
    # Private Methods
    # ----------------------------------------------------------------
    # Just keeping up appearances

    def _propose_disc(self):
        '''
        This method takes the 'advice' string from self._test_disc and proposes
        a new disc. 
        '''
        self.previous_disc = self.active_disc.copy() #?WHY? "previous_disc" is used nowhere. Also, whatever list it has copied will change when active_disc changes

        # Possible actions
        proposition = {
                        'grow'       : self._grow_bubble,
                        'move_grow'  : self._move_and_grow_bubble,
                        'commit+new' : self._commit_and_propose,
                        'move_shrink': self._move_and_shrink_bubble, #? not used
                        'move'       : self._move_bubble,
                        'shrink'     : self._shrink_bubble,
                        'none'       : 1+1,
                        }

        # Perform proposed action
        #print self.advice
        proposition[self.advice]()

    def _check_progress(self):
        '''
        Perform additional mid-loop checks on progress. Breaks if the disc
        becomes greater or smaller than its limits
        '''
        if self.active_disc['r'][0]*2 <= self.minimal_width:
            self.breakflag = True
            self._breakcode = 5
            self.breakdirections = [None, None]
            self.breakwidth = [None, None]
            self.breakmessage = 'Disc smaller than minimal width ('+str(self.active_disc['r'][0]*2)+'<'+str(self.minimal_width)
        if self.active_disc['r'][0]*2 > self.maximal_width:
            self.breakflag = True
            self._breakcode = 3
            self.breakmessage = 'Disc greater than maximal width ('+str(self.active_disc['r'][0]*2)+'>'+str(self.maximal_width)+')'

            
    def _edge_detection(self):
        '''
        This function does multiple things:
        
        1. 
        Determine the intersection points of the active disc
        with its edges. The edges are defined by the parameters (a,b,xmin,xmax)
        where y = ax+b for xmin<x<xmax. 
        
        2. 
        Build polylines by finding out which edges are connected. The polylines
        serve as the 'tangent lines' of the disc. 
        
        '''
#        self._update_disc()
    
        # Declarations
        intersection_points = dict()
        part_numbers = []

        # Get the polygon parts near the disc
        for edge in self.active_disc['edges']:
            part_numbers.append(edge['part'])           

        part_numbers = np.unique(part_numbers)

        # Every part gets its own dictionary entry
        for part in part_numbers:    
            intersection_points[str(part)] = [[],[],[],[]] #x,y,p,pi
        
        # Outer Ring parameters
        xc = self.active_disc['x'];
        yc = self.active_disc['y'];
        r = self.active_disc['r'][1];

        #==============================================================================
        # Loop through edges and detect intersection points with the disc
        #==============================================================================
        counter = 0
        for edge in self.active_disc['edges']:
            counter+=1
            
            # Edge parameters (y = ax+b for xmin<x<xmax)
            a=edge['a']
            b=edge['b'] 
            xmin=edge['xmin']
            xmax=edge['xmax']
            if abs(a) > 5000: #? Shouldn't it be abs(a) > 5000 ? changed it
                # Edge cannot be used to determine intersection points, assume 
                # perfectly vertical line
                  
                x = xmin
                p = edge['part']
                pi = edge['partindex']
                y = self.active_disc['ring']['y'][(self.active_disc['ring']['p']==p) & (self.active_disc['ring']['pi']==pi)]    
                
                #??? What are these continuous if's?
                if len(y) == 0:
                    pi = edge['partindex2']
                    y = self.active_disc['ring']['y'][(self.active_disc['ring']['p']==p) & (self.active_disc['ring']['pi']==pi)]    
                if len(y) == 0:
                    pi = edge['partindex']
                    y = self.active_disc['inner']['y'][(self.active_disc['inner']['p']==p) & (self.active_disc['inner']['pi']==pi)] 
                if len(y) == 0:
                    pi = edge['partindex2']
                    y = self.active_disc['inner']['y'][(self.active_disc['inner']['p']==p) & (self.active_disc['inner']['pi']==pi)]

                    
                print '----------------------------'
                print a
                print 'from '+str(edge['partindex'])+' to '+str(edge['partindex2'])
                print y

                theta = np.angle((x-xc)+(y-yc)*1j)

                print theta                
                delta_y = np.sin(theta)*self.active_disc['r'][1]*2
                print delta_y
                print '----------------------------'
                
                if np.sqrt((x-xc)**2+(y+delta_y)**2) <= self.active_disc['r'][1]:
                    intersection_points[str(edge['part'])][0].append(x)
                    intersection_points[str(edge['part'])][1].append(y+delta_y)
                    intersection_points[str(edge['part'])][2].append(p)
                    intersection_points[str(edge['part'])][3].append(pi)
                else:
                    intersection_points[str(edge['part'])][0].append(x)
                    intersection_points[str(edge['part'])][1].append(y-delta_y)
                    intersection_points[str(edge['part'])][2].append(p)
                    intersection_points[str(edge['part'])][3].append(pi)
            
            else:
                # Discriminant to find intersections between line of line segment
                # and outer ring
                A = a**2+1
                B = 2*(a*b-a*yc-xc)
                C = xc**2+yc**2-r**2-2*b*yc+b**2
                discriminant = B**2-4*A*C
           
                if discriminant < 0: #no intersection
                    pass 
                elif discriminant == 0: #1 intersection == tangent
                    pass
                elif discriminant > 0: #2 intersections
                    # Points of intersections. 
                    x1 = (-B+np.sqrt(discriminant))/(2*A)
                    y1 = a*x1+b
                    x2 = (-B-np.sqrt(discriminant))/(2*A)
                    y2 = a*x2+b
                    
                    #TODO: ? review again - resolve cases?                    
                    if (x1>xmin)&(x1<xmax): #if first point is within the line segment (within the circle as well)
                        intersection_points[str(edge['part'])][0].append(x1)
                        intersection_points[str(edge['part'])][1].append(y1)
                        intersection_points[str(edge['part'])][2].append(edge['part'])
                        intersection_points[str(edge['part'])][3].append(edge['partindex'])
                    if (x2>xmin)&(x2<xmax): #if second point is within the line segment (within the circle as well)
                        intersection_points[str(edge['part'])][0].append(x2)
                        intersection_points[str(edge['part'])][1].append(y2)
                        intersection_points[str(edge['part'])][2].append(edge['part'])
                        intersection_points[str(edge['part'])][3].append(edge['partindex'])
        
                    #? what if the linestring is fully within the ring?
        self.intersection_points = intersection_points
        
        # The intersection points are stored to the object's polygon. This is 
        # necessary to navigate sparsely pointed polygons. In this way the disc
        # 'touches the walls' whilst trying to find a way through the river
        # corridor. The stepsize should be smaller than the discs radius to avoid
        # losing touch.
        if intersection_points:
            # Add intersection points to ring and collection (touching the wall)
            for part in intersection_points:
                for i in reversed(range(len(intersection_points[part][0]))):
                    x = intersection_points[part][0][i]
                    y = intersection_points[part][1][i]
                    pi = intersection_points[part][3][i]
                    d = np.sqrt((x-self.active_disc['x'])**2+(y-self.active_disc['y'])**2)
                    theta = np.angle((x-self.active_disc['x'])+(y-self.active_disc['y'])*1j)
                    
                    self.active_disc['ring']['x']=np.append(self.active_disc['ring']['x'],x)
                    self.active_disc['ring']['y']=np.append(self.active_disc['ring']['y'],y)
                    self.active_disc['ring']['d']=np.append(self.active_disc['ring']['d'],d)
                    self.active_disc['ring']['theta']=np.append(self.active_disc['ring']['theta'],theta)
                    self.active_disc['ring']['p']=np.append(self.active_disc['ring']['p'],int(part))
                    self.active_disc['ring']['pi']=np.append(self.active_disc['ring']['pi'],pi)

                    self._append_to_collection({'x':x,'y':y,'pi':pi,'p':int(part)})
                    
        self._update_disc()
        #==============================================================================
        #  Determine which edges are connected and build polylines
        #==============================================================================        
        polylines = []
        current_polyline = -1
        # Loop through intersection points
        for part in part_numbers:
            current_polyline +=1
            polylines.append({'pi':[],'x':[],'y':[]})
            
            sortindex   = np.argsort(self.active_disc['ring']['pi'][self.active_disc['ring']['p']==part])
            partx       = self.active_disc['ring']['x'][self.active_disc['ring']['p']==part][sortindex]
            party       = self.active_disc['ring']['y'][self.active_disc['ring']['p']==part][sortindex]
            partindexes = self.active_disc['ring']['pi'][self.active_disc['ring']['p']==part][sortindex]

            maxpi = np.max(self.pi_collection[self.p_collection==part])

            if len(partindexes)==1:
                # Now there is only 1 point (tangent point, no tangent line)
                polylines[current_polyline]['y']=party
                polylines[current_polyline]['x']=partx
                polylines[current_polyline]['pi']=partindexes
            
            else:
                for j in range(1,len(partindexes)):
                    polylines[current_polyline]['x'].append(partx[j-1])
                    polylines[current_polyline]['y'].append(party[j-1])
                    polylines[current_polyline]['pi'].append(partindexes[j-1])
                    if np.abs(partindexes[j]-partindexes[j-1]) <= 1:
                        polylines[current_polyline]['x'].append(partx[j])
                        polylines[current_polyline]['y'].append(party[j])
                        polylines[current_polyline]['pi'].append(partindexes[j])
                    elif (partindexes[j] == maxpi) & (np.min(polylines[current_polyline]) == 0):
                        polylines[current_polyline]['x'].append(partx[j])
                        polylines[current_polyline]['y'].append(party[j])
                        polylines[current_polyline]['pi'].append(partindexes[j])
                    else:
                        # Check if the gap between the polylines is big enough
                        # to think that they should be two seperate lines. 
                        distance = np.sqrt((partx[j]-partx[j-1])**2 + (party[j]-party[j-1])**2)
                        if distance > self.minimal_width:                       
                            current_polyline +=1
                            polylines.append({'pi':[partindexes[j]],'x':[partx[j]],'y':[party[j]]})

              
        self.polylines = polylines   

    def _append_to_collection(self,appendDict):
        '''
        Append to the collection (the 'polygon point cloud')
        '''
        
        # Where in the list to insert the new point   
#        insert_at_index = np.where((self.pi_collection == appendDict['pi']) & (self.p_collection==appendDict['p']))
#        insert_at_index = insert_at_index[0][0]
        # Give other indexes a higher index number
        pi_plus_one = (self.pi_collection >= appendDict['pi']) & (self.p_collection==appendDict['p'])
        self.pi_collection[pi_plus_one] = self.pi_collection[pi_plus_one] + 1
        
        # Insert new point in collection
        self.x_collection = np.append(self.x_collection,appendDict['x'])
        self.y_collection = np.append(self.y_collection,appendDict['y'])
        self.p_collection = np.append(self.p_collection,appendDict['p'])
        self.pi_collection = np.append(self.pi_collection,appendDict['pi'])

        
        for key in self.transient_collection:
            length_test = len(self.transient_collection[key]) - self.maximum_transient_collection_size
            if length_test > 0:
                self.transient_collection[key] = self.transient_collection[key][length_test:-1]
        

        self.transient_collection['x'] = np.append(self.transient_collection['x'],appendDict['x'])
        self.transient_collection['y'] = np.append(self.transient_collection['y'],appendDict['y'])
        self.transient_collection['p'] = np.append(self.transient_collection['p'],appendDict['p'])
        self.transient_collection['pi'] = np.append(self.transient_collection['pi'],appendDict['pi']+0.001)

        # Insert new point in transient
        
#        self.x_collection = np.insert(self.x_collection,insert_at_index,appendDict['x'])
#        self.y_collection = np.insert(self.y_collection,insert_at_index,appendDict['y'])
#        self.p_collection = np.insert(self.p_collection,insert_at_index,appendDict['p'])
#        self.pi_collection = np.insert(self.pi_collection,insert_at_index,appendDict['pi'])

    def _test_disc(self):
        '''
        This method tests the proposed disc (self.active_disc) and returns an 
        advice (grow, move and grow, shrink, etc)
        '''
        # Aliases to be used in function
        x_ring      = self.active_disc['ring']['x']
        x_inner     = self.active_disc['inner']['x']

        # Are there no vertices in the inner disc?
        if not x_inner.size:
            # Are there no vertices in the ring?
            if not x_ring.size:
                # no vertices in outer ring either
                # Are we dealing with an outlier?
                if self.flag_outlier: #? will never enter this if because "_check_for_outlier()" is no longer used
                    self.breakmessage = 'Outlier, discard partial ma'
                    self._breakcode = 2
                    self.breakflag = True
                else:
                    self.advice = 'grow'
            else:
                # Vertices in outer ring, not in inner ring (yes!)
                # How many theta blobs? #?what are these???
                if len(self.active_disc['ring']['theta']) > len(self.active_disc['ring']['x']):
                    raise Exception('Active Disc corrupt (more theta than x values)')

                number_of_blobs = self._count_tangents('ring')
                if number_of_blobs == 1:
                    # alpha of the tangent
                    self._blobtangent= self._blob_tangents()
                    self._proposed_direction,forget,forget,forget= self._determine_direction(self._blobtangent)
                    self.advice = 'move_grow'

                elif number_of_blobs == 2:
                    self._blobtangent = self._blob_tangents()
                    # Check if blobs are behind current direction (circle of 
                    # +/- 110 degrees in either direction)
                    self._check_for_outlier() #? useless - "pass"-only function
                    # New direction
                    new_direction,direction_change,new_width,forget = self._determine_direction(self._blobtangent)

                    # Direction change
                    if  direction_change > 50: #? why is this the threshold and not a smaller one?
                        # Sudden course direction, grow first. try to find bifurcation or smth
                        # Advice to not follow this direction and grow instead
                        self._proposed_direction = np.mean(self.active_disc['ring']['theta'])
                        self.advice = 'move_grow'
                    else:
                        self.advice = 'commit+new'
                        self._proposed_direction = new_direction

                elif number_of_blobs > 2:
                    # Bifurcation candidate
                    self._blobtangent = self._blob_tangents()
                    new_direction,direction_change,break_width,break_coordinates = self._determine_direction(self._blobtangent)

                    self._proposed_direction = new_direction
                    self._blobtangent = self._blob_tangents()
                    self._commit_to_library()#x,y,p,pi

                    self.breakdirections = new_direction
                    self.breakwidth = break_width
                    self.breakx = break_coordinates[0]+self.active_disc['x']
                    self.breaky = break_coordinates[1]+self.active_disc['y']
                    self.breakflag = True
                    self._breakcode = 1
                    self.breakmessage = 'Bifurcation candidate detected'


        else:
            # There are vertices in the inner disc
            
            number_of_blobs = self._count_tangents('inner')
            self._blobtangent = self._blob_tangents(area = 'inner')
            number_of_blobs = len(self._blobtangent)
            
            if number_of_blobs == 1:
                self.advice = 'move'
            elif number_of_blobs > 1:
                self.advice = 'shrink'

#==============================================================================
#     Actions
#==============================================================================

    def _find_width(self):
        #region definitions
        midx = self.active_disc['x']
        midy = self.active_disc['y']
        x = list()
        y = list()
        #endregion

        for line in self.polylines:
            x.append(line['x'][len(line['x'])/2])
            y.append(line['y'][len(line['x'])/2])

            dis2cent = [np.sqrt((line['x'][i]-midx)**2 + (line['y'][i]-midy)**2) for i in range(len(line['x']))]

        return min(dis2cent)*2

    def _commit_and_propose(self):
        '''
        Commit the active disc to the disc library, propose new disc
        '''

        # Before committing, check real with
        width = self._find_width()

        # Commit to medial axis dict
        self.disc_library['medial_axis']['x'].append(self.active_disc['x'])
        self.disc_library['medial_axis']['y'].append(self.active_disc['y'])
        #self.disc_library['medial_axis']['width'].append(self.active_disc['r'][0]*2.)
        self.disc_library['medial_axis']['width'].append(width)
        self.disc_library['medial_axis']['theta'].append(self.current_direction)

        # The newly proposed disc receives the new direction        
        self.active_disc['d'] = self._proposed_direction

        # Commit the disc to library          
        self.disc_library['discs'].append(self.active_disc.copy())

        # First disc clearly does not apply anymore
        self.firstdisc = False

        # Propose new disc
        self.active_disc['x'] = self.active_disc['x']+np.cos(np.radians(self._proposed_direction))*self.active_disc['r'][0]*self.step_size
        self.active_disc['y'] = self.active_disc['y']+np.sin(np.radians(self._proposed_direction))*self.active_disc['r'][0]*self.step_size
        self.current_direction = self._proposed_direction
        self._update_disc()

    def _commit_to_library(self):
        '''
        Commit the active disc to the disc library, propose new disc
        '''
        self.disc_library['medial_axis']['x'].append(self.active_disc['x'])
        self.disc_library['medial_axis']['y'].append(self.active_disc['y'])
        self.disc_library['medial_axis']['width'].append(self.active_disc['r'][0]*2.)
        self.disc_library['medial_axis']['theta'].append(self.current_direction)

        # The newly proposed disc receives the new direction        
        self.active_disc['d'] = self._proposed_direction

        # Commit the disc to library          
        self.disc_library['discs'].append(self.active_disc.copy())

        # Firstdisc clearly does not apply anymore
        self.firstdisc = False

    def _update_disc(self, getset=True, xset = [], yset=[], pset = [], piset = []):
        '''
        This function 'updates the active disc' by looking what is in or near 
        the active disc. 
        
        1.
        This function selects vertices from the 'polygon point cloud' that are
        near the active disc and finds out if they are within the disc ('inner')
        or in the outer ring ('ring')
        
        2. 
        The function constructs edges based on these points. 
        
        '''
        #? how about using shapely? :)
        # Declarations
        self.active_disc['edges'] = [] #clear any existing stored edges
        
        # Aliases
        inner_radius = self.active_disc['r'][0]
        outer_radius = self.active_disc['r'][1] #inner_radius+tolerance
        x_center   = self.active_disc['x']
        y_center   = self.active_disc['y']

        # Perform rest of calculation on a subset
        # The subset is formed by the closest points to the disc center, in a
        # square area with half-edge equal to the search_set_size and centroid
        # the disc center
        x_set,y_set,p_set,pi_set = self._get_subset()

        # Distances of points in set to active disc center
        set_distances = np.sqrt((x_set-x_center)**2+(y_set-y_center)**2)

        #Out of the subset, choose only points that lie within the ring
        part_numbers_in_disc = p_set[(set_distances < outer_radius)]
        part_indexes_in_disc= pi_set[(set_distances < outer_radius)]
        part_numbers = np.unique(p_set[(set_distances < outer_radius)])
     
#==============================================================================
#         Add vertices just outside disc (for extended edge detection).  #? what is meant here?
#==============================================================================
        p_selection = np.array([])
        pi_selection = np.array([])
        for i in part_numbers:
            sorted_indexes = np.sort(part_indexes_in_disc[part_numbers_in_disc==i])

            pimax = np.max(self.pi_collection[self.p_collection==i])
            pimin = np.min(self.pi_collection[self.p_collection==i])

            #? Why this way?
            add_to_indexes = []
            for j in range(len(sorted_indexes)):
                if sorted_indexes[j] == pimax: #assumes closed polygon
                    add_to_indexes.append(0)
                else:
                    add_to_indexes.append(sorted_indexes[j]+1)
                if sorted_indexes[j] == pimin:
                    add_to_indexes.append(pimax)
                else:
                    add_to_indexes.append(sorted_indexes[j]-1)
            sorted_indexes = np.append(sorted_indexes,add_to_indexes)
            sorted_indexes = np.unique(sorted_indexes)
            pi_selection = np.append(pi_selection,sorted_indexes)
            p_selection = np.append(p_selection,i*np.ones(len(sorted_indexes)))
            
        # Find the x,y,p values for these points
        x_selection = np.array([])
        y_selection = np.array([])
        for i in range(len(p_selection)):
            x_selection = np.append(x_selection,x_set[(p_set==p_selection[i])&(pi_set==pi_selection[i])])
            y_selection = np.append(y_selection,y_set[(p_set==p_selection[i])&(pi_set==pi_selection[i])])

#==============================================================================
#         Make edges
#==============================================================================
        #? Why constantly duplicating the same arrays?
        for i in part_numbers:
            pi = pi_selection[p_selection==i]
            xset = x_selection[p_selection==i]
            yset = y_selection[p_selection==i]
            
            for j in range(1,len(pi)):
                if (np.abs(pi[j]-pi[j-1]) <= 1) or (pi[j]==pimax):
                    if pi[j] == 0:
                        x1 = xset[pi == pi[j]][0]
                        x2 = xset[pi == pi[-1]][0]
                        y1 = yset[pi == pi[j]][0]
                        y2 = yset[pi == pi[-1]][0]
                        second_partindex = pi[-1]

                    else:
                        x1 = xset[pi == pi[j-1]][0]
                        x2 = xset[pi == pi[j]][0]
                        y1 = yset[pi == pi[j-1]][0]
                        y2 = yset[pi == pi[j]][0]
                        second_partindex = pi[j-1]

                    if x1 == x2: 
                        # Steep climber. The formula y=ax+b does not apply for
                        # such cases. So we do a little cheat by moving one
                        # of the x coordinates a unit to the right (east)
                        #? you could still assign a=float("inf") etc, and it would still be 'caught' as special case in _edge_detection()
                        x2 = x1+1.0

                    a = (y1-y2)/(x1-x2)
                    b = y1-(a*x1)
                    xmin = np.min([x1, x2])
                    xmax = np.max([x1, x2])
                    #? maybe explain why this notation is good/needed?
                    self.active_disc['edges'].append({'a':a,'b':b,'xmin':xmin,'xmax':xmax,'part':int(i),'partindex':int(pi[j]),'partindex2':int(second_partindex)})
    
        # Vertices within inner circle
        p_inner = p_set[(set_distances < inner_radius)]
        pi_inner= pi_set[(set_distances < inner_radius)]
        x_inner = x_set[(set_distances < inner_radius)]
        y_inner = y_set[(set_distances < inner_radius)]
        inner_dist = set_distances[(set_distances < inner_radius)]
        inner_theta = np.angle((x_inner-x_center)+(y_inner-y_center)*1j)

        # Vertices in outer ring
        x_ring = x_set[(set_distances >= inner_radius) & (set_distances <= outer_radius) ]
        y_ring = y_set[(set_distances >= inner_radius) & (set_distances <= outer_radius) ]
        if self.usepolygon:
            p_ring = p_set[(set_distances >= inner_radius) & (set_distances <= outer_radius) ]
            pi_ring= pi_set[(set_distances >= inner_radius) & (set_distances <= outer_radius) ]
        ring_dist = set_distances[(set_distances >= inner_radius) & (set_distances <= outer_radius) ]
        ring_theta = np.angle((x_ring-x_center)+(y_ring-y_center)*1j)

#==============================================================================
#         Assign to active disc
#==============================================================================
        self.active_disc['set'] = {}
        self.active_disc['set']['x'] = x_selection
        self.active_disc['set']['y'] = y_selection
        self.active_disc['set']['p'] = p_selection.astype(int)
        self.active_disc['set']['pi'] = pi_selection.astype(int)
        self.active_disc['inner'] = {}
        self.active_disc['inner']['x'] = x_inner
        self.active_disc['inner']['y'] = y_inner
        self.active_disc['inner']['d'] = inner_dist
        self.active_disc['inner']['theta'] = inner_theta
        self.active_disc['ring'] = {}
        self.active_disc['ring']['x'] = x_ring
        self.active_disc['ring']['y'] = y_ring
        if self.usepolygon:
            self.active_disc['ring']['p'] = p_ring.astype(int)
            self.active_disc['ring']['pi'] = pi_ring.astype(int)
            self.active_disc['inner']['p'] = p_inner.astype(int)
            self.active_disc['inner']['pi'] = pi_inner.astype(int)
        self.active_disc['ring']['d'] = ring_dist
        self.active_disc['ring']['theta'] = ring_theta


    def _get_subset(self):
        '''
        Get a subset of the 'polygon point cloud' (*_collection), to speed
        up calculations
        '''
        # Initialize arrays
        x_set  = np.array([])
        y_set = np.array([])
        p_set = np.array([]) # partnumber
        pi_set = np.array([]) #part-indexnumber

        # Two Clouds: Collection (= original ldb)
        distances_x= np.abs(self.x_collection-self.active_disc['x'])
        distances_y= np.abs(self.y_collection-self.active_disc['y'])

        selection = (distances_x < self.search_set_size) & (distances_y < self.search_set_size)
        #? Unsure if it will have the same functionality... (choosing only points in circle)
#        x_set = self.x_collection[selection]
#        y_set = self.y_collection[selection]
#        p_set = self.p_collection[selection]
#        pi_set = self.pi_collection[selection]
        
        #???? HERE THE CHOICE IS ALL POINTS OF A PART THAT LIES WITHIN THE CIRCLE - ALL POINTS OF THAT PART (=POLYGON)!
        #???? Instead, I wrote the above which basically selects ONLY the points within the circle -
        #???? BUT: If the edges of the polygon are too long it might happen that it finds no points...
        #???? A more effiecient way would be to check edges - shapely can help with that

        parts_in_selection = np.unique(self.p_collection[selection])
        for partnumber in parts_in_selection:
            x_set = np.append(x_set,self.x_collection[self.p_collection==partnumber])
            y_set = np.append(y_set,self.y_collection[self.p_collection==partnumber])
            p_set = np.append(p_set,self.p_collection[self.p_collection==partnumber])
            pi_set = np.append(pi_set,self.pi_collection[self.p_collection==partnumber])
        
#        # Two Clouds: Transient (= collection of intersection points)
#        distances_x= np.abs(self.transient_collection['x']-self.active_disc['x'])
#        distances_y= np.abs(self.transient_collection['y']-self.active_disc['y'])
#
#        selection = (distances_x < self.search_set_size) & (distances_y < self.search_set_size)
#        parts_in_selection = np.unique(self.transient_collection['p'][selection])
#        for partnumber in parts_in_selection:
#            x_set = np.append(x_set,self.transient_collection['x'][self.transient_collection['p']==partnumber])
#            y_set = np.append(y_set,self.transient_collection['y'][self.transient_collection['p']==partnumber])
#            p_set = np.append(p_set,self.transient_collection['p'][self.transient_collection['p']==partnumber])
#            pi_set = np.append(pi_set,self.transient_collection['pi'][self.transient_collection['p']==partnumber])

        return np.array(x_set),np.array(y_set),np.array(p_set),np.array(pi_set)


    def _check_for_outlier(self):
        '''
        This method is the 'silent guardian' checking if the disc might 
        have slipped through the veil. 
        
        This function should not be called anymore if using edge detection. 
        '''
        pass
#        self._blobtangent = self._blob_tangents()
#        minmax = [self.current_direction-self.outlier_tolerance, self.current_direction+self.outlier_tolerance]
#
#        minmax[minmax<0] = minmax[minmax<0]+360
#        minmax[minmax>360] = minmax[minmax<0]-360
#        in_slice = []
#            in_slice.append(self._degree_in_slice(t,minmax[0],minmax[1]))
#
#        in_slice = np.array(in_slice)
#        if in_slice.all():
#            self.flag_outlier= True
#        else:
#            self.flag_outlier = False

    def _degree_in_slice(self,degree,slicemin,slicemax): #? not used
        if slicemin > slicemax:
            if degree > slicemin:
                verdict = True
            elif degree < slicemax:
                verdict = True
            else:
                verdict = False
        else:
            if slicemin < degree < slicemax:
                verdict = True
            else:
                verdict = False
        return verdict

    def _theta2alpha(self,theta):   #? not used
        alpha = theta/(2*np.pi)*360
        alpha[alpha<0] = alpha[alpha<0] + 360
        return alpha

    def _count_tangents(self,area = 'ring'):
        '''
        Find number of tangents. 
        
        note: this used to be a convulated function but is trimmed down when 
        switching to polylines. In future maybe re-add the 'area' functionality. 
        '''
        number_of_blobs = len(self.polylines)
        return number_of_blobs

    def _shrink_bubble(self):
        # Number of blobs in ring
        inner_dist = self.active_disc['inner']['d']
        new_radius = np.min(inner_dist)-self.shrink_parameter
        if new_radius <= 0.0:
            print "The bubble has collapsed on itself while shrinking"
            raise ValueError
        self.active_disc['r'] = [new_radius,new_radius+self.search_tolerance,self.search_set_size]
        
        # Update the active (= proposed) disc
        self._update_disc()

    def _grow_bubble(self):
        new_radius = 1.5*self.active_disc['r'][0]
        self.active_disc['r'] = [new_radius,new_radius+self.search_tolerance,self.search_set_size]
        
        # Update the active (= proposed) disc
        self._update_disc()

    def _move_bubble(self):
        '''
        Move disc in opposite direction of one blob in inner circle
        '''
        move_dist = np.max([self.active_disc['r'][0]- np.min(self.active_disc['inner']['d']),0.1])
        tangent = self._blobtangent[0]
        if tangent >= 0:
            reverse_alpha = tangent-180
        else:
            reverse_alpha = 180+tangent
            
        #print 'move distance: '+str(move_dist) + ' in ' +str(reverse_alpha)
        self.active_disc['x'] = self.active_disc['x'] + np.cos(np.radians(reverse_alpha))*move_dist
        self.active_disc['y'] = self.active_disc['y'] + np.sin(np.radians(reverse_alpha))*move_dist
        
        # Update the active (= proposed) disc
        self._update_disc()

    def _move_and_grow_bubble(self):
        tangent = self._blobtangent[0]
        old_radius = self.active_disc['r'][0]
        new_radius = old_radius+0.05*self.active_disc['r'][0]
        if tangent >= 180:
            reverse_alpha = tangent-180
        else:
            reverse_alpha = tangent+180

        self.active_disc['x'] = self.active_disc['x'] + np.cos(np.radians(reverse_alpha))*(new_radius-old_radius)
        self.active_disc['y'] = self.active_disc['y'] + np.sin(np.radians(reverse_alpha))*(new_radius-old_radius)
        self.active_disc['r'] = self.active_disc['r'] = [new_radius,new_radius+self.search_tolerance,self.search_set_size]
        
        # Update the active (= proposed) disc
        self._update_disc()


    def _move_and_shrink_bubble(self):  #? not used
        '''
        Note: this function cannot be called anymore since rev. 60. It is kept 
        until it's clear it is not needed anymore
        '''
        
        tangent = self._blobtangent[0]

        old_radius = self.active_disc['r'][0]
        inner_blob_distance = old_radius-np.min(np.array(self.active_disc['inner']['d']))

        move_dist = 0.5*inner_blob_distance+0.25*self.search_tolerance
        new_radius = old_radius-0.5*inner_blob_distance

        if tangent >= 180:
            reverse_alpha = tangent-180
        else:
            reverse_alpha = tangent+180

        self.active_disc['x'] = self.active_disc['x'] + np.cos(np.radians(reverse_alpha))*move_dist
        self.active_disc['y'] = self.active_disc['y'] + np.sin(np.radians(reverse_alpha))*move_dist
        self.active_disc['r'] = self.active_disc['r'] = [new_radius,new_radius+self.search_tolerance,self.search_set_size]
        
        # Update the active (= proposed) disc
        self._update_disc()

    def _base_length(self,beta):
        '''
        Calculate the base length of an isosceles triangle. 
        beta = the angle in degrees between the isosceles
        
        returns the length of the base of the triangle
        '''
        delta = 2*self.active_disc['r'][0]*np.sin(np.radians(beta))
        return delta
        
    def _blob_tangents(self,area='ring'):
        '''
        Returns the tangents of polylines in the disc ring        
        '''
        if area == 'ring':
            tangent = []
            for line in self.polylines:
                dx =  np.mean([np.min(line['x']),np.max(line['x'])]) - self.active_disc['x']
                dy =  np.mean([np.min(line['y']),np.max(line['y'])]) -self.active_disc['y']
                alpha = np.angle(dx+dy*1j, True)
                if alpha < 0:
                    alpha += 360
                tangent.append(alpha)
        
        elif area == 'inner':
            polylines = self.make_polylines('inner')
            tangent = []
            for line in polylines:
                if len(line['x']) > 1:
                    linex = np.mean([np.min(line['x']),np.max(line['x'])])
                    liney = np.mean([np.min(line['y']),np.max(line['y'])]) 
                else:
                    linex = self.active_disc['inner']['x'][0]
                    liney = self.active_disc['inner']['y'][0]

                dx =  linex - self.active_disc['x']
                dy =  liney - self.active_disc['y']
                alpha = np.angle(dx+dy*1j, True)
                if alpha < 0:
                    alpha += 360
                tangent.append(alpha)
        
        return tangent
        
    def make_polylines(self,area):
        polylines = {}               
        partnumbers = np.unique(self.active_disc[area]['p'])
        polylines = []
        current_polyline = -1
        # Loop through intersection points
        for part in partnumbers:
            current_polyline +=1
            polylines.append({'pi':[],'x':[],'y':[]})
            # The connection number groups intersection points to the same or different connected polyline
            
            sortindex = np.argsort(self.active_disc[area]['pi'][self.active_disc[area]['p']==part])
            partx  = self.active_disc[area]['x'][self.active_disc[area]['p']==part][sortindex]
            partindexes =  self.active_disc[area]['pi'][self.active_disc[area]['p']==part][sortindex]
            party =  self.active_disc[area]['y'][self.active_disc[area]['p']==part][sortindex]
            
    
            for j in range(1,len(partindexes)):
                if np.abs(partindexes[j]-partindexes[j-1]) <= 1:
                    polylines[current_polyline]['x'].append(partx[j])
                    polylines[current_polyline]['y'].append(party[j])
                    polylines[current_polyline]['pi'].append(partindexes[j])
                else:
                    current_polyline +=1
                    polylines.append({'pi':[],'x':[],'y':[]})
                        
        return polylines

    def _determine_direction(self,alpha):
        '''
        Returns the direction and width between two blobs. Two directions are calculated,
        the one closest to the original direction is picked. 
        '''
        alpha = np.array(alpha)

        if len(alpha)>2:
            # More than 3 alphas, so bifurcation candidates
            candidates = []
            alpha = np.sort(alpha)
            alpha = np.append(alpha,alpha[0]+360)

            # Calculate width between blobs
            AlphaBetweenBlobs = np.diff(alpha)
            Radius = self.active_disc['r'][0]#np.mean(self.active_disc['r'])
            WidthBetweenBlobs = 2*Radius*np.sin(np.radians(0.5*AlphaBetweenBlobs))

            # Calculate coordinates between blobs
            dXBetweenBlobs = Radius*np.cos(np.radians(AlphaBetweenBlobs))
            dYBetweenBlobs = Radius*np.sin(np.radians(AlphaBetweenBlobs))

            for i in range(0,len(alpha)-1):
                candidates.append((alpha[i]+alpha[i+1])/2)

            # Now delete the one which is most likely the direction we're coming from (so opposite from current direction)
            CurrentDirection = np.array(self.current_direction)
            OppositeDirection = np.copy(CurrentDirection)

            if OppositeDirection < 180:
                OppositeDirection = OppositeDirection +180
            else:
                OppositeDirection = OppositeDirection - 180

            DirectionCandidates = np.abs(np.array(OppositeDirection)-np.array(candidates))


            # 3 and 360 degrees are very near. The maximal distance of one candidate is always 180 degrees
            DirectionCandidates[DirectionCandidates > 180] = 360 - DirectionCandidates[DirectionCandidates > 180]
            DirectionComeFrom = np.argmin(DirectionCandidates)
            OutputDirection = np.delete(candidates,DirectionComeFrom)
            OutputWidth = np.delete(WidthBetweenBlobs,DirectionComeFrom)
            OutputX = np.delete(dXBetweenBlobs,DirectionComeFrom)
            OutputY = np.delete(dYBetweenBlobs,DirectionComeFrom)
            OutputCoordinates = [OutputX,OutputY]
            ChangeInDirection = 0
        else:
            mean_alpha = np.mean(alpha)
            # Determine the alternative candidate (opposite direction)            
            if mean_alpha > 180:
                alternative_candidate = mean_alpha - 180
            else:
                alternative_candidate = mean_alpha + 180

            # Calculate width between blobs
            AlphaBetweenBlobs = np.diff(alpha)
            Radius =self.active_disc['r'][0] #np.mean(self.active_disc['r'])
            OutputWidth = 2*Radius*np.sin(0.5*AlphaBetweenBlobs)
            # Calculate coordinates between blobs
            dXBetweenBlobs = Radius*np.cos(0.5*AlphaBetweenBlobs)
            dYBetweenBlobs = Radius*np.sin(0.5*AlphaBetweenBlobs)

            OutputCoordinates = [dXBetweenBlobs,dYBetweenBlobs]

            # Determine which candidate is closest to current direction
            candidates = np.array([mean_alpha,alternative_candidate])
            # Distance (in degrees) between candidates and current direction)
            DirectionCandidates = np.abs(self.current_direction-candidates)
            # 3 and 360 degrees are very near. The maximal distance of one candidate is always 180 degrees
            DirectionCandidates[DirectionCandidates > 180] = 360 - DirectionCandidates[DirectionCandidates > 180]
            pick = np.argmin(DirectionCandidates)
            ChangeInDirection = DirectionCandidates[pick]
            OutputDirection =candidates[pick]


        return OutputDirection,ChangeInDirection,OutputWidth,OutputCoordinates