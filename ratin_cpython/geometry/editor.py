"""
.. module:: editor
    :platform: Windows
    :synopsis: 
    

The editor class supplies a tool to 
* draw river branch on empty canvas
* edit river branches
* generate grid
* save grid   
matplotlib canvas.  

Dependencies:

    * numpy (1.7+)
    * SciPy (0.12+)
    * Matplotlib
    

"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.artist import Artist
from matplotlib.mlab import dist_point_to_segment
from scipy import interpolate

class Editor:
    """
    ---------------------------------------------
    Launches the river editor in given figure/axes. Use:
    
    from rattool.geometry import rivereditor
    fig,ax = plt.subplots(1,1)
    rv = rivereditor.RiverEditor(fig,ax)

    ---------------------------------------------    
    Default key-bindings
    
      'g' grab and move selected vertex
      
      'v' begin drawing of new line
      
      '.' increase width of selected vertex
      
      'c' reset view
      
      'r' box zoom
      
      't' toggle grid visibility
      
      'm' generate grid

      'd' delete selected vertex

      'e' extrude vertex (only first and last of a line)
      
      ---------------------------------------------

    """

    def __init__(self,fig=None,ax=None):
        print 'initializing'
        '''
        // Default settings
        '''
        self.default_width = 0.002  # Initial width of the river
        self.grid_n = 8           # Number of n-lines for grid generation
        self.spline_factor = 10    # Number vertices in spline compared to lines (= number of m-lines)
        self.epsilon = 15  # max pixel distance to count as a vertex hit
        self.figure_color = [0.9,0.9,0.9,1]
        self.riverwidth_modifier = 1.1

        '''
        // Key bindings
        '''
        # RAT key bindings
        self.key_bindings = dict()
        self.key_bindings['add_vert'] = 'v'
        self.key_bindings['move_vert'] = 'g'
        self.key_bindings['delete_vert'] = 'd'
        self.key_bindings['insert_vert'] = 'i'
        self.key_bindings['multiselect'] = 'shift'
        self.key_bindings['increase_width'] = '.'
        self.key_bindings['generate_grid'] = 'm'
        self.key_bindings['extrude_vert'] = 'e'
        
        # Overwrite matplotlib keybindings
        mpl.rcParams['keymap.grid'] = 'h'
        mpl.rcParams['keymap.home'] = 'c'
        mpl.rcParams['keymap.zoom'] = 'r'

        '''
        // Initialize figure (do not edit!)
        '''
        if fig is None or ax is None:
            self.launch_standalone()
        else:
            # Import figure to class structure 
            self.fig = fig
            self.ax = ax

        # Initialization of internal parameters
        self._ind = [0]                 
        self._multiselect = False
        self._selected_action = None
        self._draw_geometry = False
        self._plot_grid = False
        self._panview = False
        self._initiate_panning = True
        
        self.fig.set_facecolor(self.figure_color)
        self.canvas = self.fig.canvas
        

        # Connect callbacks to methods
        self.canvas.mpl_connect('draw_event', self.draw_callback)
        self.canvas.mpl_connect('button_press_event', self.button_press_callback)
        self.canvas.mpl_connect('key_press_event', self.key_press_callback)
        self.canvas.mpl_connect('button_release_event', self.button_release_callback)
        self.canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)
        self.canvas.mpl_connect('key_release_event', self.key_release_callback)
        
    '''
    // -----------------------------------------------------------------------
    // Callbacks
    // -----------------------------------------------------------------------
    '''
    def draw_callback(self, event):
        if self._draw_geometry:
            self.background = self.canvas.copy_from_bbox(self.ax.bbox)
            self.ax.draw_artist(self.line)
            self.ax.draw_artist(self.spline)
            self.ax.draw_artist(self.left_spline)
            self.ax.draw_artist(self.selected_vert)
            self.draw_grid()
            
    def button_press_callback(self, event):
        '''
        Callback for mouse button press (1 = left, 2 = middle, 3 = right)
        '''
        if event.inaxes==None: return
        if event.button == 1: 
            self._selected_action = None

            self.selected_action = 'm'
        elif event.button == 2:
            return
            #self._panview = True
            #self.initiate_vert_movement(event)
            #self.pan_view(event)
        # Right mouse click
        elif event.button == 3:
            if self._multiselect:
                ind = self.get_ind_under_point(event)
                if ind is not None and self._ind[0] is not None:
                    self._ind = np.append(self._ind,self.get_ind_under_point(event))
                elif ind is not None and self._ind[0] is None:
                    self._ind = [self.get_ind_under_point(event)]
                self.select_vert()
            else:
                
                self._ind = [self.get_ind_under_point(event)]
                self.select_vert()            

    def button_release_callback(self, event):
        '''
        whenever a mouse button is released (not needed for now)
        '''
        self._panview = False
        self._initiate_panning = True
        return

    def key_press_callback(self, event):
        '''
        This methods binds the key to a method
        '''
        print event.key
        if not event.inaxes:  # If key pressed outside the canvas, do nothing
            return
        if event.key=='t':   self.toggle_grid_visibility()
        elif event.key==self.key_bindings['add_vert']: self.draw_line(event)
        elif event.key==self.key_bindings['move_vert']: self.grab_vert(event)   
        elif event.key==self.key_bindings['delete_vert']: self.delete_vert()
        elif event.key==self.key_bindings['insert_vert']: self._selected_action = 'i'
        elif event.key==self.key_bindings['extrude_vert']: self.extrude_vert(event)
        elif event.key == self.key_bindings['multiselect']: self._multiselect=True
        elif event.key == self.key_bindings['increase_width']: self.increase_width(event)
        elif event.key == self.key_bindings['generate_grid']: self.generate_grid()        
        else: self.selected_action = None
     
    def key_release_callback(self,event):
        if event.key == self.key_bindings['multiselect']: self._multiselect = False
            
    def motion_notify_callback(self, event):
        '''
        On mouse movement
        '''
        if self._ind[0] is None: return
        if event.inaxes is None: return
        if self._selected_action is 'm': self.move_vert(event)
        elif self._panview: self.pan_view(event) 
           

            
    '''
    // -----------------------------------------------------------------------
    // Line and grid edit methods
    // -----------------------------------------------------------------------
    '''        
    def draw_line(self,event):
        '''
        Method to begin drawing a new line
        '''
        x,y = event.xdata, event.ydata
        
        xs = np.array([x,x+0.01])
        ys = np.array([y,y+0.01])
        self._draw_geometry = True
        self._ind = [1]
        # Add initial line
        self.river_widths = list()
        self.river_widths.append(self.default_width*np.ones(len(xs)))
        self.line = Line2D(xs, ys, color = [0.6,0.6,0.6],linestyle='--',marker='o', markerfacecolor='r', animated=True)
        self.selected_vert = Line2D([xs[self._ind],xs[self._ind]+0.000001],
                                    [ys[self._ind],ys[self._ind]+0.000001],
                                    color = [0.6,0.6,0.6],linestyle='none',
                                    marker='o', markerfacecolor='c', animated=True)
        
        # Create smoothed spline through points
        self.create_spline(xs,ys)

        # Plot line
        self.ax.add_line(self.line)
        self.ax.add_line(self.selected_vert)
        self._selected_action = 'm'
        
    def toggle_grid_visibility(self):
        self._plot_grid = not self._plot_grid
        self.canvas.draw()
        
    def draw_grid(self):
        '''
        Method to draw (as in plot) the grid
        '''
        if self._plot_grid:
            for line in self.grid[0]['n']:
                self.ax.draw_artist(line)
            for line in self.grid[0]['m']:
                self.ax.draw_artist(line)
        self.canvas.blit(self.ax.bbox)
        
    def subdivide(self,x,factor):
        '''
        Subdivide vector x with factor 
        '''
        
        xnew = [float(x[0])]
        xnold = xnew
        for xn in x[1:]:
            stepsize = (xn-xnold)/factor
            for i in range(0,factor-1):
                xnew = np.append(xnew,xnew[-1]+stepsize)
            xnew = np.append(xnew,xn)
            xnold = float(xn)
        return xnew
        
    def create_spline(self,xs,ys):
        self.spline = Line2D(np.zeros(len(xs)),np.zeros(len(xs)),
                             color='r',animated=True)
        self.left_spline= Line2D(np.zeros(len(xs)),np.zeros(len(xs)),
                                 color='g',animated=True) 
        self.update_spline(xs,ys)   
        self.ax.add_line(self.spline)
        self.ax.add_line(self.left_spline)   
        
    def generate_grid(self):
        # N lines
        self.grid = []
        grid = {'n':[],'m':[]}
        x = self.spline.get_xdata()
        y = self.spline.get_ydata()
        xm = []
        ym= []
        # Tangent using complex j for x,y coordinates
        theta = np.angle(np.diff(x)+np.diff(y)*1j)        
        theta = np.concatenate((np.array([theta[0]]),theta))
        
        perpendicular_theta = theta + 0.5*np.pi
      
        
        widths = self.spline_widths
        for i in range(0,self.grid_n+1):
            n_widths = i*(widths/self.grid_n)
            xl = np.array(x)+n_widths*np.array(np.cos(perpendicular_theta))
            yl = np.array(y)+n_widths*np.array(np.sin(perpendicular_theta))
            
            for j in range(0,len(xl)):
                try:
                    xm[j].append(xl[j])
                    ym[j].append(yl[j])
                except IndexError:
                    xm.append([xl[j]])
                    ym.append([yl[j]])
                                        
                   
            grid['n'].append(Line2D(xl,yl, color = [0,0,0],linestyle='-',animated=True))
            self.ax.add_line(grid['n'][i])
        
        
        for j in range(0,len(x)):
            grid['m'].append(Line2D(xm[j],ym[j],color = [0.3,0.3,0.3], animated=True))
            self.ax.add_line(grid['m'][j])
        
        self.grid.append(grid)
        self._plot_grid = True
        self.canvas.draw()
        
    def update_spline(self,xs,ys):
        if not(len(xs)>2): return
        # Calculate spline's x,y coordinates
        tck,u = interpolate.splprep([xs,ys],s=0,k=2,nest = -1)
        out = interpolate.splev(np.linspace(0,1,self.spline_factor*len(xs)),tck,der=0)
        
        # Calculate spline's w,s coordinates
        s  = np.sqrt(np.diff(xs)**2+np.diff(ys)**2)
        s = np.concatenate((np.array([0]),s))
        s = np.cumsum(s)
        w = self.river_widths[0] 

        
        snew = np.sqrt(np.diff(out[0])**2+np.diff(out[1])**2)
        snew = np.concatenate((np.array([0]),snew))
        snew = np.cumsum(snew)
  
        wnew = np.interp(snew,s,w)
        
        self.spline_widths = wnew

        
        self.update_spline_widths(out[0],out[1],wnew)
        self.spline.set_data(out[0],out[1])
        self.canvas.draw()
    
    def update_spline_widths(self,x,y,widths):
        # Tangent using complex j for x,y coordinates
        alpha = np.angle(np.diff(x)+np.diff(y)*1j)        
        alpha = np.concatenate((np.array([alpha[0]]),alpha))
        theta = alpha
        perpendicular_theta = theta + 0.5*np.pi
             
        
        xl = np.array(x)+widths*np.array(np.cos(perpendicular_theta))
        yl = np.array(y)+widths*np.array(np.sin(perpendicular_theta))
        self.left_spline.set_data(xl,yl)
        self.canvas.draw()
        
        
    def pan_view(self,event):
        if self._initiate_panning:
            self._lastposition = [event.xdata,event.ydata]
            self._initiate_panning = False
            
        dx = -(event.xdata - self._lastposition[0])*0.85
        dy = -(event.ydata - self._lastposition[1])*0.85
        current_ax_limits = self.ax.axis()
        self.ax.axis([current_ax_limits[0]+ dx,
                        current_ax_limits[1]+dx,
                        current_ax_limits[2]+dy,
                        current_ax_limits[3]+dy])
                        
        self._lastposition = [event.xdata,event.ydata]                
        self.canvas.draw()
        
    '''
    // -----------------------------------------------------------------------
    // Vertex functions
    // -----------------------------------------------------------------------
    '''
    def grab_vert(self,event):
        self._selected_action = 'm'
        self.initiate_vert_movement(event)   
        
    def select_vert(self):
        '''
        This method highlights the selected vertex (does not actually 'select'
        it. See get_ind_under_point for this)
        '''
        if not(self._ind[0] is None):
            self.selected_vert.set_visible(True)
            x = self.line.get_xdata()
            y = self.line.get_ydata()
            if len(self._ind) == 1:      
                self.selected_vert.set_data([x[self._ind[0]],x[self._ind[0]]+0.00000001],
                                            [y[self._ind[0]],y[self._ind[0]]+0.00000001])
            else:
                self.selected_vert.set_data(x[self._ind],y[self._ind])
                                            
            self.canvas.draw()
        else:
            self.selected_vert.set_visible(False)
            self.canvas.draw()

    def get_ind_under_point(self, event):
        '''
        Get the index of the vertex under point if within epsilon tolerance
        '''
        # display coords
        xy = np.asarray([self.line.get_xdata(), self.line.get_ydata()])
        xy = xy.transpose()
        xyt = self.line.get_transform().transform(xy)
        xt, yt = xyt[:, 0], xyt[:, 1]
        d = np.sqrt((xt-event.x)**2 + (yt-event.y)**2)
        indseq = np.nonzero(np.equal(d, np.amin(d)))[0]
        ind = indseq[0]

        if d[ind]>=self.epsilon:
            ind = None

        return ind
   
    def move_vert(self,event):
        x,y = event.xdata, event.ydata
        xline = self.line.get_xdata()
        yline = self.line.get_ydata()
        if len(self._ind) > 1:
            # Multiple vertices selected, do not drag to mouse but move relative
            # to mouse. 
            xline[self._ind] = x-self._distance_at_initiation[0]
            yline[self._ind] = y-self._distance_at_initiation[1]
        else:
            try:
                yline[self._ind] = y
                xline[self._ind] = x
            except:
                print self._ind
        self.line.set_data(xline,yline)
        self.update_spline(xline,yline)
        self.select_vert()

        self.canvas.restore_region(self.background)

        self.ax.draw_artist(self.line)
        self.ax.draw_artist(self.spline)
        self.ax.draw_artist(self.left_spline)
        self.ax.draw_artist(self.selected_vert)
        self.draw_grid()
        self.canvas.blit(self.ax.bbox)
        
    def delete_vert(self):
        '''
        Delete the selected vertex or vertices. (x,y and width)
        '''
        if self._ind[0] is not None:
            for ind in self._ind:
                x = [tup for i,tup in enumerate(self.line.get_xdata()) if i!=ind]
                y = [tup for i,tup in enumerate(self.line.get_ydata()) if i!=ind]
                w = [tup for i,tup in enumerate(self.river_widths[0]) if i!=ind]    
                self.river_widths[0] = w                 
                self.line.set_data(np.array(x),np.array(y)) 
            
            self._ind = [None]
            self.update_spline(x,y)
            self.select_vert()                   
            self.canvas.draw()
    
    def extrude_vert(self,event):
        '''
        This methods 'extrudes' from the selected vertices. 
        
        >> only from the first and last vertex for now
        
        '''
        if len(self._ind) is not 1: return
        
        x = self.line.get_xdata()
        y = self.line.get_ydata()

        if self._ind[0] == len(x)-1:
            x = np.append(x,event.xdata+1e-8)
            y = np.append(y,event.ydata+1e-8)
            self.river_widths[0] = np.append(self.river_widths[0],self.river_widths[0][-1])
            self._ind[0] += 1
        elif self._ind[0] == 0:
            x = np.append(x[0]-0.1,x)
            y = np.append(y[0]-0.1,y)
            self.river_widths[0] = np.append(self.river_widths[0][0],self.river_widths[0])
           
        self.line.set_xdata(x)
        self.line.set_ydata(y)
        self.update_spline(x,y)
        
        self.select_vert()
        self._selected_action = 'm'
        self.canvas.draw()
            
    def initiate_vert_movement(self,event):
        '''
        At start of movement (after pressing 'grab vert' key, the difference
        between the mouse and the vertices is calculated. Vertices are moved 
        according to their relative position from the mouse pointer.)
        '''
        if self._ind[0] is None: return
        xv,yv = self.line.get_xdata(),self.line.get_ydata()
        xm,ym = event.xdata,event.ydata
        self._distance_at_initiation = [xm-xv[self._ind],ym-yv[self._ind]]
        
    def increase_width(self,event):
        '''
        Increase the width of the selected vertex. 
        '''
        self.river_widths[0][self._ind] = self.river_widths[0][self._ind]*\
                                          self.riverwidth_modifier
        self.update_spline(self.line.get_xdata(),self.line.get_ydata())
        self.canvas.draw()

    def launch_standalone(self):
        self.fig,self.ax = plt.subplots(1,1)
        
        

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    fig,ax = plt.subplots(1,1) 
    p = Editor(fig,ax)
    plt.show()
