import numpy as np
import cPickle as pk
import time

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib.widgets import Button
from matplotlib.widgets import Slider

class Animation:
    """
    Takes a list of PySnap and launch an interactive animation.

    SliderTime: Boolean, if True, the slider shows times. If false,
                it shows snap number.
    """
    def __init__(self,snaplist, SliderTime=True, **kwargs):
        self.snaplist=snaplist
        self.nsnap=len(snaplist)
        self.SliderTime = SliderTime
        self.time=[]
        for S in snaplist:
            self.time.append(S.t)
        self.n=0
        self.delay= 20 # delay between frames, in ms
        self.symbol= "bo" 
        self.markersize= 2
        self.pause_switch=True # The run is paused by default
        self.BackgroundColor='white'
        self.kwargs=kwargs
    def create_frame(self):
        """Build the window, create the line object"""
        self.fig = plt.figure(figsize=(10,10))
        self.ax = self.fig.add_subplot(111, projection='3d',
                                       adjustable='box',
                                       axisbg=self.BackgroundColor)
        self.ax.set_aspect('equal')
        plt.tight_layout()
        self.ax.set_xlabel("x")
        self.ax.set_ylabel("y")
        self.ax.set_zlabel("z")
        X = self.snaplist[self.n].x
        Y = self.snaplist[self.n].y
        Z = self.snaplist[self.n].z
        max_range = np.array([X.max()-X.min(), 
                              Y.max()-Y.min(), 
                              Z.max()-Z.min()]).max() / 3.0
        #mean_x = X.mean(); mean_y = Y.mean(); mean_z = Z.mean()
        mean_x =0; mean_y = 0; mean_z = 0
        self.fig.subplots_adjust(bottom=0.06) # Adjustment to make room for slider
        self.line,=self.ax.plot(X, Y, Z, self.symbol, markersize=self.markersize,
                                **self.kwargs )
        self.canvas=self.ax.figure.canvas
        self.ax.set_xlim(mean_x - max_range, mean_x + max_range)
        self.ax.set_ylim(mean_y - max_range, mean_y + max_range)
        self.ax.set_zlim(mean_z - max_range, mean_z + max_range)
        ax_pauseB=plt.axes([0.04, 0.02, 0.06, 0.025])
        self.pauseB=Button(ax_pauseB,'Play')
        self.pauseB.on_clicked(self.Pause_button)
        slider_ax = plt.axes([0.18, 0.02, 0.73, 0.025])
        if self.SliderTime:
            self.slider_time = Slider(slider_ax, 
                                      "Time", 
                                      self.snaplist[0].t, 
                                      self.snaplist[self.nsnap-1].t, 
                                      valinit = self.snaplist[0].t, 
                                      color = '#AAAAAA')
            self.slider_time.on_changed(self.slider_time_update)
        else:     
            self.slider_snap = Slider(slider_ax, 
                                      "Nsnap", 
                                      0, 
                                      self.nsnap, 
                                      valinit = 0, 
                                      color = '#AAAAAA')
            self.slider_snap.on_changed(self.slider_snap_update)
    def update_line(self):
        """Using self.n as a reference, we update the line and the slider"""
        self.line.set_data(self.snaplist[self.n%self.nsnap].x,
                           self.snaplist[self.n%self.nsnap].y)
        self.line.set_3d_properties(self.snaplist[self.n%self.nsnap].z)
        if self.SliderTime:
            self.slider_time.set_val(self.snaplist[self.n%self.nsnap].t)        
        else:
            self.slider_snap.set_val(self.n%self.nsnap)   
        self.canvas.draw()
    def timer_update(self,line):
        if not self.pause_switch:
            self.n=self.n+1
            self.update_line()
    def stop_timer_on_close(self,event):
        self.timer.stop()
    def Pause_button(self,event):
        self.pause_switch = not self.pause_switch
    def slider_time_update(self,val):
        for i in range(self.nsnap-1):
            if self.time[i]<=val and self.time[i+1]>val:
                self.n=i
        if val==self.time[self.nsnap-1]:
            self.n=self.nsnap-1
        self.line.set_data(self.snaplist[self.n].x,
                           self.snaplist[self.n].y)
        self.line.set_3d_properties(self.snaplist[self.n].z)
        self.canvas.draw()
    def slider_snap_update(self,val):
        for i in range(self.nsnap):
            if i<=val and i+1>val:
                self.n=i
        self.line.set_data(self.snaplist[self.n].x,
                           self.snaplist[self.n].y)
        self.line.set_3d_properties(self.snaplist[self.n].z)
        self.canvas.draw()
    def launch(self):
        # We create the figure, axis, slider and button, and plot the initial state
        self.create_frame()
        # When the figure is closed, the timer is stopped
        self.fig.canvas.mpl_connect('close_event', self.stop_timer_on_close)
        # The timer is setup with the right delay
        self.timer=self.fig.canvas.new_timer(interval=self.delay)
        args=[self.line]
        # We tell the timer to call the update function at each "tick"
        self.timer.add_callback(self.timer_update,*args)
        self.timer.start()



def plot3d(x,y,z,*args,**kwargs):
    """
    Simple 3d plotting function. If a figure is provided through
    the keyword "fig", it is used, and no plt.show() will be performed
    inside the function. If an ax object is provided through "ax", it 
    will be used:
       plot3d( x, y, z, fig = figure12 )
       plot3d( x, y, z, ax = ax3 )
    Any other arguments or keyword arguments will be sent to the normal
    plot function:
       plot3d(x, y, z, "bo" )
       plot3d(x, y, z, linewidth=4 )
    """
    if not kwargs.has_key("ax"):
        if kwargs.has_key("fig"):
            fig =  kwargs.pop("fig")
            Show = False
        else:
            fig = plt.figure()
            Show = True
        ax=fig.add_subplot(111,projection="3d")
    else:
        ax = kwargs.pop("ax")
        Show=False
    ax.plot(x,y,z,*args,**kwargs)
    if Show:
        plt.show()
    return ax



class AddAnimation():
    """
    Hitch a ride on an existing Animation instance and add data to be displayed.
    ( Remember these are better used within the ReadRun class environment )    
    
    # snaplist1 and snaplist2 have to have the same number of snap
    Anim1 = Animation( snaplist1 )
    Anim1.launch()
    Anim2 = AddAnimation( Anim, snaplist2 )
    Anim2.launch()

    The time on the slider is the one of snaplist1, the time contained 
    in the PySnaps of snaplist2 is not taken into account.
    """
    
    def __init__(self,OldAnim, snaplist, **kwargs):
        if len(snaplist) is not OldAnim.nsnap:
            raise Exception("Number of snap doesn't fit with main animation.")
        self.n = OldAnim.n
        self.nsnap = len(snaplist)
        self.ax = OldAnim.ax
        self.delay = OldAnim.delay
        self.previous_n = OldAnim.n
        self.OldAnim = OldAnim
        self.snaplist = snaplist
        self.symbol = "bo"
        self.markersize = 2
        self.kwargs = kwargs
        self.fig = OldAnim.fig
    def create_frame(self):
        X = self.snaplist[self.n].x
        Y = self.snaplist[self.n].y
        Z = self.snaplist[self.n].z
        max_range = np.array([X.max()-X.min(), 
                              Y.max()-Y.min(), 
                              Z.max()-Z.min()]).max() / 3.0
        mean_x = X.mean(); mean_y = Y.mean(); mean_z = Z.mean()
        self.line,=self.ax.plot(X, Y, Z, self.symbol, markersize=self.markersize,
                                **self.kwargs )
    def update_line(self,new_n):
        self.line.set_data(self.snaplist[new_n%self.nsnap].x,
                           self.snaplist[new_n%self.nsnap].y)
        self.line.set_3d_properties(self.snaplist[new_n%self.nsnap].z)
    def timer_update(self,line):
        if self.OldAnim.n != self.previous_n:
            self.update_line(self.OldAnim.n)
    def launch(self):
        self.create_frame()
        self.timer=self.fig.canvas.new_timer(interval=self.delay)
        args=[self.line]
        # We tell the timer to call the update function every 100ms
        self.timer.add_callback(self.timer_update,*args)
        self.timer.start()


















