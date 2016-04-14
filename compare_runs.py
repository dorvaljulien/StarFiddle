import numpy as np
import cPickle as pk
import matplotlib.pyplot as plt
import time

from matplotlib.widgets import Slider
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib.widgets import Button

from quick import *


class CompareAnimation:
    """
    Launch two parallel animations of runs, so the user can 
    easily compare the structures. Example:
        R1 = ReadRun("fake/path/run_1/")
        R2 = ReadRun("fake/path/run_2/")
        new_anim=CompareAnimation(R1.S,R2.S)
        new_anim.launch()
        plt.show()
    """
    def __init__(self,snaplist1,snaplist2, symbol="bo", dt = None, markersize=2, **kwargs):
        if snaplist1[0].t < snaplist2[0].t:
            S = deepcopy(snaplist1[0])
            S.t = snaplist2[0].t
            snaplist1.reverse()
            snaplist1.append(S)
            snaplist1.reverse()
        if snaplist2[0].t < snaplist1[0].t:
            S = deepcopy(snaplist2[0])
            S.t = snaplist1[0].t
            snaplist2.reverse()
            snaplist2.append(S)
            snaplist2.reverse()
        self.snaplists = [ snaplist1, snaplist2 ]
        self.times = [ [s.t for s in snapl] for snapl in self.snaplists ]
        #self.nsnap=len(snaplist)
        self.n = [ 0, 0 ]
        self.t=0.
        self.tmax = max(self.snaplists[0][-1].t, self.snaplists[1][-1].t)
        self.dt = self.tmax/200. if dt is None else dt
        self.delay=1
        self.symbol= symbol
        self.markersize= markersize
        self.pause_switch=True
        self.BackgroundColor='white'
        self.kwargs=kwargs
    def create_frame(self):
        self.fig = plt.figure(figsize=(15,10))
        print "fig created"
        self.ax = []
        for i in range(2):
            ax = self.fig.add_subplot(121+i, projection='3d',
                                       adjustable='box',
                                       axisbg=self.BackgroundColor) 
            ax.set_aspect('equal')
            plt.tight_layout()
            ax.set_xlabel("x (pc)")
            ax.set_ylabel ("y (pc)")
            ax.set_zlabel ("z (pc)")
            self.ax.append(ax)
        X,Y,Z = [],[],[]
        for snapl in self.snaplists:
            X.append( snapl[0].x )
            Y.append( snapl[0].y )
            Z.append( snapl[0].z )
        max_range = np.array([X[0].max()-X[0].min(), 
                              Y[0].max()-Y[0].min(), 
                              Z[0].max()-Z[0].min()]).max() / 3.0
        mean_x = X[0].mean(); mean_y = Y[0].mean(); mean_z = Z[0].mean()
        self.fig.subplots_adjust(bottom=0.06)#, left=0.1)
        self.line, self.canvas = [],[]
        for (x,y,z,ax) in zip(X,Y,Z,self.ax):
            self.line.append(ax.plot(x, y, z, self.symbol, markersize=self.markersize,
                                    **self.kwargs )[0])
            self.canvas.append(ax.figure.canvas)
            ax.set_xlim(mean_x - max_range, mean_x + max_range)
            ax.set_ylim(mean_y - max_range, mean_y + max_range)
            ax.set_zlim(mean_z - max_range, mean_z + max_range)
        ax_pauseB=plt.axes([0.04, 0.02, 0.06, 0.025])
        self.pauseB=Button(ax_pauseB,'Play')
        self.pauseB.on_clicked(self.Pause_button)
        slider_ax = plt.axes([0.18, 0.02, 0.73, 0.025])
        self.slider_time = Slider(slider_ax, 
                                  "Time", 
                                  self.snaplists[0][0].t, 
                                  self.snaplists[0][-1].t, 
                                  valinit = self.snaplists[0][0].t, 
                                  color = '#AAAAAA')
        self.slider_time.on_changed(self.slider_time_update)
    def update_lines(self):
        nsnaps = [0,0]        
        for j,snaplist,time in zip(range(2),self.snaplists,self.times):
            ind = np.nonzero(time < self.t)[0]
            nsnaps[j] = max(ind) if len(ind)!=0 else 0
        for (line,snaplist,n) in zip(self.line, self.snaplists, nsnaps):
            line.set_data(snaplist[n].x, snaplist[n].y)
            line.set_3d_properties(snaplist[n].z)
        for canv in self.canvas:
            canv.draw()
    def timer_update(self,lines):
        if not self.pause_switch:
            self.t = self.t+self.dt
            if self.t > self.tmax:
                self.t = self.t - self.tmax
            self.update_lines()
            self.slider_time.set_val(self.t)        
    def Pause_button(self,event):
        self.pause_switch = not self.pause_switch
    def slider_time_update(self,val):
        self.t = val
        self.update_lines()
    def launch(self):
        self.create_frame()
        self.timer=self.fig.canvas.new_timer(interval=self.delay)
        args=[self.line]
        # We tell the timer to call the update function every 100ms
        self.timer.add_callback(self.timer_update,*args)
        self.timer.start()


