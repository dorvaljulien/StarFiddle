import time
import numpy as np
import cPickle as pk

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import Slider
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib.widgets import Button

class PlotAnimation:
    """
    Takes a list of PySnap and launch an interactive animation.
    """
    def __init__(self, Anim, DataControl, **kwargs):
        self.data = DataControl
        self.previous_n = 0
        self.Anim= Anim
    def timer_update(self):
        if self.Anim.n != self.previous_n:
            self.data.Update(self.Anim.n)
        self.previous_n = self.Anim.n
    def launch(self):
        self.data.Initialize()
        self.timer=self.data.fig.canvas.new_timer(interval=self.Anim.delay)
        args=[]
        # We tell the timer to call the update function every 100ms
        self.timer.add_callback(self.timer_update,*args)
        self.timer.start()



if __name__ is "__main__":
    from quick import *
    R =  ReadRun("/home/dorval/work/amuse/clump_finding/p10k_fragmentation/")
    R.Animation()
    A = R.Anim
    class Data():
        def __init__(self,):
            self.X = [ np.random.random(20) for i in range(len(R))]
            self.Y = [ np.random.random(20) for i in range(len(R))]
        def Initialize(self):
            X, Y = self.X[0], self.Y[0]
            self.fig = plt.figure()
            self.ax = self.fig.add_subplot(111)
            self.line, = plt.plot(X, Y, "b")
            self.canvas=self.ax.figure.canvas
        def Update(self, nx, ny):
            X, Y = self.X[n], self.Y[n]
            self.line.set_data(X, Y)
            self.canvas.draw()
    D= Data()
    P = PlotAnimation(A,D)
    P.launch()
    plt.show()
