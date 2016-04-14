import numpy as np
import cPickle as pk
import time
from copy import deepcopy

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib.widgets import Button
from matplotlib.widgets import Slider

class SlidersAnimation:
    """
    Launch an animation with an arbitrary number of sliders.
    See the source code for Examplemultipleslidersdata to see
    what Data should be.
    """
    def __init__(self, Data, **kwargs ):
        self.D = Data
        self.kwargs = kwargs
        self.Ns =  self.D.dim
        self.ns = np.zeros(self.Ns)
        self.previous_ns = deepcopy(self.ns)
        self.vals = [ V[0] for V in self.D.Variables ]
    def create_frame(self,Ns):
        fig, ax = self.D.fig, self.D.ax
        fig.subplots_adjust(bottom=0.06+Ns*0.04)#, left=0.1)
        self.canvas=self.D.ax.figure.canvas
        self.sliders = []
        for i in range(Ns):
            slider_ax = plt.axes([0.15, (Ns-i)*0.04, 0.73, 0.025])
            self.sliders.append(Slider(slider_ax, 
                                      self.D.Names[i], 
                                      self.D.Variables[i][0], 
                                      self.D.Variables[i][-1], 
                                      self.D.Variables[i][0], 
                                      color = '#AAAAAA'))
            self.sliders[-1].on_changed(self.SliderUpdate)
    def SliderUpdate(self,val):
        for i in range(self.D.dim):
            ind = np.nonzero(self.D.Variables[i] <= self.sliders[i].val)[0]
            if len(ind) == 0:
                self.ns[i] = 0
            else:
                self.ns[i] = np.max(ind)
        if sum(self.ns!=self.previous_ns) is not 0:
            self.D.Update(self.ns)   
            self.D.canvas.draw()
            self.previous_ns = deepcopy(self.ns)
    def launch(self):
        self.D.Initialize()
        self.create_frame(self.D.dim)


class ExampleMultipleSlidersData():
    """
    This is what Data should be when fed to SlidersAnimation.
    This instance is fed to SlidersAnimation, that will call its
    Initialize and Update methods (so don't change method names !
    Or do, and suffer the consequences, you fool.)

    This example is tested after this definition.

    In this example, I'm manipulating a sine wave by changing its
    frequency, phase and amplitude. The names of the parameters
    are provided in self.Names, and the ranges in self.Variables.
    self.X is just the base array on which I build the sine wave.

    In the Initialize method, I just build the figure, making room
    for enough sliders thanks to self.dim. Afterwards, I plot a basic 
    sine wave and store the line object in self.line. I also store
    the canvas in self.canvas.

    These are modified in the Update method. The ns array is the
    array containing the indices of the values at which sliders 
    currently are, we extract them and name them for clarity.
    Then we modify the line object and update the canvas.

    This is very flexible, anything that can be plotted and updated
    later can be put into this form ! Have fun !
    """
    def __init__(self):
        # You can change this
        self.Names = ["Frequency", "Phase", "Amplitude"]
        self.Variables = [ np.linspace(0.5,10,100),
                           np.linspace(0.,15,100),
                           np.linspace(0.1,1.,20) ]
        # Important, don't touch this
        self.dim = len(self.Names)
        # This is particular to this example
        self.X = 0.1 * np.arange(1000) 

    def Initialize(self):
        # Probably shouldn't change this
        #  ( or maybe the figure size )
        self.fig = plt.figure(figsize=(10,10+self.dim*0.06))
        self.ax = self.fig.add_subplot(111)
        self.canvas=self.ax.figure.canvas
        # This is particular to this example
        X = self.X
        Y = np.sin( X )
        # Creation of an updatable object
        self.line, = plt.plot(X, Y)

    def Update(self, ns):
        # This whole method is very flexible: you just get
        # the ns array, retrieve the corresponding values 
        # in self.Variables and update your object the way you 
        # like.
        nfreq,nphase,namp = int(ns[0]), int(ns[1]), int(ns[2])
        freq = self.Variables[0][nfreq]
        phase = self.Variables[1][nphase]
        amp = self.Variables[2][namp]
        X, Y = self.X, amp * np.sin( freq*self.X + phase )
        self.line.set_data(X, Y)
        self.canvas.draw()

if __name__ == "__main__":
    # Testing the example
    D = ExampleMultipleSlidersData()
    A = SlidersAnimation(D)
    A.launch()

