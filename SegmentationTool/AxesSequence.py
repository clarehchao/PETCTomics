# -*- coding: utf-8 -*-
"""
Created on 5/2/16 3:02 PM

@author: shuang Shih-ying Huang

source code: http://stackoverflow.com/questions/13443474/matplotlib-sequence-of-figures-in-the-same-window

"""

import matplotlib.pyplot as plt

class AxesSequence(object):
    """Creates a series of axes in a figure where only one is displayed at any
    given time. Which plot is displayed is controlled by the arrow keys."""
    def __init__(self):
        self.fig = plt.figure()
        self.axes = []
        self._i = 0 # Currently displayed axes index
        self._n = 0 # Last created axes index
        self._first_slice = 0
        self._second_slice = 0
        self.fig.canvas.mpl_connect('key_press_event', self.on_keypress)

    @property
    def last_display_index(self):
        """ """
        return self._i

    @property
    def selected_slices(self):
        """ """
        return [self._first_slice, self._second_slice]

    @property
    def last_createax_index(self):
        """ """
        return self._n

    @property
    def GetFig(self):
        """ """
        return self.fig

    @property
    def GetAxes(self):
        """ """
        return self.axes


    def __iter__(self):
        while True:
            yield self.new()

    def new(self):
        # The label needs to be specified so that a new axes will be created
        # instead of "add_axes" just returning the original one.
        ax = self.fig.add_axes([0.15, 0.1, 0.8, 0.8],
                               visible=False, label=self._n)
        self._n += 1
        self.axes.append(ax)
        return ax

    def on_keypress(self, event):
        if event.key == 'right':
            # print 'right key is pressed'
            self.next_plot()
        elif event.key == 'left':
            # print 'left key is pressed'
            self.prev_plot()
        elif event.key == '1':
            print 'detect 1'
            self.set_slice(1)
        elif event.key == '2':
            print 'detect 2'
            self.set_slice(2)
        else:
            return
        self.fig.canvas.draw()

    def next_plot(self):
        if self._i < len(self.axes):
            self.axes[self._i].set_visible(False)
            self.axes[self._i+1].set_visible(True)
            self._i += 1

    def prev_plot(self):
        if self._i > 0:
            self.axes[self._i].set_visible(False)
            self.axes[self._i-1].set_visible(True)
            self._i -= 1
    def set_slice(self,islice):
        if islice == 1:
            self._first_slice = self._i
            print 'set slice to {}'.format(self._first_slice )
        elif islice == 2:
            self._second_slice = self._i
            print 'set slice to {}'.format(self._second_slice )

    def show(self):
        self.axes[0].set_visible(True)
        plt.show()