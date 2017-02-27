import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import time

class Annotate(object):
    def __init__(self):
        self.ax = plt.gca()
        # self.ax = ax
        # self.fig = fig
        self.rect = Rectangle((0,0), 1, 1, facecolor='None', edgecolor='green')
        self.x0 = None
        self.y0 = None
        self.x1 = None
        self.y1 = None
        self.is_pressed = False
        self.ax.add_patch(self.rect)
        # self.fig.canvas.mpl_connect('button_press_event',self.on_press)
        # self.fig.canvas.mpl_connect('button_release_event', self.on_release)
        # self.fig.canvas.mpl_connect('motion_notify_event', self.on_motion)
        # self.fig.canvas.draw()
        
        self.ax.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.ax.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.ax.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.ax.figure.canvas.draw()
    def on_press(self, event):
        print 'press'
        self.is_pressed = True
        # if self.is_pressed is True:
        #     print 'on_press: press is true'
        self.x0 = event.xdata
        self.y0 = event.ydata
        print self.x0, self.y0
        #self.rect.set_width(self.x1 - self.x0)
        #self.rect.set_height(self.y1 - self.y0)
        #self.rect.set_xy((self.x0, self.y0))
        #self.rect.set_linestyle('dashed')
        #self.fig.canvas.draw()
        #self.ax.figure.canvas.draw()
    def on_motion(self,event):
        print 'motion'
        if self.is_pressed is False:
            # print 'on_motion: press is False'
            return
        self.x1 = event.xdata
        self.y1 = event.ydata
        self.rect.set_width(self.x1 - self.x0)
        self.rect.set_height(self.y1 - self.y0)
        self.rect.set_xy((self.x0, self.y0))
        self.rect.set_linestyle('dashed')
        #self.fig.canvas.draw()
        self.ax.figure.canvas.draw()
        print 'on_motion:after draw figure'
    def on_release(self, event):
        self.is_pressed = False
        self.x1 = event.xdata
        self.y1 = event.ydata
        print 'on_release: {},{},{},{}'.format(self.x0,self.y0,self.x1,self.y1)
        self.rect.set_width(self.x1 - self.x0)
        self.rect.set_height(self.y1 - self.y0)
        self.rect.set_xy((self.x0, self.y0))
        self.rect.set_linestyle('solid')
        #self.fig.canvas.draw()
        self.ax.figure.canvas.draw()
        # print 'on_release: after draw figure'
        # time.sleep(3)
        # print 'after 3-sec pause'
        # plt.close()
        #print self.x0,self.y0,self.x1,self.y1
        return [self.x0,self.x1,self.y0,self.y1]
