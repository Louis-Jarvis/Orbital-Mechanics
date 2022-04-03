from re import M
from tkinter import W, Y
from turtle import shape
from scipy.integrate import solve_ivp
import numpy as np
from typing import Tuple
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from collections import deque

from sympy import Interval

#TODO:
# equations of motion in:
#    - radial coordinates
#    - polar coordinates
#    - make this a class
# solver function/code?
# set tolerance

# fX = np.array([
#     [0, 0, 0, 1, 0, 0],
#     [0, 0, 0, 0, 1, 0],
#     [0, 0, 0, 0, 0, 1],
#     [1, 0, 0, 0, 0, 0],
#     [0, 1, 0, 0, 0, 0],
#     [0, 0, 1, 0, 0, 0]
# ])

def magnitude(vector: np.ndarray) -> float: 
    return np.sqrt(sum([x**2 for x in vector])) ## TODO use np.dot ??


## TODO ? make this a class???? maybe?
def ODE_2_body_system(t,X,mu,dof = 6):
    '''
    dX/dt = f(t,X) = f(X) where f(X) is a system of NxN nonlinear equations
    dX/dt = [  dxdt,  dydt,   dzdt,   d^2(x)/dt^2,   d^2(y)/dt^2,   d^2(z)/dt^2 ]
    and X = [ x y z dx/dt dy/dt dz/dt]
    '''
    A = np.eye(dof,dtype='float')[[3,4,5,0,1,2],]
    i,j = zip(*[(3,0), (4,1), (5,2)])
    A[i,j] = mu/( magnitude(X[:3])**3 )
    #np.fill_diagonal(fX,1)
    return np.dot(A, X) ## equates to dX/dt

#_f = lambda M,i,j: M[i,j] = 1
#map( , 0:6, [3, 4, 5, 0, 1, 2] )

def solve2body(ode: callable, ic: np.ndarray, ti: float, tf : float) -> Tuple[np.ndarray, np.ndarray]:
    sol = solve_ivp(fun=ode, t_span=(ti,tf), y0=ic, args = (1e6, 6))
    return sol.t, sol.y[0:3],sol.y[3:6]

def calcEnergy(r,v,mu):
    Ek = 0.5 * sum(np.dot(v,v)) ## (1/2)*m*(v^2 + )
    Ep = -mu/r ## use magnitude function??
    return Ek,Ep

##X0 = np.array([1,1,1,0,0,0])
##t,r,v = solve2body( ode = ODE_2_body_system, ic = X0, ti = 0, tf = 20 )
##x,y,z = r

class twoBody():
    def __init__(self, ICs: np.ndarray, ti: float, tf: float ,mu: float) -> None:
        ## solve system, initialise params
        self.t,r,v= solve2body(ode = ODE_2_body_system, ic = ICs, ti = ti, tf = tf)
        
        ## position vector components: xi, yj, zk
        self.x, self.y, self.z = r

        ## velocity vector components: ui, vj, wk
        self.u, self.v, self.w = v

        return
    def Re_solve(self):
        pass

    def plot_position_vs_time(self, vector = 'x', textsize =  8):
        _,ax = plt.subplots()

        if vector == 'xyz':
            ax.plot(self.t/3600, self.x, label = 'x') ## TODO use map with this?
            ax.plot(self.t/3600, self.y, label = 'y')
            ax.plot(self.t/3600, self.z, label = 'z')
        elif ( ('x'== vector) | ('y'==vector) | ('z'==vector) ):
            ax.plot(self.t/3600, getattr(self,vector), label = vector)

        ax.set_title('{:s} vector vs time'.format(vector), fontsize = textsize*2, fontweight = 'bold')
        ax.set_xlabel('Time (hours)', fontsize = textsize)
        ax.set_ylabel('{:s} component of position (m)'.format(vector), fontsize = textsize)
        ax.ticklabel_format(useOffset=False,style= 'plain')
        return ax

    def plot_velocity_vs_time(self, vector = 'u' ):
        _,ax = plt.subplots()
        ax.plot(self.t, getattr(self,vector))
        ax.set_xlabel()
        ax.set_ylabel()
        return ax

    def plot_magnitude_vs_time():
        pass

    def animateOrbit(self,**kwargs):
        
        xx = self.x
        yy = self.y

        history_len = len(self.t)
        fig, ax = plt.subplots() ## Initialise empty plots
        
        xdata, ydata = [], []
        ln, = plt.plot([], [], 'o', **kwargs)
        trace, = ax.plot([], [], '.-', lw=1, ms=2, **kwargs)
        history_x, history_y = deque(maxlen=history_len), deque(maxlen=history_len)

        def init():
            ax.set_xlim(min(xx), max(xx))
            ax.set_ylim(min(yy), max(yy))
            return ln,

        def motion(i): ## Fun Fact: this is a closure! It can access the variables in the outer function
            
            if i == 0:
                history_x.clear()
                history_y.clear()

            history_x.appendleft(xx[i])
            history_y.appendleft(yy[i])

            #xdata.append(xx[i])
            #ydata.append(yy[i])

            ln.set_data(xx[i], yy[i])
            trace.set_data(history_x, history_y)
            return ln, trace

        #a,_ = motion(1)
        #b,_ = motion(2)

        ani = FuncAnimation(fig, motion, len(self.y), interval = 20, blit=True, init_func=init)
        return ani

r0 = np.array([-5.11196475628565e-12, -12460.3473585584, -24882.7386997828]) * 1e3
v0 = np.array([4.38104966252361, -3.60350140900501e-16, -7.19602603237097e-16]) * 1e3

X0 = np.concatenate((r0,v0),axis = 0)

system = twoBody(X0,0,3600,1e6)
#system.plot_position_vs_time(vector='xyz')
system.animateOrbit()
plt.show()