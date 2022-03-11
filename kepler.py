from asyncio.windows_events import NULL
from email import generator
import numpy as np
from scipy.optimize import newton
import matplotlib.pyplot as plt

#from utils import setattrs

def trueAnomaly(E,e):
    '''
    True anomaly is the actual measured angle *theta* in the orbital plane between a vector from the focus to periapsis (closest point)
    and a vector joining the focus and an objects actual position.
    
    :param: float E eccentric anomaly (radians) - see below
    :param: float e eccentricity

    :return: float theta the true anomaly
    '''
    theta = 2 * np.arctan( np.tan(.5 * E) / np.sqrt( (1 - e)/(1 + e) ) )    
    return theta if (theta >= 0.0 ) else theta + np.pi

def radialDistance(theta,a,e):
    '''
    Radial distance from the focus
    Derived using equations:  
        (1) p = a*(1-e^2)
        (2) r = p/(1-e)
    
    :param: float a semi-major axis (m)
    :param: float e eccentricity of the orbit
    
    :return: float r the radius of orbit measured from the body which is being orbited.
    '''
    r = ( a*(1 - e**2) )/(1 + e*np.cos(theta))
    return r

def meanAnomaly(t,t0,mu,a):
    '''
    This calculates the mean anomaly
    Mean anomaly is the angle between the periapsis point and the *imagined* position of an object for the same elapsed time 
    since periapsis for a **circular orbit** around the same body with the same orbital period.
    
    M = n(t-t0) where n = root(mu/a^3)

    :param: float t time (s)
    :param: float t0 time of periapsis passage (s)
    :param: float gravitational parameter
    :param: float a semi-major axis (m)

    :return: M the mean anomaly
    '''
    return np.sqrt( mu/a**3 )*(t-t0)

def degToRad(deg):
    return deg*(np.pi/180)

def Kepler(M,e,tol = 1e-3,N = 500,v = False):
    """
    Solves Kepler's 1st equation M - e*sin(E) = E, giving the Eccentric anomaly.
    Use Newton's method to solve f(E) = 0 = M - E - e*sin(E), i.e. find the root E.
    The Eccentric anomaly is the angle between the **centre** of an ellipse and a point 
    on a circle with radius a for a **circular orbit** around the same body with the same orbital period.

    :param: float M
    :param: float e
    :param: float tol
    :return: float E
    """
    E = newton(func = lambda E: M - e*np.sin(E) - E, x0 = M, tol = tol, full_output=True, maxiter= N, disp = False)

    if v:
        print('Converged: {}\nNumber of iterations for tol {:.2e} : {:03d}'.format(E[1].converged,tol,E[1].iterations))

    return E[0]


class Orbit():
    def __init__(self,time: np.array, e: float, a: float, mu: float) -> None:
        ### Constants - corresponding to orbital parameters
        self.e = e                      ## eccentricity
        self.a = a                      ## semi-major axis (m)
        self.mu = mu                    ## gravitational parameter of central body
        self.time = time                ## time array
    
    def __repr__(self) -> str:
        pass

    def __update_orbital_params(self,ti: float) -> tuple: ## private function
        """
        :param: float ti
        :return: tuple
        """
        Mk = meanAnomaly(ti, self.time[0], self.mu,self.a)
        Ek = Kepler(Mk,self.e)
        thetak = trueAnomaly(Ek,self.e)
        r1 = radialDistance(thetak,self.a, self.e)
        return  (Mk, Ek, thetak, r1)

    def __iter__(self):                 ## generator function (for efficiency purposes)
        for i in range(len(self.time)): ## iteratively update params
            yield self.__update_orbital_params(self.time[i])


class visOrbit(Orbit):

    def __init__(self,time: np.array,**kw) -> None:
        super().__init__(time,**kw)        ## access all the attributes from the inherited class
        vals = (params for params in self) ## unpack the tuple values from: ((a1,b1,c1,d1),(a2,b2,c2,d2),...) into ((a1,a2,a3,a4),(b1,b2,b3,b4),...)
        (self.mean_anomaly, self.eccentric_anomaly, self.true_anomaly, self.radial_distance) = tuple(zip(*vals))

    def __repr__(self) -> str:
        ##return super().__repr__()
        pass

    def plotAltitude(self,
                     textsize = 8,
                     fontweight = 'bold'
                     ):
        _,ax = plt.subplots()
        ax.plot(
            np.asarray(self.time)/3600,            ## time in hours
            np.asarray(self.radial_distance)/1e3   ## altitude to km
            )
        
        ax.set_title('Altitude vs time', fontweight = fontweight, loc = 'center', fontsize = textsize*2)
        ax.set_ylabel('Altitude (km)')
        ax.set_xlabel('Time (Hours)')
        ax.ticklabel_format(useOffset=False,style = 'plain')
        return ax
        
    def plotAnomaly(self,
                    which = "true",
                    textsize = 8,
                    fontweight = 'bold',
                    fill = True,
                    alpha = 0.2,
                    fillcolour = 'blue'
                    ): ##TODO option to choose which anomaly is plotted
        
        ## Choose anomaly to plot
        angle = self.true_anomaly if which == "true" else self.mean_anomaly if (which == "mean") else self.eccentric_anomaly

        _,ax = plt.subplots()
        ax.plot(
            np.asarray(self.time/3600),
            np.asarray(angle)
            )

        ax.set_title('{:s} Anomaly vs Time'.format(which.title()), fontweight = fontweight, loc = 'center', fontsize = textsize*2)
        ax.set_ylabel('{:s} Anomaly (radians)'.format(which.title()), fontsize = textsize)
        ax.set_xlabel('Time (Hours)', fontsize = textsize)
        ax.ticklabel_format(useOffset=False,style= 'plain')
        
        if fill is True:
            plt.fill_between(self.time/3600,0, angle, color=fillcolour, alpha=alpha)
        return ax
    
    def animateOrbit(self):
        pass
    

Q1 = False
###############################################################################################
if Q1 :
    # Question 1.)
    # a.) above ^^

    # b.)
    kwargs = {'e': 0.74,'tol': 1e-3,'v': True} ## why? just coz B)
    E = Kepler(M = degToRad(230),N = 1, **kwargs)
    print('E is: {:.3f} \n'.format(E))

    ## iterate over multiple tolerances for comparison
    Es = list( map( lambda eps: Kepler(M = degToRad(230),e = 0.74, tol = eps), [1e-3,1e-12]) ) ## why? just coz B)
    print(*Es)

    ## Get the radial distance and true anomaly
    theta = trueAnomaly(Es[1], 0.74)
    r = radialDistance(theta,a = 245e5, e = 0.74)
    print('230 degrees\ntrueAnomaly: {:3f}\nradial distance: {:.3f}\n'.format(theta,r))

    ## repeat for M = 180
    E180 = Kepler(M = degToRad(180),e = 0.74, tol = 1e-3)
    theta180 = trueAnomaly(E180, 0.74)
    r = radialDistance(theta180,a = 245e5, e = 0.74)
    print('180 degrees\ntrueAnomaly: {:3f}\nradial distance: {:.3f}\n'.format(theta,r))

###############################################################################################
#Question 2.)

const = {
    'a': 24000e3,
    'e': 0.735,
    'mu': 3.986e5
    }

## orbital periods
T = (lambda a,mu : 2*np.pi*np.sqrt(a**3/mu))(const['a'],const['mu'])
##time = np.arange(0, 2*T, 1e4)

timeSmall = np.arange(0,36e3,60)
myPlot = visOrbit(timeSmall,**const)

myPlot.plotAltitude()
plt.show()

myPlot.plotAnomaly()
plt.show()