import numpy as np
from scipy.optimize import newton
import matplotlib.pyplot as plt

def trueAnomaly(E,e): ## theta
    return 2 * np.arctan( ( 1/np.sqrt(( (1 - e)/(1 + e) )) ) * np.tan(.5 * E) )

def radialDistance(theta,a,e):
    '''
    Derived using equations:  
        (1) p = a*(1-e^2)
        (2) r = p/(1-e)
    
    :param: float a semi-major axis (m)
    :param: float e eccentricity of the orbit
    
    :return: r
    '''
    r = a*(1 - e**2)/(1 + e*np.cos(theta))
    return r

def meanAnomaly(t,t0,mu,a):
    return np.sqrt( (mu/a**3)*(t-t0) )

def degToRad(deg):
    return deg*(np.pi/180)

def Kepler(M,e,tol = 1e-3,N = 500,v = False):
    """
    Solves Kepler's 1st equation M - e*sin(E) = E
    Use Newton's method to solve f(E) = 0 = M - E - e*sin(E), i.e. find the root E.

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
    def __init__(self,time,e,a,M,mu,t0 = 0):
        Ei = Kepler(M,e) ## trueAnomaly

        ## Attributes - corresponding to orbital parameters
        ### Constants
        self.e = e                      ## eccentricity
        self.a = a                      ## semi-major axis (m)
        self.mu = mu                    ## gravitational parameter of central body

        ## Initial Orbital motion params
        self.M = M                      ## mean anomaly (radians)
        self.E = Ei                     ## eccentric Anomaly (radians)
        self.theta = trueAnomaly(Ei,e)  ## true Anomaly (radians)
        self.r = radialDistance(theta,a,e)    ## radial distance

        ## ARRAYS TO STORE INFO 
        self.time = time
        self.altitude = np.zeros(len(time))
        self.trueAnomaly = np.zeros(len(time))

    def __update_orbital_params(self,ti): ## private function
        """
        :param: M0 float the mean anomaly at time t
        """
        M1 = meanAnomaly(ti, self.time[0], self.mu,self.a)
        E1 = Kepler(M1,self.e)
        theta1 = trueAnomaly(E1,self.e)
        r1 = radialDistance(theta1,self.a, self.e)
        return  E1, r1

    def calcOrbitalMotion(self,time):
        
        ## 1st iteration
        E1, r1 = self.__update_orbital_params(self.time[0])

        for i in range(0,len(time)-1):
            ## subsequent update step
            E1, r1 = self.__update_orbital_params(time[i])
            self.trueAnomaly[i] = E1
            self.altitude[i] = r1
        return

    def plotAltitude(self):
        _,ax = plt.subplots()
        ax.plot(
            self.time/3600,          ## hours
            self.altitude/1e3        ## altitude to km
            )
        ax.set_ylabel('Altitude (km)')
        ax.set_xlabel('Time (Hours)')
        ax.ticklabel_format(useOffset=False)
        return ax
        

    def plotAnomaly(self):
        _,ax = plt.subplots()
        ax.plot(
            self.time/3600,
            self.trueAnomaly
            )
        ax.set_ylabel('Anomaly (radians)')
        ax.set_xlabel('Time (Hours)')
        ax.ticklabel_format(useOffset=False)
        plt.fill_between(np.arange(0,len(self.trueAnomaly)), self.trueAnomaly, color='blue', alpha=0.3)
        return ax
    


###############################################################################################
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
    'M': 0,
    'a': 24000e3,
    'e': 0.735,
    'mu': 3.986e5,
    't0': 0.0,
    }

T = (lambda a,mu : 2*np.pi*np.sqrt(a**3/mu))(const['a'],const['mu']) ## orbital periods

time = np.arange(0,2*T,15)

timeSmall = np.arange(0,3.6e5,60)

MEO = Orbit(timeSmall,**const)
MEO.calcOrbitalMotion(timeSmall)

MEO.plotAltitude()
plt.show()

MEO.plotAnomaly()
plt.show()