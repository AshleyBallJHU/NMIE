# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 11:08:21 2015

@author: user
"""

import numpy as np

#%%

# CH7 - INITAL VALUE PROBLEMS

## Taylor
'''define a function taylore which implements the taylor series method of 
integration of order four. the user must supply the function deriv that return
the 4xn array'''


'''
x, y = initial conditions
xStop = terminal value of x
h = increment of x used in integration
deriv = user-supplied function returning array
    [y'[0]    y'[1]    y'[2] ...  y'[n-1]
     y''[0]   y''[1]   y''[2] ... y''[n-1]
     ...
     ...]
'''

from numpy import array
def taylor(deriv, x, y, xStop, h):
    X = []
    Y = []
    X.append(x)
    Y.append(y)
    while x < xStop:
        h = min(h, xStop - x)
        D = deriv(x,y)
        H = 1.0
        for j in range(4):
            H = H*h/(j+1)
            y = y + D[j]*H
        x = x + h
        X.append(x)
        Y.append(y)
    return array(X), array(Y)
    

#%%

##Print Solutino
'''define a function to pring X and Y obtained from taylor integration. The
amount of data is controlled by the freq parameter.'''

'''
freq = prints every ___ timestep
'''

def printSoln(X,Y,freq):
    
    def printHead(n):
        print "\n      x ",
        for i in range(n):
            print "    y[",i,"]",
        print
    
    def printLine(x,y,n):
        print "%13.4e"%x,
        for i in range(n):
            print "%13.4e"%y[i],
        print
    
    m = len(Y)
    try: n = len(Y[0])
    except TypeError: n =1
    if freq == 0: freq = m
    printHead(n)
    for i in range(0,m,freq):
        printLine(X[i], Y[i], n)
    if i != m-1: printLine(X[m-1],Y[m-1],n)
        
        
#%%
        
###example 7_2
        
def deriv(x,y):
    D = np.zeros((4,2))
    D[0] = [y[1], -0.1*y[1] - x]
    D[1] = [D[0,1], 0.01*y[1] + 0.1*x - 1.]
    D[2] = [D[1,1], -0.001*y[1] - 0.01*x + 0.1]
    D[3] = [D[2,1], 0.0001*y[1] + 0.001*x - 0.01]
    return D

x = 0.0
xStop = 2.0
y = np.array([0.0, 1.0])
h = 0.25
freq = 1
X,Y = taylor(deriv, x, y, xStop, h)
printSoln(X,Y,freq)
raw_input("\nPress return to exit")


        
#%%
        
##Fourth order Runge-Kutta Method
''' defines a function integrate which uses the Runge-Kutta method of order 
four. A function F(x,y) must be provided that defines the eq y'=F(x,y)'''

'''
x, y = initial conditions
xStop = terminal value of x
h = increment of x used in integration
F = user-supplied function that returns the array
    F(x,y) = {y'[0], y'[1],...,y'[n-1]}
'''

def integrate(F,x,y,xStop,h):
    
    def run_kut4(F,x,y,h):
        K0 = h*F(x,y)
        K1 = h*F(x+h/2., y + K0/2.)
        K2 = h*F(x+h/2., y + K1/2.)
        K3 = h*F(x+h, y + K2)
        return (K0 + 2.*K1 + 2.*K2 + K3)/6.
    
    X = []
    Y = []
    X.append(x)
    Y.append(y)
    while x < xStop:
        h = min(h, xStop - x)
        y = y + run_kut4(F,x,y,h)
        x = x+h
        X.append(x)
        Y.append(y)
    return array(X), array(Y)
    
#%%
    
###example 7_4

def F(x,y):
    F = np.zeros((2))
    F[0] = y[1]
    F[1] = -0.1*y[1] - x
    return F

x = 0.
xStop = 2.
y = np.array([0.,1.])
h = 0.25
freq = 1
X,Y = integrate(F,x,y,xStop,h)
printSoln(X,Y, freq)
raw_input("Press return to exit")

#%%

###example 7_5

def F(x,y):
    F = np.zeros((1))
    F[0] = 3.0*y[0]-4.0*np.exp(-x)
    return F

x = 0.
xStop = 10.
y = np.array([1.])
h = 0.1
freq = 20
X,Y = integrate(F,x,y,xStop,h)
printSoln(X,Y, freq)
raw_input("Press return to exit")


#%%

###example 7_6

def F(x,y):
    F = np.zeros((4))
    F[0] = y[1]
    F[1] = y[0]*(y[3]**2) - 3.9860e14/(y[0]**2)
    F[2] = y[3]
    F[3] = -2.*y[1]*y[3]/y[0]
    return F

x = 0.
xStop = 1200.
y = np.array([7.15014e6, 0., 0., 0.937045e-3])
h = 50.
freq = 2

X,Y = integrate(F,x,y,xStop,h)
printSoln(X,Y,freq)
raw_input("Presse return to exit")


#%%

##RUNGE-KUTTA 5th ORDER

'''defines a function for the adaptive R-K integration method'''

'''
x,y = initial conditions
xStop = terminal value of x
h = increment of integration
tol = per step error tolerance
F = user-suppliued function
'''

from numpy import sqrt, zeros

def integrate(F,x,y,xStop,h,tol=1.e-6):
    
    def run_kut5(F,x,y,h):
        C = array([37./378, 0., 250./621, 125./594, 0., 512./1771])
        D = array([2825./27648, 0., 18575./48384, 13525./55296, 
                   277./14336, 1./4])
        n = len(y)
        K = zeros((6,n))
        K[0] = h*F(x,y)
        K[1] = h*F(x + 1./5*h, y + 1./5*K[0])
        K[2] = h*F(x + 3./10*h, y + 3./40*K[0] + 9./40*K[1])
        K[3] = h*F(x + 3./5*h, y + 3./10*K[0] - 9./10*K[1] + 6./5*K[2])
        K[4] = h*F(x + h, y - 11./54*K[0] + 5./2*K[1] - 70./27*K[2] + 
               35./27*K[3])
        K[5] = h*F(x + 7./8*h, y + 1631./55296*K[0] + 175./512*K[1] + 
               575./13824*K[2] + 44275./110592*K[3] + 253./4096*K[4])
        E = zeros((n))
        dy = zeros((n))
        for i in range(6):
            dy = dy + C[i]*K[i]
            E = E + (C[i] - D[i])*K[i]
        e = sqrt(sum(E**2)/n)
        return dy, e
        
    X = []
    Y = []
    X.append(x)
    Y.append(y)
    stopper = 0
    
    for i in range(10000):
        dy,e = run_kut5(F,x,y,h)
        if e <= tol:
            y = y + dy
            x = x + h
            X.append(x)
            Y.append(y)
            if stopper == 1: break
        if e != 0.0:
            hNext = 0.9*h*(tol/e)**2
        else: hNext = h
        if (h>0.0) == ((x + hNext) >= xStop):
            hNext = xStop - x
            stopper = 1
        h = hNext
    return array(X), array(Y)
    

#%%

###example 7_8

from numpy import exp

def F(x,y):
    F = zeros((2,))
    F[0] = y[1]
    F[1] = -9.80665 + 65.351e-3 * y[1]**2 * exp(-10.53e-5*y[0])
    return F
    
x = 0.
xStop = 10.
y = array([9000., 0.])
h = 0.5
freq = 10
X,Y = integrate(F,x,y,xStop,h,tol=1.0e-2)
printSoln(X,Y, freq)
raw_input("\nPress return to exit")


#%%

##MIDPOINT METHOD
'''defines a function that combines the midpoint method and RIchardson 
extrapolation. Starts with two integration steps, and in each subsequent step
integration steps are doubled, followed by RE until successive solutions 
differ by less than a designated tolerance.'''

'''
x,y = initial conditions
xStop = terminal value of x
yStop = y(xStop)
F = user-supplied function
'''

def integrate(F,x,y,xStop, tol):
    
    def midpoint(F,x,y,xStop,nSteps):
        h = (xStop - x)/nSteps
        y0 = y
        y1 = y0 + h*F(x,y0)
        for i in range(nSteps-1):
            x = x + h
            y2 = y0 + 2.*h*F(x,y1)
            y0 = y1
            y1 = y2
        return 0.5*(y1+y0+h*F(x,y2))
    
    def richardson(r,k):
        for j in range(k-1,0,-1):
            const = 4.**(k-j)
            r[j] = (const*r[j+1] - r[j])/(const - 1.)
        return
    
    kMax = 51
    n = len(y)
    r = zeros((kMax,n))
    nSteps = 2
    r[1] = midpoint(F,x,y,xStop,nSteps)
    r_old = r[1].copy()
    
    for k in range(2, kMax):
        nSteps = nSteps*2
        r[k] = midpoint(F,x,y,xStop,nSteps)
        richardson(r,k)
        
        e = sqrt(sum((r[1]-r_old)**2)/n)
        
        if e < tol: return r[1]
        r_old = r[1].copy()
    print "Midpoint method did not converge"
    

#%%

##BULIRCH-STOER METHOD
'''applies the midpoint method piecewise to get a more accurate soln.'''

'''
x,y = initial conds
xStop = terminal value of x
H = increment of x at which results are stored
F = user-supplied function
'''

def bulStoer(F,x,y,xStop,H,tol=1.0e-6):
    X = []
    Y = []
    X.append(x)
    Y.append(y)
    while x < xStop:
        H = min(H, xStop-x)
        y = integrate(F,x,y,x+H,tol)
        x = x + H
        X.append(x)
        Y.append(y)
    return array(X), array(Y)
    
#%%
    
###example 7.11
    
def F(x,y):
    F = zeros((2))
    F[0] = y[1]
    F[1] = (-y[1] -y[0]/0.45 + 9.)/2.
    return F

H = 0.5
xStop = 10.
x = 0.
y = array([0.,0.])
X,Y = bulStoer(F,x,y,xStop,H)
printSoln(X,Y,1)
raw_input("\nPress return to exit")

import matplotlib.pyplot as plt

plt.figure(711)
plt.plot(X, Y[:,1])
plt.xlabel('time (s)')
plt.ylabel('current (A)')
