# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 14:34:51 2015

@author: user
"""

#CHAPTER 7 - PROBELM SET 1

import numpy as np
from numpy import array
import matplotlib.pyplot as plt

#%%

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
   
##3
from numpy import sin, cos

def deriv(x,y):
    D = np.zeros((4,1))
    D[0] = sin(y[0])
    D[1] = cos(y[0])*sin(y[0])
    D[2] = -sin(y[0])**3 + cos(y[0])**2 * sin(y[0])
    D[3] = -5*sin(y[0])**3 * cos(y[0]) + cos(y[0])**3 * sin(y[0])
    return D

x0 = 0.
xStop = 0.5
freq = 0
y0 = np.array([1.])
h = 0.1
X,Y = taylor(deriv, x0, y0, xStop, h)
printSoln(X,Y, freq)

if False:
    plt.plot(X,Y[:,0], color='black', label='approx.sol.')
    plt.hold('on')
    plt.xlabel('x'); plt.ylabel('y(x)')
    plt.grid('on')
    plt.legend()
    plt.show()
    

#%%

##4

from numpy import zeros

def F(x,y):
    F = zeros((1))
    F[0] = y[0]**(1./3.)
    return F
    
x0 = 0.
xStop = 5.
freq = 0

print "Part(a): y(0) = 0"
y0 = array([0.])
h = 0.01
X,Y = integrate(F,x0,y0,xStop,h)
printSoln(X,Y,freq)

print "Part (b): y(0) = 10^(-16)"
y0 = array([0.])
h = 0.01
X,Y = integrate(F,x0,y0,xStop,h)
printSoln(X,Y,freq)

Y_truth = ((2./3.)*X)**(3./2.)
print h, np.linalg.norm(Y_truth - Y[:,0])


#%%

##7

def F(x,y):
    F = zeros((2,))
    F[0] = y[1]
    F[1] = -sin(y[0])
    return F
    
x0 = 0.
xStop = 10.
y0 = array([1., 0.])
h = 0.05
freq = 20

X,Y = integrate(F,x0,y0,xStop,h)
printSoln(X,Y,freq)

#%%

##8

g = 9.80665
c_D = 0.2028
m = 80

def F(x,y):
    F = zeros((2,))
    F[0] = y[1]
    F[1] = - (g - (c_D/m)*y[1]**2)
    return F
    
x0 = 0.
xStop = 14.
y0 = array([500., 0.])
h = 0.1
freq = 20

X,Y = integrate(F, x0, y0, xStop, h)
printSoln(X,Y,freq)

y_sort_inds = np.argsort(Y[:,0])
y_sorted = Y[y_sort_inds,0]
x_sorted = X[y_sort_inds]
print "Time when y = 0 is: ", np.interp(0, y_sorted, x_sorted)

#%%

##9

k = 75
m = 2.5

def P(t):
    if t<2:
        return 10.*t
    else:
        return 20.
    
def F(x,y):
    F = zeros((2,))
    F[0] = y[1]
    F[1] = P(x)/m - (k/m)*y[0]
    return F
    
x = 0.
xStop = 6.
y0 = array([0.,0.])
h = 0.01
freq = 20

X,Y = integrate(F,x0,y0,xStop,h)
printSoln(X,Y,freq)

if False: 
    plt.plot( X, Y[:,0] )
    plt.xlabel('time (t)'); plt.ylabel('y(t)');
    plt.axhline(y=0,color='green')
    plt.grid('on')
    plt.show() 

#%%

##10

a = 16.
g = 9.80665

def F(x,y):
    F = zeros((2,))
    F[0] = y[1]
    F[1] = g*(1.-a*(y[0]**3))
    return F
    
x0 = 0.
xStop = 2.
y0 = array([0.1, 0.])
h = 0.1
freq = 20

X,Y = integrate(F,x0,y0,xStop,h)
printSoln(X,Y,freq)