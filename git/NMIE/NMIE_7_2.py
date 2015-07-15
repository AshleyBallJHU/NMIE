# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 12:33:34 2015

@author: user
"""

#CHAPTER 7 - PROBELM SET 2

#import numpy as np
from numpy import array, zeros, sin, cos, pi, exp, sqrt, random
import matplotlib.pyplot as plt

#%%

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


def integrate4(F,x,y,xStop,h):
    
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
    
    
def integrate5(F,x,y,xStop,h,tol=1.e-6):
    
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
    

def evalPoly(a,x):
    n = len(a) -1
    p = a[n]
    dp = 0. + 0.j
    ddp = 0. + 0.j
    for i in range(1,n+1):
        ddp = ddp*x + 2.*dp
        dp = dp*x + p
        p = p*x + a[n-1]
    return p,dp,ddp
    
def polyRoots(a, tol=1.e-12):
    
    def laguerre(a,tol):
        x = random()
        n = len(a)-1
        for i in range(30):
            p,dp,ddp = evalPoly(a,x)
            if abs(p) < tol: return x
            g = dp/p
            h = g*g - ddp/p
            f = sqrt((n-1)*(n*h - g*g))
            if abs(g+f) > abs(g-f): dx = n/(g+f)
            else: dx = n/(g-f)
            x = x - dx
            if abs(dx) < tol*max(abs(x),1.): return x
        print "Too many iterations in Laguerre"
    
    def deflPoly(a, root):
        n = len(a) -1
        b = [(0.+0.j)]*n
        b[n-1] = a[n]
        for i in range(n-1,-1,-1):
            b[i] = a[i+1]+root*b[i+1]
        return b
    
    n = len(a) -1
    roots = zeros((n))
    for i in range(n):
        x = laguerre(a,tol)
        if abs(x.imag) < tol: x = x.real
        roots[i] = x
        a = deflPoly(a,x)
    return roots
    
    
    
    
    
#%%
    
##3    
    
def F(x,y):
    F = zeros((1,))
    F[0] = x - 10*y[0]
    return F

x0 = 0.
xStop = 5.
freq = 0

y0 = array([10.])

y_truth = lambda x: 0.1*x - 0.01 + 10.01*exp(-10*x)
print "y_truth(%5.3f) = %5.3f" % (xStop,y_truth(xStop))

plt.hold('on')

for h in [0.1, 0.25, 0.5]:
    X,Y = integrate4(F,x0,y0,xStop,h)
    print "h = ", h
    printSoln(X,Y,freq)
    plt.plot(X,Y[:,0], label = 'h = %5.2f' % h)

plt.xlabel('x')
plt.ylabel('y(x)')
plt.grid('on')
plt.legend(loc='best')

#%%

##4

def F(x,y):
    F = zeros((1,))
    F[0] = x - 10*y[0]
    return F 

x0 = 0.
xStop = 10.
freq = 0

y0 = array([10.])

y_truth = lambda x: 0.1*x - 0.01 + 10.01*exp(-10*x)
print "y_truth(%5.3f) = %5.3f" % (xStop,y_truth(xStop))

plt.hold('on')

for h in [0.1, 0.25, 0.5]:
    X,Y = integrate5(F,x0,y0,xStop,h)
    print "h = ", h
    printSoln(X,Y,freq)
    plt.plot(X,Y[:,0], label = 'h = %5.2f' % h)

plt.xlabel('x')
plt.ylabel('y(x)')
plt.grid('on')
plt.legend(loc='best')


#%%

##5

m = 2.
c = 460.
k = 450.

#print "What are the decay rates for this diff eq?"
#print polyRoots[k/m, c/m, 1.]

def F(x,y):
    F = zeros((2,))
    F[0] = y[1]
    F[1] = -(k/m) * y[0] - (c/m) * y[1]
    return F

x0 = 0.
xStop = 0.2
freq = 0

y0 = array([0.01,0.])

h = 0.001
if True:
    X,Y = integrate4(F,x0,y0,xStop,h)
else:
    X,Y = integrate5(F,x0,y0,xStop,h)

printSoln(X,Y,freq)

if True:
    plt.plot(X,Y[:,1])
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid('on')


#%%

##7

def F(x,y):
    F = zeros((2,))
    F[0] = y[1]
    F[1] = 16.81*y[0]
    return F
    
x = 0.
xStop = 2.
freq = 0

y01 = array([1., -4.1])
y02 = array([1., -4.11])

h = 0.1
X1,Y1 = integrate5(F,x0,y01,xStop,h)
printSoln(X1,Y1,freq)
X2,Y2 = integrate5(F,x0,y02,xStop,h)
printSoln(X2,Y2,freq)

