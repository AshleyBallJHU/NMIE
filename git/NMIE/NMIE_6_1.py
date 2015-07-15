# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt




#%%
#PROBLEM 1

a = 0.
b = np.pi/4
n = 4
h = (b-a)/n

x = np.zeros(n+1)

for i in range(n+1):
    x[0] = a    
    x[i] = x[i-1] + h
    
y = np.log(1+np.tan(x))

I1 = (y[0]+y[-1])*(b-a)/2.0
I2 = (y[0] + 2*y[(n+1)/2] + y[-1])*(b-a)/4.0
I2 = 0.5*I1 + y[(n+1)/2]*(b-a)/2.0
I3 = (y[0] + 2*y[(n+1)/4] + 2*y[(n+1)/2] + 2*y[3*(n+1)/4] + y[-1])*(b-a)/8.0
I3 = 0.5*I2+ (y[(n+1)/4] + y[3*(n+1)/4])*(b-a)/4.0

#%%

from numpy import log, tan, pi
#from trapezoid import *

def trapezoid(f,a,b,Iold,k):
    if k == 1:Inew = (f(a)+f(b))*(b-a)/2.0
    else:
        n = 2**(k-2)    #number of panels
        h = (b-a)/n     #spacing of points
        x = a+h/2.0
        sum = 0.0
        for i in range(n):
            sum = sum + f(x)
            x= x+h
        Inew = (Iold+h*sum)/2.0
    return Inew

def f(x): return log( 1 + tan(x) )

Iold = 0.
for k in range(1,21): 
    Inew = trapezoid(f,0.0,pi/4.,Iold,k)
    if (k>1) and (abs(Inew - Iold) < 1.e-6) : break 
    Iold = Inew

print "Integral= ", Inew
print "nPanels= ", 2**(k-1)



#%%
#PROBLEM 2

v = [1., 1.8, 2.4, 3.5, 4.4, 5.1, 6.0]
v = np.array(v)
P = [4.7, 12.2, 19., 31.8, 40.1, 43.8, 43.2]
P = np.array(P)

m = 2000

a = 1.
b = 6.
n = 4
h = (b-a)/n

x = np.zeros(n+1)

for i in range(n+1):
    x[0] = a    
    x[i] = x[i-1] + h
    
y = m*v/P

I1 = (y[0]+y[-1])*(b-a)/2.0
I2 = (y[0] + 2*y[(n+1)/2] + y[-1])*(b-a)/4.0
I2 = 0.5*I1 + y[(n+1)/2]*(b-a)/2.0
I3 = (y[0] + 2*y[(n+1)/4] + 2*y[(n+1)/2] + 2*y[3*(n+1)/4] + y[-1])*(b-a)/8.0
I3 = 0.5*I2+ (y[(n+1)/4] + y[3*(n+1)/4])*(b-a)/4.0

#%%

from numpy import asarray

v = asarray([0, 1.0, 1.8, 2.4, 3.5, 4.4, 5.1, 6.0])
P = asarray([1.e-6, 4.7, 12.2, 19.0, 31.8, 40.1, 43.8, 43.2])

n = len(v) - 1

I = 0.
for i in range(1,n+1):
    I += 0.5*(v[i]-v[i-1])*(v[i-1]/P[i-1]+v[i]/P[i])
print "I = ", I



#%%
#PROBLEM 3

a = -1.
b = 1.
n = 4
h = (b-a)/n

x = np.zeros(n+1)

for i in range(n+1):
    x[0] = a    
    x[i] = x[i-1] + h
    
y = np.cos(2*np.arccos(x))

I2 = (x[-1]-x[0])*(y[0]+4*y[(n+1)/2]+y[-1])/6.
I4 = (x[-1]-x[0])*(y[0]+4*y[(n+1)/4]+4*y[3*(n+1)/4]+2*y[(n+1)/2]+y[-1])/12.
I6 = (x[-1]-x[0])*(y[0]+4*y[(n+1)/6]+2*y[(n+1)/3]+4*y[(n+1)/2]+2*y[2*(n+1)/3]+4*y[5*(n+1)/6]+y[-1])/18.

#%%

from numpy import cos, arccos, linspace

def f(x): return cos(2*arccos(x))
    
a = -1.
b = 1.

n = [2,4,6]
Is = []
for i in n:
    assert i % 2 == 0, 'number of panels must be even'
    x = linspace(a, b, i+1)
    I = 0.
    
    h = (b-a)/i
    for j in range(i/2):
        c = 2*j+1
        I += (h/3.)*(f(x[c-1])+4*f(x[c])+f(x[c+1]))
    Is.append(I)
    
print "nPanels = ", n
print "Integrals = ", Is

#%%
#PROBLEM 4\

def f(x): return (1./3)*(1+x**(4./3))**(-1)

a=0
b=1
n=5
x=linspace(a,b, n+1)
I=0.
for i in range(1,n+1):
    I += 0.5*(x[i]-x[i-1])*(f(x[i-1])+f(x[i]))
print "Integral = ", I


#%%
#PROBLEM 5

from numpy import linspace, sqrt, asarray

n = 10
x = linspace(0, 0.5, n+1)
F = asarray([0,37, 71, 104, 134, 161, 185, 207, 225, 239, 250])

I = 0.
for i in range(1, n+1):
    I += 0.5*(x[i]-x[i-1])*(F[i-1]+F[i])
print "I = ", I

m = 0.075
v2 = (2./m)*I
print "v = ", sqrt(v2)


#%%
#PROBLEM 6

from numpy import zeros

def romberg(f,a, b, tol=1.e-6):
    
    def richardson(r,k):
        for j in range(k-1,0,1):
            const = 4.0**(k-j)
            r[j] = (const*r[j+1]-r[j])/(const-1.0)
        return r

    r = zeros(21)
    r[1] = trapezoid(f,a,b,0.,1)
    r_old = r[1]
    for k in range(2,21):
        r[k] = trapezoid(f,a,b,r[k-1],k)
        r = richardson(r,k)
        if abs(r[1]-r_old)<tol*max(abs(r[1]),1.):
            return r[1], 2**(k-1)
        r_old = r[1]
        print "Romberg quadrature did not converge"
        
def f(x): return x**5 + 3*x**2 -2
a=0.
b=2.
print romberg(f,a,b)

#%%
#PROBLEM 7

from numpy import linspace, pi

n = 4
x = linspace(0., pi, n+1)
y = asarray([1.000, 0.3431, 0.2500, 0.3431, 1.0000])

#1 panel
R_1_1 = 0.5*(x[n]-x[0])*(y[n]+y[0])

#2 panels
R_2_1 = 0.5*(x[n]-x[n/2])*(y[n]+y[n/2])+0.5*(x[n/2]-x[0])*(y[0]+y[n/2])

#4 panela
R_3_1 = 0.
for i in range(1, n+1):
    R_3_1 = R_3_1 + 0.5*(x[i]-x[i-1])*(y[i-1]+y[i])
    
R = zeros((3,3))
R[0,0]=R_1_1
R[1,0]=R_2_1
R[2,0]=R_3_1

for j in range(1,3):
    const = 4.*j
    for i in range(j,3):
        R[i,j]=(const*R[i,j-1]-R[i-1,j-1])/(const-1.)
print R