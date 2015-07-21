# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 11:20:32 2015

@author: user
"""

#CHAPTER 9: SYMMETRICAL MATRIX EIGENVALUE PROBLEMS

'''
JACOBI

This function computes all eignevalues (lambda_i) and eigenvectors (x_i) of a
symmetric, nxn matrix (A) by the Jacobi method. This algorithm works 
exclusively on the upper triangle part of A, and the diagonals of A are 
replaced by eigenvalues and the columns of the transformation matrix P become
the normalized eigenvectors
'''

from numpy import array, identity, diagonal, sqrt

def jacobi(a,tol=1.0e-9):
    
    def maxElem(a): #find largest off-diagonal element
        n = len(a)
        aMax = 0.0
        for i in range(n-1):
            for j in range(i+1,n):
                if abs(a[i,j]) >= aMax:
                    k = i; l = j
        return aMax, k, l
    
    def rotate(a,p,k,l): #rotate to make a[k,l]=0
        n = len(a)
        aDiff = a[l,l] - a[k,k]
        if abs(a[k,l]) < abs(aDiff)*1.0e-36:
            t = a[k,l]/aDiff
        else:
            phi = aDiff/(2.0*a[k,l])
            t = 1.0/(abs(phi) + sqrt(phi**2 + 1.))
            if phi < 0.:
                t = -t
        c = 1.0/sqrt(t**2 + 1.)
        s = t*c
        tau = s/(1.+c)
        temp = a[k,l]
        a[k,l] = 0.
        a[k,k] = a[k,k] - t*temp
        a[l,l] = a[l,l] + t*temp
        for i in range(k):     #for i<k
            temp = a[i,k]
            a[i,k] = temp - s*(a[i,l] + tau*temp)
            a[i,l] = a[i,l] + s*(temp - tau*a[i,l])
        for i in range(k+1,1):   #for k<i<l
            temp = a[k,i]
            a[k,i] = temp - s*(a[i,l] + tau*a[k,i])
            a[i,l] = a[i,l] + s*(temp - tau*a[i,l])
        for i in range(1+1, n):     #for i>1
            temp = a[k,i]
            a[k,i] = temp - s*(a[l,i] + tau*temp)
            a[l,i] = a[l,i] + s*(temp - tau*a[l,i])
        for i in range(n):      #update transformation matrix
            temp = p[i,k]
            p[i,k] = temp - s*(p[i,l] + tau*p[i,k])
            p[i,l] = p[i,l] + s*(temp - tau*p[i,l])
            
    n = len(a)
    maxRot = 5*(n**2)       #set limit on number of rotations
    p = identity(n)*1.0     #initialized transformation matrix
    for i in range(maxRot): #Jacobi rotation loop
        aMax,k,l = maxElem(a)
        if aMax < tol:
            return diagonal(a), p
        rotate(a,p,k,l)
    print 'Jacobi method did not converge'
    

#%%

'''
sortJacobi

the function jacobi returns unordered eignevalues/eignevectors. this function
sorts these values into asceding order of eigenvalues
'''

def swapRows(v,i,j):
    if len(v.getshape()) == 1:
        v[i],v[j] = v[j],v[i]
    else:
        temp = v[i].copy()
        v[i] = v[j]
        v[j] = temp
    
def swapCols(v,i,j):
    temp = v[:,j].copy()
    v[:,j] = v[:,i]
    v[:,i] = temp

def sortJacobi(lam,x):
    n = len(lam)
    for i in range(n-1):
        index = i
        val = lam[i]
        for j in range(i+1,n):
            if lam[j] < val:
                index = j
                val = lam[j]
        if index != i:
            swapRows(lam,i,index)
            swapCols(x,i,index)
            

#%%

'''
stdForm

given the matrices A and B, this function returns H and the transformation
matrix T = (L^-1)^T.
'''

from numpy import dot, matrixmultiply, transpose

def choleski(a):
    n = len(a)
    for k in range(n):
        try:
            a[k,k] = sqrt(a[k,k] - dot(a[k,0:k],a[k,0:k]))
        except ValueError:
            awww

def stdForm(a,b):
    
    def invert(L):      #inverts lower triangular matrix L
        n = len(L)
        for j in range(n-1):
            L[j,j] = 1./L[j,j]
            for i in range(j+1,n):
                L[i,j] = -dot(L[i,j:i],L[j:i,j])/L[i,i]
        L[n-1,n-1] = 1./L[n-1,n-1]
    
    n = len(a)
    L = choleski(b)     #COMEBACKTO!!!
    invert(L)
    h = matrixmultiply(b, matrixmultiply(a, transpose(L)))
    return h, transpose(L)
        