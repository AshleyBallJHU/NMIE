# -*- coding: utf-8 -*-
"""
Created on Wed Jul 08 11:24:21 2015

@author: user
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

storage_data = pd.read_csv('chow_ch8.csv', sep=',')#, header = 0)

H = storage_data['H']
Q = storage_data['Q']
S = H*43560.

delta_t = 600.

S_Q_func = 2.*S/delta_t + Q


plt.figure(1)
plt.clf()
plt.plot(S_Q_func, Q, 'o-')
plt.xlabel('Storage Discharge Function [cfs]')
plt.ylabel('Outflow [cfs]')
plt.show()

#%%

inflow_data = pd.read_csv('chow_ch8_2.csv', sep=',')
t = inflow_data['t']
I = inflow_data['inflow']
add_I = np.zeros(len(I))
for i in range(len(I+1)):
    add_I[0] = 0.
    add_I[i+1] = I[i] + I[i+1]
    #delta_I = np.append(delta_I[i], delta_I)
print add_I

#%%

S_Q_func_1 = np.zeros(len(I))
Q_new = np.zeros(len(I))  

for i in range(len(I)):
    S_Q_func_1[0] = 0.
    S_Q_func_1[i] = I[i+1]+I[i]+2*S[i]/delta_t - Q_new[i]
print S_Q_func_1

#%%

S_Q_func_2 = np.ones(len(I))

for i in range(len(I)):
    S_Q_func_2[0] = 0.
    S_Q_func_2[i] = add_I[i] + S_Q_func_1[i-1]
print S_Q_func_2
        
#%%        

S_new = S_Q_func_2
iprev = np.zeros(len(I))
inext = np.zeros(len(I))

for i in range(len(I)):
    if S_new[i] < S_Q_func[0]: 
        iprev[i] = 0

    elif S_new[i] > S_Q_func[-1]:
        iprev[i] = len(S_new) - 1
    else:
        iprev[i] = np.nonzero(S_Q_func <= S_new[:-1])[0][-1]
        
        '''the issue is here. this function is comparing all of the values in
        S_Q_func with all of the values in S_new according to their indices. 
        Instead we want for a single value in S_new to be compared to all the 
        values in S_Q_func...'''

    inext[i] = iprev[i] + 1       

    Q_new[0] = 0.
    Q_new[i] = Q[iprev[i]] + ((Q[inext[i]]-Q[iprev[i]])/(S_Q_func[inext[i]]-
                   S_Q_func[iprev[i]]))*(S_new[inext[i]] - S_Q_func[iprev[i]])
print Q_new


#%%

plt.figure(9)
plt.clf()
plt.plot(t, I, 'bo-', label = 'inflow [cfs]')
plt.plot(t, Q_new, 'ro-', label = 'outflow [cfs]')
plt.legend()
