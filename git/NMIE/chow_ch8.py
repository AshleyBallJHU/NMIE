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
I = inflow_data['inflow']
add_I = np.zeros(len(I))
for i in range(len(I)):
    add_I[0] = 0.
    add_I[i+1] = I[i] + I[i+1]
    #delta_I = np.append(delta_I[i], delta_I)

#%%

S_Q_func_1 = np.zeros(len(I))
for i in range(len(I)):
    if i == 0:
        S_Q_func_1[i] = 0.
    else:
        S_Q_func_1[i] = I[i+1]+I[i]+2*S[i]/delta_t - Q[i]

