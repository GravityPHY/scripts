import sys
import os
import numpy as np
import pandas as pd

f = 'bands_data.dat'

a=np.loadtxt(f)
pos=[]
neg=[]
k_p=[]
k_n=[]
for i in range(len(a)):
	if a[i,1]>0:
		pos.append(a[i,1])
		k_p.append(a[i,0])
	else:
		neg.append(a[i,1])
		k_n.append(a[i,0])
bandgap=np.min(pos)-np.max(neg)
print(k_p[np.argmin(pos)])
print(k_n[np.argmax(neg)])
print(bandgap)
