from matplotlib import pyplot as plt
import numpy as np
from numpy import random as rn
import copy
import cmath

import pickle

def ztrans(i, N):
    z = complex(1, 1)/np.sqrt(2) * complex(0, 1) * np.exp(2*np.pi/N * i * complex(0,1))/2 + complex(0.5, 0.5)
    return z.real, z.imag



#init conditions

N = 30 # how many iterations
Nsides = 3
Nangle = []
for i in range(Nsides):
    Nangle.append(ztrans(i, Nsides))
x = 0.5 # init value
y = 0.5 # init value
Frac = 0.51




NSide = 100
curr = np.zeros((NSide, NSide)).astype("int")
prev = np.zeros((NSide, NSide)).astype("int")
cumu = np.zeros((NSide, NSide)).astype("int")


prev[50,50] = 1

#init conditions

N = 300 # how many iterations
Nsides = 4
Nangle = []
for i in range(Nsides):
    Nangle.append(ztrans(i, Nsides))
x = 0.5 # init value
y = 0.5 # init value
Frac = 0.55






lenNangle = len(Nangle)
randarr = rn.choice(lenNangle, N)


scatterdotsx = []
scatterdotsy = []

for iir,ir in enumerate(randarr):
    # find distance between random point
    for i in range(NSide):
        for j in range(NSide):
            if prev[i,j] == 0:
                continue
            x = i / float(NSide) # prev[i,j]
            y = j / float(NSide) # prev[i,j]

            newx, newy = x + (Nangle[ir][0] - x) * Frac, y + (Nangle[ir][1] - y) * Frac
            
            newx = np.min([1, np.max([0, newx])])
            newy = np.min([1, np.max([0, newy])])
             
            newi = int(np.floor(newx * NSide))
            newj = int(np.floor(newy * NSide))
            
            curr[newi, newj] += prev[i,j]
    for i in range(NSide):
        for j in range(NSide):
            
            
            cumu[i,j] += curr[i,j] 
            prev[i,j] = copy.copy(curr[i,j])
            #x = newpt[0]
            #y = newpt[1]
            #if iir < 20:
            #    print(newpt)
            #scatterdotsx.append(copy.copy(x))
            #scatterdotsy.append(copy.copy(y))
    print(iir)


    outtuple = (iir, prev, cumu)
    f = open('c:\\quant\\blah_' + str(iir) + '.pkl', 'wb')
    pickle.dump(outtuple, f)
    
    
import pdb; pdb.set_trace()
plt.contour(cumu)
plt.show()
    
"""
xangle = []
yangle = []
for ipt in Nangle:
    xangle.append(ipt[0])
    yangle.append(ipt[1])

# x-axis and y-axis values of the N-agon
xangle.append(Nangle[0][0])
yangle.append(Nangle[0][1])



"""


#for i in range(NSide):
#  for j in range(NSide):
#    for ih in range(cumu[i,j]):
#        scatterdotsx.append(i/float(NSide))
#        scatterdotsy.append(j/float(NSide))

# Creating bins
x_min = 0 # np.min(scatterdotsx)
x_max = 1 # np.max(scatterdotsx)
  
y_min = 0 # np.min(scatterdotsy)
y_max = 1 # np.max(scatterdotsy)
  
x_bins = np.linspace(x_min, x_max, 100)
y_bins = np.linspace(y_min, y_max, 100)
  
#fig, ax = plt.subplots(figsize =(10, 7))
# Creating plot
plt.hist2d(scatterdotsx, scatterdotsy, bins =[x_bins, y_bins], cmap = plt.cm.nipy_spectral)
plt.title("Changing the color scale and adding color bar")
  
# Adding color bar
plt.colorbar()
  
#ax.set_xlabel('X-axis') 
#ax.set_ylabel('X-axis') 
  
# show plot
plt.tight_layout() 
#plot.show()

