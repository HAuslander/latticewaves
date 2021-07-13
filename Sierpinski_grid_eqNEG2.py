from matplotlib import pyplot as plt
import numpy as np
from numpy import random as rn
import copy
import cmath

import pickle

rn.seed(50)

def ztrans(i, N, bOpp = False):
    
    z = complex(1, 1)/np.sqrt(2) * complex(0, 1) * np.exp(2*np.pi/N * i * complex(0,1))/2 + complex(0.5, 0.5)
    if bOpp:
        if N % 2 != 0:
            print("ERROR: bOpp can only be True if nuber of sides of polygon is even")
        z = z * complex()
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


#init conditions

N = 30000 # how many iterations
Nsides = 4
Nangle = []
for i in range(Nsides):
    Nangle.append(ztrans(i, Nsides))
x = 0.5 # init value
y = 0.5 # init value
Frac = 0.52



curr = np.zeros((NSide, NSide, Nsides)).astype("int")
prev = np.zeros((NSide, NSide, Nsides)).astype("int")
prevrand = np.zeros((NSide, NSide, Nsides)).astype("int")
cumu = np.zeros((NSide, NSide)).astype("int") # accumulates along directions only
timecumu = np.zeros((NSide, NSide)).astype("int")

prev[50,50,0] = 1
prev[50,50,0] = 1


lenNangle = len(Nangle)
randarr = rn.choice(lenNangle, N)
randarr2 = rn.choice(lenNangle, N)


scatterdotsx = []
scatterdotsy = []

minnonzero = 10000
maxnonzero = 0
iir_glob = 0

def printnonzero(x, bPrint = True):
    global minnonzero
    global maxnonzero
    nnonzero = 0
    nonzeroval = []
    
    allval = 0
    for i in range(NSide):
        for j in range(NSide):
            totval = 0
            for kdir in range(Nsides):
                if x[i,j,kdir] == 0:
                    continue
                #print("ijkdir ", i, j, kdir, x[i,j,kdir])
                nnonzero += 1
                totval += x[i,j,kdir]
                allval += x[i,j,kdir]
            if totval != 0:
                nonzeroval.append(totval)
    if bPrint:
        plt.hist(nonzeroval, bins=np.arange(-120,151,2))
        plt.show()
    print ("Nnonzero ", nnonzero, " allval ", allval, " min ", minnonzero, " max ", maxnonzero)
    
    
    if iir_glob > 200:
        #import pdb; pdb.set_trace()
        if nnonzero < minnonzero:
            minnonzero = copy.copy(nnonzero)
            print("maxnonzero ", nnonzero, " min ", minnonzero, " max ", maxnonzero)
        if nnonzero > maxnonzero:
            maxnonzero = copy.copy(nnonzero)
            print("maxnonzero ", nnonzero, " min ", minnonzero, " max ", maxnonzero)
            
    


for iir,ir in enumerate(randarr):
    iir_glob = iir
    # find distance between random point
    for i in range(NSide):
        for j in range(NSide):
            for kdir in range(Nsides):
                if prev[i,j,kdir] == 0:
                    continue
                x = i / float(NSide) # prev[i,j]
                y = j / float(NSide) # prev[i,j]
    
                newx, newy = x + (Nangle[ir][0] - x) * Frac, y + (Nangle[ir][1] - y) * Frac
                
                
                newx = np.min([1, np.max([0, newx])])
                newy = np.min([1, np.max([0, newy])])
                
                newi = int(np.floor(newx * NSide))
                newj = int(np.floor(newy * NSide))
                
                curr[newi, newj, ir] += prev[i, j, kdir]
                #if prev[i, j, kdir] != 0:
                #    import pdb; pdb.set_trace()
                
                
                
                # 2nd random choice (because the kernel is not 1/n it is 2/n)
                ir2 = randarr2[iir]
                newx2, newy2 = x + (Nangle[ir2][0] - x) * Frac, y + (Nangle[ir2][1] - y) * Frac
                
                
                newx2 = np.min([1, np.max([0, newx2])])
                newy2 = np.min([1, np.max([0, newy2])])
                
                newi2 = int(np.floor(newx2 * NSide))
                newj2 = int(np.floor(newy2 * NSide))
                
                curr[newi2, newj2, ir2] += prev[i, j, kdir]
                #if (newi2, newj2) == (31, 68):
                #    import pdb; pdb.set_trace()
                #if prev[i, j, kdir] != 0:
                #    import pdb; pdb.set_trace()
                
                
                # 
                koppdir = (kdir + Nsides//2) % Nsides
                negnewx, negnewy = x + (Nangle[kdir][0] - x) * Frac, y + (Nangle[kdir][1] - y) * Frac
               
                
                negnewx = np.min([1, np.max([0, negnewx])])
                negnewy = np.min([1, np.max([0, negnewy])])
                
                negnewi = int(np.floor(negnewx * NSide))
                negnewj = int(np.floor(negnewy * NSide))
                
                #if (negnewi, negnewj) == (31, 68):
                #    import pdb; pdb.set_trace()
                
                
                curr[negnewi, negnewj, koppdir] += -prev[i, j, kdir]
                #if prev[i, j, kdir] != 0:
                #    print ("yo")
                #    import pdb; pdb.set_trace()

    #import pdb; pdb.set_trace()
                
    #if True: #iir > 25 and iir % 5 == 0:
    #    printnonzero(prev)
        
            
    cumu = np.zeros((NSide, NSide)).astype("int")        
    for i in range(NSide):
        for j in range(NSide):
            for kdir in range(Nsides):
            
            
                cumu[i,j] += curr[i, j, kdir] 
                timecumu[i,j] += curr[i, j, kdir] 
                prev[i,j, kdir] = copy.copy(curr[i,j, kdir])
                curr[i, j, kdir] = 0
                #x = newpt[0]
                #y = newpt[1]
                #if iir < 20:
                #    print(newpt)
                #scatterdotsx.append(copy.copy(x))
                #scatterdotsy.append(copy.copy(y))
    print(iir)


    #outtuple = (iir, prev, cumu)
    #f = open('c:\\quant\\blah_' + str(iir) + '.pkl', 'wb')
    #pickle.dump(outtuple, f)
    
    if iir > 25 and iir % 100 == 0:
        printnonzero(prev)
        #printnonzero(curr)

        plt.contour(cumu, levels = np.arange(-120, 121, 20))
        plt.show()
        
        
        
        plt.contour(timecumu, levels = np.arange(-120, 121, 20))
        plt.show()

        
        #plt.plot( np.sum(cumu, 0) )
        #plt.show()
        #plt.plot( np.sum(cumu, 1) )

        #print([np.sum(cumu[:,x]) for x in range(NSide)])
        #print([np.sum(cumu[x,:]) for x in range(NSide)])
        
        #print ("x", np.sum([np.sum(cumu[:,x]) for x in range(NSide)]))
        #print ("y", np.sum([np.sum(cumu[x,:]) for x in range(NSide)]))
        
        #import pdb; pdb.set_trace()
        plt.show()
        #import pdb; pdb.set_trace()
        
    else:    
        printnonzero(prev, False)
  
        
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

