#!/cygdrive/c/Python27/python

print("""

DEPRECATED: this is not only much slower than latticegas_youngsdoubleslitexp_penalty but the
need for regularization is much more pronounced, i.e. even though the sum of the amplitudes are
conserved, it's just the sum of increasingly large negative and positive grid values.


""")


print("""

perform basic Brownian-Huygens (in an effort to see whether 
total number of particles stabilize, maybe with a 
Pauli Exclusion Principle, or else, demanding that
Brownian particles always take oppositely oriented
steps)

END RESULT: in 2dim (and presumably all higher ones, the absolute number of particles keeps increasing)

""")


import sys
import os
import optparse
import copy

import numpy as np

import os.path
import datetime, time

import glob

import pickle 
import pandas as pn    

import matplotlib.pyplot as plt      
from matplotlib.backends.backend_pdf import PdfPages

import matplotlib.pyplot as plt
import matplotlib.image as mpimg

#import statgrid2dclassical # don't need this cython routine after all


bUseLatticeGasForRand = False
bCplusplus = False

f = None #open('c:\\quant\\blah4.txt', 'r')

def print2(x):
  if not(f is None):
    f.write(x)
  print(x)

def randamp(shapetuple, amp):
  if amp == 0:
      return np.zeros(shapetuple)
  absamp = np.abs(amp * np.ones(shapetuple))
  sgnamp = np.sign(amp * np.ones(shapetuple))
  
  if np.max(absamp) > 1.0:
      import pdb; pdb.set_trace()
      print("ERROR: randamp only works for functions whose values are in [-1, 1]" + str(absamp))
  
  rnd = np.random.random_sample(shapetuple)
  #import pdb; pdb.set_trace()
  return sgnamp * (rnd < absamp).astype("int")
  #if np.random.rand() < absamp:
  #    return sgnamp 
  #else:
  #    return 0



def ReadParams():

  p = optparse.OptionParser()
    
  p.add_option("-f", "--file", default="",
                  action="store", type="string", dest='file',
                  help="Can be either a single file (in which case it is assumed to be a list of filenames) or a comma-delimited list of files")

  p.add_option("-o", "--output", default="",
                  action="store", type="string", dest='output',
                  help="output file")

  p.add_option("--startdt", default="",
                  action="store", dest='startdt',
                  help="yyyy-mm-dd")

  p.add_option("--enddt", default="",
                  action="store", dest='enddt',
                  help="yyyy-mm-dd")


  p.add_option("--dir", default="/Users/joe/quant/Series/",
                  action="store", dest='dir',
                  help="default is /Users/joe/quant/Series/ ")


  
  p.add_option("--minSR", "--minsr", "--cutoffSR", "--cutoffsr", default=0.0,
                  action="store", dest='minSR', type='float',
                  help="default is /Users/joe/quant/Series/ ")


  p.add_option("--targetSD", "--targetSTD", "--targetstdev",  default=1.0,
                  action="store", dest='targetSD', type='float',
                  help="default is /Users/joe/quant/Series/ ")


  
  p.add_option("--templateHO", default="/Users/joe/quant/BRtemplate_current_LIVE.txt",
                  action="store", dest='templatefileHO', type='string',
                  help="default is /Users/joe/quant/BRtemplate_current.txt ")

  
  p.add_option("--templateEODC", default="/Users/joe/quant/BRtemplateEODC_LIVE.txt",
                  action="store", dest='templatefileEODC', type='string',
                  help="default is /Users/joe/quant/BRtemplateEODC.txt ")


  
  p.add_option("--tag", default="LIVE",
                  action="store", dest='tag', type='string',
                  help="will be added to outputfiles instead of @@@1 ")


  p.add_option("--maxgap", default=10,
                  action="store", dest='maxgap', type='int',
                  help="will be put into maxgap slot ")


  p.add_option("--dumptag", default="",
                  action="store", dest='dumptag', type='string',
                  help="THE FIRST PART OF THIS MUST BE A PATHNAME, e.g. /home/joe/joe_feathers/feathers-client/input/DUMP, and  will be used to create --dumpfile and --dumpbracketfile args ")


  p.add_option("--loadtag", default="",
                  action="store", dest='loadtag', type='string',
                  help="THE FIRST PART OF THIS MUST BE A PATHNAME, e.g. /home/joe/joe_feathers/feathers-client/input/DUMP, and will be used to create --loadfile and --loadbracketfile args ")


  p.add_option("--OvNtextrastring", default="",
                  action="store", dest='OvNtextrastring', type='string',
                  help="gives you extra restrictions on the OvNt trades ")


  p.add_option("--EODCextrastring", default="",
                  action="store", dest='EODCextrastring', type='string',
                  help="gives you extra restrictions on the EODC trades ")


  p.add_option("--picklefile", default="",
                  action="store", dest='picklefile', type='string',
                  help="pickle file, to which will be added the random seed ")


  p.add_option("--minnonzero", default=20,
                  action="store", dest='minnonzero', type='int',
                  help="minimum nonzero days -- added this because DBEF only had 2 non zero trades in 2015 ")

  p.add_option("--seed", default=84848484,
                  action="store", dest='seed', type='int',
                  help="rand number seed ")

  p.add_option("--discrete", default=False,
                  action="store_true", dest='bDiscrete', 
                  help="are values integers or not ")

  p.add_option("--discreet", default=False,
                  action="store_true", dest='bDiscrete', 
                  help="are values integers or not-- by the way, the correct word is DISCRETE, idiot ")


  
                
  return p.parse_args()                

def zero_out_parity(grid, zeroORone):
  for i in range(grid.shape[1]):
    for j in range(grid.shape[2]):
      if (i + j) % 2 == zeroORone:
        grid[:,i,j,:] = 0

def zero_out_paritySUB(grid, zeroORone):
  for i in range(grid.shape[1]):
    for j in range(grid.shape[2]):
      if (i + j) % 2 == zeroORone:
        grid[i,j,:] = 0

def sumpar(grid, zeroORone):
  sumgrid = 0
  for i in range(grid.shape[1]):
    for j in range(grid.shape[2]):
      if (i + j) % 2 == zeroORone:
        sumgrid += np.sum(grid[i,j,:])
  return sumgrid
    
    
Comp000 = 0 # completely arbitrary
Comp090 = 1 # completely arbitrary
Comp180 = 2 # completely arbitrary
Comp270 = 3 # completely arbitrary       


NTEST = 1



# this is for graphing the path taken by a single particle -- it has no other purpose
class Tracker:
  def __init__(self, i_scenario, i_init, j_init, i_dir):
    self.i = i_init
    self.j = j_init
    self.dir = i_dir
    self.scenario = i_scenario
    
    self.history = []
    self.NMaxPath = 88889888
    
    self.bPropagateInSameDirectionIfPossible = True


  def addhistory2grid(self, grid, xval):
    #print "  historylength " + str(len(self.history))
    #import pdb; pdb.set_trace()
    
    for i_hist in range(len(self.history)-1, np.max([-1, len(self.history) - self.NMaxPath - 1]), -1):
      itup = self.history[i_hist]
      grid[(2*itup[0]):(2*itup[0] + 2), (2*itup[1]):(2*itup[1] + 2)] = xval 
    
  # first you update the graph (so that all the arcs are OUTGOING arcs), THEN you call this thing
  def postupdate(self, grid):
  
    old_i = copy.copy(self.i)
    old_j = copy.copy(self.j)
    old_dir = copy.copy(self.dir)
  
  	# use old dir to get new i,j coordinates, then randomly select new dir from the non-zero arcs there
    if self.dir == Comp000:
      self.i = (self.i + 1) % grid.shape[1]
      #self.j = self.j
      
    if self.dir == Comp090:
      #self.i = i 
      self.j = self.j + 1 % grid.shape[2]
      
    if self.dir == Comp180:
      self.i = self.i - 1
      if self.i < 0:
        self.i += grid.shape[1]
      #j = j
      
    if self.dir == Comp270:
      #i = i 
      self.j = self.j - 1
      if self.j < 0:
        self.j += grid.shape[1]
    
    if self.bPropagateInSameDirectionIfPossible:
      if grid[self.scenario, self.i, self.j, self.dir] != 0:
        self.history.append((self.i, self.j))
        return
        

















      
    nonzeropresencelist = []
    nnonzero = 0
    for idir in range(4):
      if grid[self.scenario, self.i, self.j, idir] != 0:
        #import pdb; pdb.set_trace()
        for jdir in range(np.abs(int(grid[self.scenario, self.i, self.j, idir]))):
        	nonzeropresencelist.append(idir)
        	nnonzero += 1
        	
    if nnonzero != 0:   
      #import pdb; pdb.set_trace() 	
      randint = int(np.random.random() * nnonzero)
      self.dir = nonzeropresencelist[randint]
      if len(self.history) > 0:
        lastitup = self.history[-1]
        if lastitup[0] != self.i or lastitup[1] != self.j:
          self.history.append((self.i, self.j))
      else:
        self.history.append((self.i, self.j))
    else:
      self.i = copy.copy(old_i)
      self.j = copy.copy(old_j)
      self.dir = copy.copy(old_dir)
    

    return
    
    
    






def flatten2dgridandcsv(grid, bAddBorder, t, FileName=None):
    # will only flatten the first two scenarios
    #import pdb; pdb.set_trace()
    intAddBorder = int(bAddBorder)
    flatgrid = np.zeros((grid.shape[1]*2 + intAddBorder, grid.shape[2]*2 + intAddBorder)).astype("int")
    for ix in range(grid.shape[1]):
        for iy in range(grid.shape[2]):
        
            #if (ix,iy)==(0,12):
            #    import pdb; pdb.set_trace()
            #if ix >= 30:
            #    import pdb; pdb.set_trace()
            #if (ix + iy + t) % 2 == 1:
            #    continue
            # bra vals (scenario 0)
            
            #if np.sum(np.abs(grid[0, ix, iy, :])) != 0 and ix == 2 and t == 0:
            #    import pdb; pdb.set_trace()
            #if(ix==0 and 0!=np.sum(np.abs(grid[0, ix, iy,:]))):
            #    import pdb; pdb.set_trace()
            flatgrid[intAddBorder + 2*ix + 1, intAddBorder + 2*iy + 1 + 0] = grid[0, ix, iy, 0]
            flatgrid[intAddBorder + 2*ix + 0, intAddBorder + 2*iy + 1 + 0] = grid[0, ix, iy, 1]
            flatgrid[intAddBorder + 2*ix + 0, intAddBorder + 2*iy + 0 + 0] = grid[0, ix, iy, 2]
            flatgrid[intAddBorder + 2*ix + 1, intAddBorder + 2*iy + 0 + 0] = grid[0, ix, iy, 3]
            
            #if flatgrid[2,25] != 0 and t == 0:
            #    import pdb; pdb.set_trace()

    
            if ix == grid.shape[1] - 1:
                newix = 0
            else:
                newix = ix
            if iy == grid.shape[2] - 1:
                newiy = 0
            else:
                newiy = iy
            
            #flatgrid[intAddBorder + 2*newix + 1, intAddBorder + 2*newiy + 1 + 2] = grid[1, ix, iy, 0]
            #flatgrid[intAddBorder + 2*newix + 0, intAddBorder + 2*newiy + 1 + 2] = grid[1, ix, iy, 1]
            #flatgrid[intAddBorder + 2*newix + 0, intAddBorder + 2*newiy + 0 + 2] = grid[1, ix, iy, 2]
            #flatgrid[intAddBorder + 2*newix + 1, intAddBorder + 2*newiy + 0 + 2] = grid[1, ix, iy, 3]
            #if flatgrid[2,25] != 0 and t == 0:
            #    import pdb; pdb.set_trace()
            
            
    if bAddBorder:
        flatgrid[0,0] = t
        for ix in range(grid.shape[1]):
            flatgrid[intAddBorder + 2*ix, 0] = ix
        for iy in range(grid.shape[2]):
            flatgrid[0, intAddBorder + 2*iy] = iy
            
    if not(FileName == None):
        if os.path.exists(FileName):
            with open(FileName, "ab") as f:
                f.write(b"\n")
                np.savetxt(f, flatgrid, delimiter=',')
                #import pdb; pdb.set_trace()
                
        else:
            
            #import pdb; pdb.set_trace()
            with open(FileName, "w") as f:
                #f.write(b"\n")
                #import pdb; pdb.set_trace()
                np.savetxt(f, flatgrid, delimiter=',')
    #import pdb; pdb.set_trace()
    
# used for testing only    
def arbitrary_lattice_input_output(inputamp):
    
    Narc = len(inputamp)
    fac = 2.0 / Narc
    sumamp = np.sum(inputamp)
    outputamp = []
    for i,iamp in enumerate(inputamp):
        outputamp.append(sumamp * fac - iamp)
    inputarr = np.array(inputamp)
    print(outputamp, " ", np.sum(outputamp), " ", np.sum(inputarr * inputarr))
    return outputamp
      
# used for testing only     
def arbitrary_2penalty_input_output(inputamp, inputampstar, bFermion = False):
    """
    Penalize for not preserving particle number;
    (possibly) penalize for fermionicity
    """
    Narc = len(inputamp)
    fac = 2.0 / Narc
    
    outputamp = np.ones((Narc,)).astype("int")
    outputampstar = np.ones((Narc,)).astype("int")
    for (iinpvec, ioutpvec) in zip((inputamp, inputampstar), (outputamp, outputampstar)):
      ioutpvec = -iinpvec # the huygens component
      for iivec, ivec in enumerate(iinpvec):
          absivec = np.abs(ivec)
          sgnivec = np.sgn(ivec)
          for idraw in range(absivec):
              #once
              iidraw = np.random.randint(0,Narc)
              ioutpvec[iidraw] += sgnivec
              #twice
              iidraw = np.random.randint(0,Narc)
              ioutpvec[iidraw] += sgnivec
              
    # the above has reproduced statistically the wave equation with amplitudes exactly preserved; now, let's assess the errors
    
    
        
    inputarr = np.array(inputamp)
    print(outputamp, " ", np.sum(outputamp), " ", np.sum(inputarr * inputarr))
    return outputamp
      
# slower than classicalwave() but it preserves amplitube at each node   
def wavegridupdate(t, grid, #fermpenalty, n2penalty, bFermion = False, 
                              Narc=4):
    """
    assume that even parity cells will be updated
    """
    #import pdb; pdb.set_trace()
    def diffupdate(amp, grid, iscen, ix, iy):
        absamp = np.abs(amp)
        sgnamp = np.sign(amp)
        for i in range(absamp):
            
            r1 = np.random.randint(0,Narc)
            
            bNonIdenticalBrownian = False # set to true if you want Brownian particles to always pick different output arcs
            bChooseOpposite = True
            
            if bNonIdenticalBrownian:
              r2 = np.random.randint(0,Narc-1)
              if r2 >= r1:
                r2 += 1
            if bChooseOpposite:
              r2 = (r1 + 2) % 4
            
            grid[iscen, ix, iy, r1] += sgnamp
            grid[iscen, ix, iy, r2] += sgnamp
            
            
            
    for iscen in range(grid.shape[0]):
        for ix in range(grid.shape[1]):
            for iy in range(grid.shape[2]):
                if (ix + iy + t) % 2 == 0:
                    continue
                    
                    

                if ix == 0:
                    xminusoffset = grid.shape[1] -1
                else:
                    xminusoffset = ix - 1
                if ix >= grid.shape[1] - 1:
                    xplusoffset = 0
                else:
                    xplusoffset = ix + 1
                if iy == 0:
                    yminusoffset = grid.shape[2] -1
                else:
                    yminusoffset = iy - 1
                if iy >= grid.shape[2] - 1:
                    yplusoffset = 0
                else:
                    yplusoffset = iy + 1
                    
                    
                sumamp = grid[iscen,xplusoffset,iy,2] + grid[iscen,ix,yplusoffset,3] + grid[iscen,xminusoffset,iy,0] + grid[iscen,ix,yminusoffset,1]
                    
                #if (iscen, ix, iy)==(0, 1, 12):
                #    print("whoa")
                #    import pdb; pdb.set_trace()

                
                #if (abs(grid[iscen,xplusoffset,iy,2]) + 
                #    abs(grid[iscen,ix,yplusoffset,3]) +
                #    abs(grid[iscen,xminusoffset,iy,0]) +
                #    abs(grid[iscen,ix,yminusoffset,1])) != 0:
                #    import pdb; pdb.set_trace()
                    
                    
                    
                # add the Huygens contribution first
                grid[iscen,ix,iy,0] = -1 * grid[iscen,xplusoffset,iy,2]
                grid[iscen,ix,iy,1] = -1 * grid[iscen,ix,yplusoffset,3]
                grid[iscen,ix,iy,2] = -1 * grid[iscen,xminusoffset,iy,0]
                grid[iscen,ix,iy,3] = -1 * grid[iscen,ix,yminusoffset,1]
                
                #import pdb; pdb.set_trace()
                diffupdate(grid[iscen, xplusoffset, iy, 2],  grid, iscen, ix, iy)
                diffupdate(grid[iscen, ix, yplusoffset, 3],  grid, iscen, ix, iy)
                diffupdate(grid[iscen, xminusoffset, iy, 0], grid, iscen, ix, iy)
                diffupdate(grid[iscen, ix, yminusoffset, 1], grid, iscen, ix, iy)



def  ProcessFile():
  
  global bUseLatticeGasForRand
  global bCplusplus
  opts, args = ReadParams()
  #import pdb; pdb.set_trace()
  seed = abs( opts.seed )
  
  np.random.seed( seed  )
  bDiscrete = opts.bDiscrete
  
  nplotmod = 1000 #20000 #1200000 
  if not(bDiscrete):
      nplotmod = 500 # 

  # Deprecated: it's possible to process entire i,j portion of grid at once, which is much faster
  def propagate(grid,  i, j):
    grid[:, i,j,Comp000] = grid[:, i-1,j,Comp000]
    grid[:, i,j,Comp090] = grid[:, i,j-1,Comp090]
    grid[:, i,j,Comp180] = grid[:, (i+1) % grid.shape[1],j,Comp180]
    grid[:, i,j,Comp270] = grid[:, i,(j+1) % grid.shape[2],Comp270]
  
  
  
  # Deprecated
  def propagate(grid,  i, j):
    xup = np.zeros(grid.shape[:3])
    xup[:, 1:, :] = grid[:, :-1, :, Comp000]
    xup[:, 0, :] = grid[:,-1, :, Comp000]
    
    yup = np.zeros(grid.shape[:3])
    yup[:, :, 1:] = grid[:, :, :-1, Comp090]
    yup[:, :, 0] = grid[:, :, -1, Comp090]   
      
    xdn = np.zeros(grid.shape[:3])
    xdn[:, :-1, :] = grid[:, 1:, :, Comp180]
    xdn[:, -1, :] = grid[:, 0, :, Comp180]

    ydn = np.zeros(grid.shape[:3])
    ydn[:, :, :-1] = grid[:, :, 1:, Comp270]
    ydn[:, :, -1] = grid[:, :, 0, Comp270]
            



    grid[:, i,j,Comp000] = xup
    grid[:, i,j,Comp090] = yup
    grid[:, i,j,Comp180] = xdn
    grid[:, i,j,Comp270] = ydn
  
  
  
  
  # deprecated: it's possible to process entire i,j portion of grid at once, which is much faster
  def propagate_rotate90oncollision(grid, i_scenario, i, j):
    # same as above, except that a head-on collision causes outgoing particles to rotate
    # if there's more in one colliding flow than another, the excess will propagate as before  
    if np.sum(np.array([grid[i_scenario, i-1,j,Comp000], grid[i_scenario, i,j-1,Comp090], grid[i_scenario, (i+1) % grid.shape[1],j,Comp180], grid[i_scenario, i,(j+1) % grid.shape[2],Comp270]]) != 0) <= 2:
      new000 = grid[i_scenario, i-1,j,Comp000] + min([grid[i_scenario, i,j-1,Comp090], grid[i_scenario, i,(j+1) % grid.shape[2],Comp270]]) - min([grid[i_scenario, i-1,j,Comp000], grid[i_scenario, (i+1) % grid.shape[1],j,Comp180]])
      new090 = grid[i_scenario, i,j-1,Comp090] + min([grid[i_scenario, i-1,j,Comp000], grid[i_scenario, (i+1) % grid.shape[1],j,Comp180]]) - min([grid[i_scenario, i,j-1,Comp090], grid[i_scenario, i,(j+1) % grid.shape[2],Comp270]])
      new180 = grid[i_scenario, (i+1) % grid.shape[1],j,Comp180] + min([grid[i_scenario, i,j-1,Comp090], grid[i_scenario, i,(j+1) % grid.shape[2],Comp270]]) - min([grid[i_scenario, i-1,j,Comp000], grid[i_scenario, (i+1) % grid.shape[1],j,Comp180]]) 
      new270 = grid[i_scenario, i,(j+1) % grid.shape[2],Comp270] + min([grid[i_scenario, i-1,j,Comp000], grid[i_scenario, (i+1) % grid.shape[1],j,Comp180]]) - min([grid[i_scenario, i,j-1,Comp090], grid[i_scenario, i,(j+1) % grid.shape[2],Comp270]])
      
      if np.min(np.array([new000, new090, new180, new270])) < 0 or np.max(np.array([new000, new090, new180, new270])) > 1:
        import pdb; pdb.set_trace()
      
      grid[i_scenario, i,j,Comp000] = copy.copy(new000)
      grid[i_scenario, i,j,Comp090] = copy.copy(new090)
      grid[i_scenario, i,j,Comp180] = copy.copy(new180)
      grid[i_scenario, i,j,Comp270] = copy.copy(new270)
    else:
      new000 = 1-grid[i_scenario, i-1,j,Comp000] + min([1-grid[i_scenario, i,j-1,Comp090], 1-grid[i_scenario, i,(j+1) % grid.shape[2],Comp270]]) - min([1-grid[i_scenario, i-1,j,Comp000], 1-grid[i_scenario, (i+1) % grid.shape[1],j,Comp180]])
      new090 = 1-grid[i_scenario, i,j-1,Comp090] + min([1-grid[i_scenario, i-1,j,Comp000], 1-grid[i_scenario, (i+1) % grid.shape[1],j,Comp180]]) - min([1-grid[i_scenario, i,j-1,Comp090], 1-grid[i_scenario, i,(j+1) % grid.shape[2],Comp270]])
      new180 = 1-grid[i_scenario, (i+1) % grid.shape[0],j,Comp180] + min([1-grid[i_scenario, i,j-1,Comp090], 1-grid[i_scenario, i,(j+1) % grid.shape[2],Comp270]]) - min([1-grid[i_scenario, i-1,j,Comp000], 1-grid[i_scenario, (i+1) % grid.shape[1],j,Comp180]]) 
      new270 = 1-grid[i_scenario, i,(j+1) % grid.shape[1],Comp270] + min([1-grid[i_scenario, i-1,j,Comp000], 1-grid[i_scenario, (i+1) % grid.shape[1],j,Comp180]]) - min([1-grid[i_scenario, i,j-1,Comp090], 1-grid[i_scenario, i,(j+1) % grid.shape[2],Comp270]])
      
      #if np.min(np.array([new000, new090, new180, new270])) < 0 and np.max(np.array([new000, new090, new180, new270])) > 1:
      #  import pdb; pdb.set_trace()


      grid[i_scenario, i,j,Comp000] = 1-new000
      grid[i_scenario, i,j,Comp090] = 1-new090
      grid[i_scenario, i,j,Comp180] = 1-new180
      grid[i_scenario, i,j,Comp270] = 1-new270




  # deprecated: it's possible to process entire i,j portion of grid at once, which is much faster
  def propagate_rotate90oncollision_multi(grid, i, j):
      global NTEST
    # same as above, except that a head-on collision causes outgoing particles to rotate
    # if there's more in one colliding flow than another, the excess will propagate as before  
      xup = grid[:, i-1,j,Comp000] # not the name xup indicates the flow is coming FROM the xDOWN direction (and heading in the xup direction)
      yup = grid[:, i,j-1,Comp090]
      xdn = grid[:, (i+1) % grid.shape[1],j,Comp180]
      ydn = grid[:, i,(j+1) % grid.shape[2],Comp270]
      
      nt = np.logical_not
      an = np.logical_and
      #import pdb; pdb.set_trace()
      #               
      
      #xup,yup,xdn,ydn
      #new000.astype(int),new090.astype(int),new180.astype(int),new270.astype(int)
      # xup-xdn-new000.astype(int)+new180.astype(int)
      # yup-ydn-new090.astype(int)+new270.astype(int)
      #          1                                        2-oppose                               3-nonsim                        3-accrossA                         3-acrossB                        4-fold                             L-shapeA                               L-shapeB
      new000 = an(an(an(xup,nt(yup)),nt(xdn)),nt(ydn)) + an(an(an(nt(xup),yup),nt(xdn)),ydn) + an(an(an(xup,yup),nt(xdn)),ydn) + an(an(an(xup,nt(yup)),xdn),ydn) + an(an(an(xup,yup),xdn),nt(ydn))   + an(an(an(xup,yup),xdn),ydn) + an(an(an(xup,yup),nt(xdn)),nt(ydn)) + an(an(an(xup,nt(yup)),nt(xdn)),ydn)   
      new090 = an(an(an(nt(xup),yup),nt(xdn)),nt(ydn)) + an(an(an(xup,nt(yup)),xdn),nt(ydn)) + an(an(an(xup,yup),xdn),nt(ydn)) + an(an(an(nt(xup),yup),xdn),ydn) + an(an(an(xup,yup),nt(xdn)),ydn)   + an(an(an(xup,yup),xdn),ydn) + an(an(an(xup,yup),nt(xdn)),nt(ydn)) + an(an(an(nt(xup),yup),xdn),nt(ydn))     
      new180 = an(an(an(nt(xup),nt(yup)),xdn),nt(ydn)) + an(an(an(nt(xup),yup),nt(xdn)),ydn) + an(an(an(nt(xup),yup),xdn),ydn) + an(an(an(xup,nt(yup)),xdn),ydn) + an(an(an(xup,yup),xdn),nt(ydn))   + an(an(an(xup,yup),xdn),ydn) + an(an(an(nt(xup),nt(yup)),xdn),ydn) + an(an(an(nt(xup),yup),xdn),nt(ydn))    
      new270 = an(an(an(nt(xup),nt(yup)),nt(xdn)),ydn) + an(an(an(xup,nt(yup)),xdn),nt(ydn)) + an(an(an(xup,nt(yup)),xdn),ydn) + an(an(an(nt(xup),yup),xdn),ydn) + an(an(an(xup,yup),nt(xdn)),ydn)   + an(an(an(xup,yup),xdn),ydn) + an(an(an(nt(xup),nt(yup)),xdn),ydn) + an(an(an(xup,nt(yup)),nt(xdn)),ydn) 
      
      #if not(np.array_equal((new000.astype(int) + new090.astype(int) + new180.astype(int) + new270.astype(int)),(xup + yup + xdn + ydn))):
      #  import pdb; pdb.set_trace()
          
      #if np.sum((new000.astype(int) - xup) + (new180.astype(int) - xdn)) != 0:
      #  import pdb; pdb.set_trace()
      #    
      #if np.sum((new090.astype(int) - yup) + (new270.astype(int) - ydn)) != 0:
      #  import pdb; pdb.set_trace()
          
      #if np.sum([new000, new090*2, new180*4, new270*8]) == NTEST:
      #  import pdb; pdb.set_trace()
      #  NTEST += 1
      #  print "NTEST " + str(NTEST)
        
        
      if np.sum(np.abs(yup-ydn-new090.astype(int)+new270.astype(int))) != 0 or \
           np.sum(np.abs(yup-ydn-new090.astype(int)+new270.astype(int))) != 0:
        import pdb; pdb.set_trace()
        print2("mismatch")
        
        
      grid[:, i,j,Comp000] = new000.astype(int)
      grid[:, i,j,Comp090] = new090.astype(int)
      grid[:, i,j,Comp180] = new180.astype(int)
      grid[:, i,j,Comp270] = new270.astype(int)







  # use this only if you want a simple and easy lattice-gas to drive the randomness
  def propagate_rotate90oncollision_ONESWEEP(grid):
      global NTEST
    # same as above, except that a head-on collision causes outgoing particles to rotate
    # if there's more in one colliding flow than another, the excess will propagate as before  
      #xup = grid[:, i-1,j,Comp000] # not the name xup indicates the flow is coming FROM the xDOWN direction (and heading in the xup direction)
      #yup = grid[:, i,j-1,Comp090]
      #xdn = grid[:, (i+1) % grid.shape[1],j,Comp180]
      #ydn = grid[:, i,(j+1) % grid.shape[2],Comp270]
      
      
      
      xup = np.zeros(grid.shape[:3])
      xup[:, 1:, :] = grid[:, :-1, :, Comp000]
      xup[:, 0, :] = grid[:,-1, :, Comp000]
      
      yup = np.zeros(grid.shape[:3])
      yup[:, :, 1:] = grid[:, :, :-1, Comp090]
      yup[:, :, 0] = grid[:, :, -1, Comp090]   
        
      xdn = np.zeros(grid.shape[:3])
      xdn[:, :-1, :] = grid[:, 1:, :, Comp180]
      xdn[:, -1, :] = grid[:, 0, :, Comp180]

      ydn = np.zeros(grid.shape[:3])
      ydn[:, :, :-1] = grid[:, :, 1:, Comp270]
      ydn[:, :, -1] = grid[:, :, 0, Comp270]
              
              
      
      
      
      nt = np.logical_not
      an = np.logical_and
      #import pdb; pdb.set_trace()
      #               
      
      #xup,yup,xdn,ydn
      #new000.astype(int),new090.astype(int),new180.astype(int),new270.astype(int)
      # xup-xdn-new000.astype(int)+new180.astype(int)
      # yup-ydn-new090.astype(int)+new270.astype(int)
      #          1                                        2-oppose                               3-nonsim                        3-accrossA                         3-acrossB                        4-fold                             L-shapeA                               L-shapeB
      new000 = an(an(an(xup,nt(yup)),nt(xdn)),nt(ydn)) + an(an(an(nt(xup),yup),nt(xdn)),ydn) + an(an(an(xup,yup),nt(xdn)),ydn) + an(an(an(xup,nt(yup)),xdn),ydn) + an(an(an(xup,yup),xdn),nt(ydn))   + an(an(an(xup,yup),xdn),ydn) + an(an(an(xup,yup),nt(xdn)),nt(ydn)) + an(an(an(xup,nt(yup)),nt(xdn)),ydn)   
      new090 = an(an(an(nt(xup),yup),nt(xdn)),nt(ydn)) + an(an(an(xup,nt(yup)),xdn),nt(ydn)) + an(an(an(xup,yup),xdn),nt(ydn)) + an(an(an(nt(xup),yup),xdn),ydn) + an(an(an(xup,yup),nt(xdn)),ydn)   + an(an(an(xup,yup),xdn),ydn) + an(an(an(xup,yup),nt(xdn)),nt(ydn)) + an(an(an(nt(xup),yup),xdn),nt(ydn))     
      new180 = an(an(an(nt(xup),nt(yup)),xdn),nt(ydn)) + an(an(an(nt(xup),yup),nt(xdn)),ydn) + an(an(an(nt(xup),yup),xdn),ydn) + an(an(an(xup,nt(yup)),xdn),ydn) + an(an(an(xup,yup),xdn),nt(ydn))   + an(an(an(xup,yup),xdn),ydn) + an(an(an(nt(xup),nt(yup)),xdn),ydn) + an(an(an(nt(xup),yup),xdn),nt(ydn))    
      new270 = an(an(an(nt(xup),nt(yup)),nt(xdn)),ydn) + an(an(an(xup,nt(yup)),xdn),nt(ydn)) + an(an(an(xup,nt(yup)),xdn),ydn) + an(an(an(nt(xup),yup),xdn),ydn) + an(an(an(xup,yup),nt(xdn)),ydn)   + an(an(an(xup,yup),xdn),ydn) + an(an(an(nt(xup),nt(yup)),xdn),ydn) + an(an(an(xup,nt(yup)),nt(xdn)),ydn) 
      
      #if not(np.array_equal((new000.astype(int) + new090.astype(int) + new180.astype(int) + new270.astype(int)),(xup + yup + xdn + ydn))):
      #  import pdb; pdb.set_trace()
      #    
      #if not(np.array_equal(np.sum(grid,(1,2,3)),np.sum(xup + yup + xdn + ydn,(1,2)))):
      #  import pdb; pdb.set_trace()
      #    
      #if not(np.array_equal(np.sum(new000.astype(int) + new090.astype(int) + new180.astype(int) + new270.astype(int), (1,2)),np.sum(grid, (1,2,3)))):
      #  import pdb; pdb.set_trace()
          
          
          
          
          
          
      #if np.sum((new000.astype(int) - xup) + (new180.astype(int) - xdn)) != 0:
      #  import pdb; pdb.set_trace()
      #    
      #if np.sum((new090.astype(int) - yup) + (new270.astype(int) - ydn)) != 0:
      #  import pdb; pdb.set_trace()
          
      #if np.sum([new000, new090*2, new180*4, new270*8]) == NTEST:
      #  import pdb; pdb.set_trace()
      #  NTEST += 1
      #  print "NTEST " + str(NTEST)
        
        
      if np.sum(np.abs(yup-ydn-new090.astype(int)+new270.astype(int))) != 0 or \
           np.sum(np.abs(yup-ydn-new090.astype(int)+new270.astype(int))) != 0:
        import pdb; pdb.set_trace()
        print2("mismatch")
        
         
      grid[:, :, :, Comp000] = new000.astype(int)
      grid[:, :, :, Comp090] = new090.astype(int)
      grid[:, :, :, Comp180] = new180.astype(int)
      grid[:, :, :, Comp270] = new270.astype(int)










  # not used currently
  def brownian(grid):
      global NTEST
      
      randgrid = np.random.random_sample(grid.shape[:-1]) * 4
      
      randgrid = randgrid.astype(int)
        
      #xup = grid[:, i-1,j,Comp000] # not the name xup indicates the flow is coming FROM the xDOWN direction (and heading in the xup direction)
      
      xup = np.zeros(grid.shape[:3])
      xup[:, 1:, :] = grid[:, :-1, :, Comp000]
      xup[:, 0, :] = grid[:,-1, :, Comp000]
      
      yup = np.zeros(grid.shape[:3])
      yup[:, :, 1:] = grid[:, :, :-1, Comp090]
      yup[:, :, 0] = grid[:, :, -1, Comp090]   
        
      xdn = np.zeros(grid.shape[:3])
      xdn[:, :-1, :] = grid[:, 1:, :, Comp180]
      xdn[:, -1, :] = grid[:, 0, :, Comp180]

      ydn = np.zeros(grid.shape[:3])
      ydn[:, :, :-1] = grid[:, :, 1:, Comp270]
      ydn[:, :, -1] = grid[:, :, 0, Comp270]
              
              
      #import pdb; pdb.set_trace()       
      new000 = xup * (randgrid == Comp000) + yup * (randgrid == Comp000) + xdn * (randgrid == Comp000) + ydn * (randgrid == Comp000)
      new090 = xup * (randgrid == Comp090) + yup * (randgrid == Comp090) + xdn * (randgrid == Comp090) + ydn * (randgrid == Comp090)
      new180 = xup * (randgrid == Comp180) + yup * (randgrid == Comp180) + xdn * (randgrid == Comp180) + ydn * (randgrid == Comp180)
      new270 = xup * (randgrid == Comp270) + yup * (randgrid == Comp270) + xdn * (randgrid == Comp270) + ydn * (randgrid == Comp270) 
      
      #if not(np.array_equal((new000.astype(int) + new090.astype(int) + new180.astype(int) + new270.astype(int)),(xup + yup + xdn + ydn))):
      #  import pdb; pdb.set_trace()
          
      #if np.sum((new000.astype(int) - xup) + (new180.astype(int) - xdn)) != 0:
      #  import pdb; pdb.set_trace()
      #    
      #if np.sum((new090.astype(int) - yup) + (new270.astype(int) - ydn)) != 0:
      #  import pdb; pdb.set_trace()
          
      #if np.sum([new000, new090*2, new180*4, new270*8]) == NTEST:
      #  import pdb; pdb.set_trace()
      #  NTEST += 1
      #  print "NTEST " + str(NTEST)
        
        
       
        
      grid[:, :, :, Comp000] = new000.astype(int)
      grid[:, :, :, Comp090] = new090.astype(int)
      grid[:, :, :, Comp180] = new180.astype(int)
      grid[:, :, :, Comp270] = new270.astype(int)
      




  # not used currently
  def poissonbrownian(grid, probshuffle):
      # random shuffle will take place probabilistically
      # i.e. at each point, a random number is drawn, and if it "fails", particles
      # there will simply propagate as in propagate()
      
      global NTEST
      
      randgrid = np.random.random_sample(grid.shape[:-1]) * 4
        
      randgrid_shuffleprob = (np.random.random_sample(grid.shape[:-1]) < probshuffle).astype(int)
        
      randgrid = randgrid.astype(int)
        
      #xup = grid[:, i-1,j,Comp000] # not the name xup indicates the flow is coming FROM the xDOWN direction (and heading in the xup direction)
      
      xup = np.zeros(grid.shape[:3])
      xup[:, 1:, :] = grid[:, :-1, :, Comp000]
      xup[:, 0, :] = grid[:,-1, :, Comp000]
      
      yup = np.zeros(grid.shape[:3])
      yup[:, :, 1:] = grid[:, :, :-1, Comp090]
      yup[:, :, 0] = grid[:, :, -1, Comp090]   
        
      xdn = np.zeros(grid.shape[:3])
      xdn[:, :-1, :] = grid[:, 1:, :, Comp180]
      xdn[:, -1, :] = grid[:, 0, :, Comp180]

      ydn = np.zeros(grid.shape[:3])
      ydn[:, :, :-1] = grid[:, :, 1:, Comp270]
      ydn[:, :, -1] = grid[:, :, 0, Comp270]
              
              
      #import pdb; pdb.set_trace()       
      new000 = xup * (randgrid == Comp000) + yup * (randgrid == Comp000) + xdn * (randgrid == Comp000) + ydn * (randgrid == Comp000)
      new090 = xup * (randgrid == Comp090) + yup * (randgrid == Comp090) + xdn * (randgrid == Comp090) + ydn * (randgrid == Comp090)
      new180 = xup * (randgrid == Comp180) + yup * (randgrid == Comp180) + xdn * (randgrid == Comp180) + ydn * (randgrid == Comp180)
      new270 = xup * (randgrid == Comp270) + yup * (randgrid == Comp270) + xdn * (randgrid == Comp270) + ydn * (randgrid == Comp270) 
      
      #if not(np.array_equal((new000.astype(int) + new090.astype(int) + new180.astype(int) + new270.astype(int)),(xup + yup + xdn + ydn))):
      #  import pdb; pdb.set_trace()
          
      #if np.sum((new000.astype(int) - xup) + (new180.astype(int) - xdn)) != 0:
      #  import pdb; pdb.set_trace()
      #    
      #if np.sum((new090.astype(int) - yup) + (new270.astype(int) - ydn)) != 0:
      #  import pdb; pdb.set_trace()
          
      #if np.sum([new000, new090*2, new180*4, new270*8]) == NTEST:
      #  import pdb; pdb.set_trace()
      #  NTEST += 1
      #  print "NTEST " + str(NTEST)
        

       
        
      grid[:, :, :, Comp000] = new000.astype(int) * randgrid_shuffleprob + xup * (1-randgrid_shuffleprob)
      grid[:, :, :, Comp090] = new090.astype(int) * randgrid_shuffleprob + yup * (1-randgrid_shuffleprob)
      grid[:, :, :, Comp180] = new180.astype(int) * randgrid_shuffleprob + xdn * (1-randgrid_shuffleprob)
      grid[:, :, :, Comp270] = new270.astype(int) * randgrid_shuffleprob + ydn * (1-randgrid_shuffleprob)
      






  
  def classicalwavefloatingpt(grid, N_Per_Arc=1):
      global NTEST
      # same as above, except that a head-on collision causes outgoing particles to rotate
      # if there's more in one colliding flow than another, the excess will propagate as before  
      xup = np.zeros(grid.shape[:3])
      xup[:, 1:, :] = grid[:, :-1, :, Comp000]
      xup[:, 0, :] = grid[:,-1, :, Comp000]
      
      yup = np.zeros(grid.shape[:3])
      yup[:, :, 1:] = grid[:, :, :-1, Comp090]
      yup[:, :, 0] = grid[:, :, -1, Comp090]   
        
      xdn = np.zeros(grid.shape[:3])
      xdn[:, :-1, :] = grid[:, 1:, :, Comp180]
      xdn[:, -1, :] = grid[:, 0, :, Comp180]

      ydn = np.zeros(grid.shape[:3])
      ydn[:, :, :-1] = grid[:, :, 1:, Comp270]
      ydn[:, :, -1] = grid[:, :, 0, Comp270]
      
      a = 1 / (2.0 * N_Per_Arc) # note this is TWICE the "heat-equation" kernel, for which a flow of X on one incoming arc, puts an outgoing flow of X/N_arc for all the outgoing arcs; A in the above equation is always TWICE that.
      
      # the wave formula is the above factor for every "non-same" arc, whereas for the identical arc, it is (-1 + a)
      
      #import pdb; pdb.set_trace() 
      
      #bEmptyWave = True
      #if bEmptyWave:
      #  a = 0.25
      
      #for k in range(N_Per_Arc):
      #  new000[i, j, k] =   1
      
           
      new000 = (xup + yup + xdn + ydn) * a - xdn
      new090 = (xup + yup + xdn + ydn) * a - ydn
      new180 = (xup + yup + xdn + ydn) * a - xup
      new270 = (xup + yup + xdn + ydn) * a - yup 
      
      
      # first, sum the amplitudes. If it's bosonic, then output will be the truncated integer, plus a random extra dictated by the non-integer part


       
        
      grid[:, :, :, Comp000] = new000  #.astype(int) * randgrid_shuffleprob + xup * (1-randgrid_shuffleprob)
      grid[:, :, :, Comp090] = new090  #.astype(int) * randgrid_shuffleprob + yup * (1-randgrid_shuffleprob)
      grid[:, :, :, Comp180] = new180  #.astype(int) * randgrid_shuffleprob + xdn * (1-randgrid_shuffleprob)
      grid[:, :, :, Comp270] = new270  #.astype(int) * randgrid_shuffleprob + ydn * (1-randgrid_shuffleprob)
      


  #def cythonclassicalwave(oldgrid,  newgrid, Nscen, NDimX, NDimY, NDir, randlong):


  #INCOMPLETE -- note this does not preserve amplitude (except "on average")
  def classicalwave(grid,  gridlatticegas, N_Per_Arc=1):
      global NTEST
      global bUseLatticeGasForRand
      
      
      
  
      # same as above, except that a head-on collision causes outgoing particles to rotate
      # if there's more in one colliding flow than another, the excess will propagate as before  
      xup = np.zeros(grid.shape[:3])
      xup[:, 1:, :] = grid[:, :-1, :, Comp000]
      xup[:, 0, :] = grid[:,-1, :, Comp000]
      
      yup = np.zeros(grid.shape[:3])
      yup[:, :, 1:] = grid[:, :, :-1, Comp090]
      yup[:, :, 0] = grid[:, :, -1, Comp090]   
        
      xdn = np.zeros(grid.shape[:3])
      xdn[:, :-1, :] = grid[:, 1:, :, Comp180]
      xdn[:, -1, :] = grid[:, 0, :, Comp180]

      ydn = np.zeros(grid.shape[:3])
      ydn[:, :, :-1] = grid[:, :, 1:, Comp270]
      ydn[:, :, -1] = grid[:, :, 0, Comp270]
      
      a = 1 / (2.0 * N_Per_Arc) # note this is TWICE the "heat-equation" kernel, for which a flow of X on one incoming arc, puts an outgoing flow of X/N_arc for all the outgoing arcs; A in the above equation is always TWICE that.
      
      
      bPreserveGenParticleNumber = False # the sum of bra*ket 
      if bPreserveGenParticleNumber: # the sum of bra*ket 
          incomingbraketlist = list(copy.copy(grid.shape))
          incomingbraketlist[0] //= 2
          incomingbraketlist = incomingbraketlist[:-1]
          incomingbraketlist = tuple(incomingbraketlist)
          
          #import pdb; pdb.set_trace()
          incomingbraket = np.zeros(incomingbraketlist)
          for i in range(0, grid.shape[0], 2):
              #import pdb; pdb.set_trace()
              incomingbraket[i//2,:,:] = xup[i,:,:] * xup[i+1,:,:] + yup[i,:,:] * yup[i+1,:,:] + xdn[i,:,:] * xdn[i+1,:,:] + ydn[i,:,:] * ydn[i+1,:,:]
          
      
      
      
      
      
      
      
      
      # the wave formula is the above factor for every "non-same" arc, whereas for the identical arc, it is (-1 + a)
      
      #import pdb; pdb.set_trace() 
      
      
      
      #for k in range(N_Per_Arc):
      #  new000[i, j, k] =   1
      
      a = 1 / (2.0 * N_Per_Arc) # note this is TWICE the "heat-equation" kernel, for which a flow of X on one incoming arc, puts an outgoing flow of X/N_arc for all the outgoing arcs; A in the above equation is always TWICE that.
      
      # the wave formula is the above factor for every "non-same" arc, whereas for the identical arc, it is (-1 + a)
      
      #import pdb; pdb.set_trace() 
      
      
      
      #for k in range(N_Per_Arc):
      #  new000[i, j, k] =   1
      
        
      new000 = (xup + yup + xdn + ydn) * a - xdn
      new090 = (xup + yup + xdn + ydn) * a - ydn
      new180 = (xup + yup + xdn + ydn) * a - xup
      new270 = (xup + yup + xdn + ydn) * a - yup



      # first, sum the amplitudes. If it's bosonic, then output will be the truncated integer, plus a random extra dictated by the non-integer part
      
      #import pdb; pdb.set_trace()
      #                   ~~~~~this product is non-zero only where LHS is non-zero
      #                                             ~~~~~~this term makes the addend be only one of the values [-1/2,0,1/2]
      
      if bUseLatticeGasForRand:
          new000 = (new000 + (new000 * gridlatticegas[:,:,:,Comp000] - np.trunc(new000))).astype("int")
          new090 = (new090 + (new090 * gridlatticegas[:,:,:,Comp090] - np.trunc(new090))).astype("int")
          new180 = (new180 + (new180 * gridlatticegas[:,:,:,Comp180] - np.trunc(new180))).astype("int")
          new270 = (new270 + (new270 * gridlatticegas[:,:,:,Comp270] - np.trunc(new270))).astype("int")
      else:
      
          #import pdb; pdb.set_trace()
          bFlip = False
          if np.random.random_sample() < 0.5:
              bFlip = True
          if bFlip:
              new000 = ((-np.round(-(new000 + 0.01*(np.random.random_sample(new000.shape)-0.5)))).astype("int"))  #      (new000 * gridlatticegas[:,:,:,Comp000] - np.trunc(new000))).astype("int")
              new090 = ((-np.round(-(new090 + 0.01*(np.random.random_sample(new090.shape)-0.5)))).astype("int"))  #      (new090 + (new090 * gridlatticegas[:,:,:,Comp090] - np.trunc(new090))).astype("int")
              new180 = ((-np.round(-(new180 + 0.01*(np.random.random_sample(new180.shape)-0.5)))).astype("int"))  #      (new180 + (new180 * gridlatticegas[:,:,:,Comp180] - np.trunc(new180))).astype("int")
              new270 = ((-np.round(-(new270 + 0.01*(np.random.random_sample(new270.shape)-0.5)))).astype("int"))  #      (new270 + (new270 * gridlatticegas[:,:,:,Comp270] - np.trunc(new270))).astype("int")
          
          else:
              new000 = (np.round(new000 + 0.01*(np.random.random_sample(new000.shape)-0.5))).astype("int")  #      (new000 * gridlatticegas[:,:,:,Comp000] - np.trunc(new000))).astype("int")
              new090 = (np.round(new090 + 0.01*(np.random.random_sample(new090.shape)-0.5))).astype("int")  #      (new090 + (new090 * gridlatticegas[:,:,:,Comp090] - np.trunc(new090))).astype("int")
              new180 = (np.round(new180 + 0.01*(np.random.random_sample(new180.shape)-0.5))).astype("int")  #      (new180 + (new180 * gridlatticegas[:,:,:,Comp180] - np.trunc(new180))).astype("int")
              new270 = (np.round(new270 + 0.01*(np.random.random_sample(new270.shape)-0.5))).astype("int")  #      (new270 + (new270 * gridlatticegas[:,:,:,Comp270] - np.trunc(new270))).astype("int")
      
      
      
        
      grid[:, :, :, Comp000] = new000  #.astype(int) * randgrid_shuffleprob + xup * (1-randgrid_shuffleprob)
      grid[:, :, :, Comp090] = new090  #.astype(int) * randgrid_shuffleprob + yup * (1-randgrid_shuffleprob)
      grid[:, :, :, Comp180] = new180  #.astype(int) * randgrid_shuffleprob + xdn * (1-randgrid_shuffleprob)
      grid[:, :, :, Comp270] = new270  #.astype(int) * randgrid_shuffleprob + ydn * (1-randgrid_shuffleprob)
      
      
      
      
      if bPreserveGenParticleNumber:
          outgoingbraketlist = list(copy.copy(grid.shape))
          outgoingbraketlist[0] //= 2
          outgoingbraketlist = outgoingbraketlist[:-1]
          outgoingbraketlist = tuple(outgoingbraketlist)
          outgoingbraket = np.zeros(outgoingbraketlist)
          for i in range(0, grid.shape[0], 2):
              #import pdb; pdb.set_trace()
              outgoingbraket[i//2,:,:] = np.sum(grid[i,:,:,:] * grid[i+1,:,:,:], 2)
              #import pdb; pdb.set_trace()
              for ibk in range(outgoingbraket.shape[0]):
                  for ix in range(outgoingbraket.shape[1]):
                      for iy in range(outgoingbraket.shape[2]):
                          if outgoingbraket[ibk, ix,iy] != incomingbraket[ibk, ix,iy]:
                              import pdb; pdb.set_trace()
                              
                          
          
      
      
      
      
      
      










      
  def Flatten(grid):
    tmp = 0
    grid2 = np.zeros((grid.shape[1]*2, grid.shape[2]*2))
    for i in range(grid.shape[1]):
      for j in range(grid.shape[2]):
        grid2[2*i, 2*j] = np.mean(grid[:, i,j,Comp000])
        grid2[2*i, 2*j+1] = np.mean(grid[:, i,j,Comp090])
        grid2[2*i+1, 2*j+1] = np.mean(grid[:, i,j,Comp180])
        grid2[2*i+1, 2*j] = np.mean(grid[:, i,j,Comp270])
    #import pdb; pdb.set_trace()
    return grid2
    
    
  #deprecated
  def FlattenSUB(grid, grid2):
    tmp = 0
    grid2 = np.zeros(grid2.shape)
    for i in range(grid.shape[0]):
      for j in range(grid.shape[1]):
        grid2[2*i, 2*j] = grid[i,j,Comp000]
        grid2[2*i, 2*j+1] = grid[i,j,Comp090]
        grid2[2*i+1, 2*j+1] = grid[i,j,Comp180]
        grid2[2*i+1, 2*j] = grid[i, j, Comp270]
    return grid2
    
    
    
    
    
  print("Start 2d systolic lattice gas: 2 opposing incoming flows create 2 opposing outgoing flows rotated by 90 degrees. ")
  print("Since this is systolic, there's 2 independent lattices (such that X+Y+T has either odd or even parity) ")  
  print("At each time step, we will calculate the OUTGOING flows emanaging from each point (as opposed to the incoming flows converging on a point. ")

  print("Depending on different settings wc can use this to generate continuous 'floating-pt' amplitudes, but also discrete waves that converge to the continuous distribution only with a large enough set of runs; in the latter case, the randomness can be driven not just by np.random, but also by a simple lattice-gas/IsingModel particle system.")
	
	# declare the 2-d lattice (that has 2*2 outgoing flows at each point.
  NDIM = 16 # 64 # 100 # must be at least 10 or 20 
  if NDIM % 2 == 1:
      print("Make the grid length even parity to make it truly systolic")
      import pdb; pdb.set_trace()
  
  
  history_xval = 0  # this is used to track 1 specific particle
  
  NScenarios = 2  # we get smoother distributions if se generate several different gases and average the amplitube over them, but feel free to set it to 1
  if not(bDiscrete):
      NScenarios = 2
  
  halfband = 3

  if bDiscrete:
      grid = np.zeros((NScenarios, NDIM, NDIM, 4)).astype("int")
      gridlatticegas = np.zeros((NScenarios, NDIM, NDIM, 4)).astype("int") # this will be used to feed the random numbers
  else:
      grid = np.zeros((NScenarios, NDIM, NDIM, 4))
      gridlatticegas = np.zeros((NScenarios, NDIM, NDIM, 4)) # this will be used to feed the random numbers
  
  
  #import pdb; pdb.set_trace()
  gridlatticegass = np.zeros(grid.shape) # np.random.random_sample(grid.shape) < 0.5
  gridlatticegass = gridlatticegass.astype("int")

  #import pdb; pdb.set_trace()
  ##Lagrance multipliers
  #gridlatticegass_amplambda = np.zeros(grid.shape) # Lagrange lambda for penalizing non-conservation of amplitudes
  #gridlatticegass_xmomlambda = np.zeros(grid.shape) # Lagrange lambda for penalizing non-conservation of x-momentum
  #gridlatticegass_ymomlambda = np.zeros(grid.shape) # Lagrange lambda for penalizing non-conservation of x-momentum
  #gridlatticegass_ymomlambda = np.zeros(grid.shape) # Lagrange lambda for penalizing non-conservation of x-momentum
  
  #disable particle count lambda initially -- 
  #gridlatticegass_nlambda = np.zeros(grid.shape) # Lagrange lambda for penalizing non-conservation of particle number (only the even scenarios are used)
  
  #bFermion = False
  #if bFermion:
  #    gridlatticegass_paulilambda = np.zeros(grid.shape) # Lagrange lambda for penalizing violation of pauli exclusion principle (more than one particle per arc -- only the even scenarios are used)
  

  bEnableSingleParticleTracking = False
  
  
  if bEnableSingleParticleTracking:
    tracker = Tracker(0, NDIM//2, NDIM//2, 0)  

  
	
	
  #########################################
  #########################################
  ##
  ##  INITIAL CONDITIONS (begin)
  ##
  #########################################
  #########################################

  #                    ^ #
  #                 dA | #   
  #                    | #   
  #                    v #   
  #                      ^    
  #                 dB   |   
  #                      v    
  #                      #   
  #                      #   
  #                      #   
  #                      #   
  #                      # <====d2screen====>  
  #                      #   
  #                      #   
  #                      #   
  #                      #   
  #                      ^    
  #                 dB   |   
  #                      v      
  #                    ^ #
  #                    | #   
  #                 dA | #   
  #                    v #   



  d2split = 0
  d2screen = [1,2,3,4,6,8,10,12]
  dA = NDIM//2-4 # make it even parity if you want everything to be systolic
  dB = 1
  pulsetime = 8# 12
  
  #grid = np.random.random_sample(grid.shape) <= 0.0

 
  grid[:, d2split, :dA ,:] = 0 #  
  grid[:, d2split, -dA: ,:] = 0 #  
  grid[:, d2split, dA+dB:-(dA+dB) ,:] = 0 #  

  for i in range(len(d2screen)):
      grid[:, d2split + d2screen[i], : ,:] = 0 #  
  
  #startampL = []
  #startampR = []
  
  if bDiscrete:
      screen = np.zeros((len(d2screen), NDIM,)).astype("int")
      screenalt = np.zeros((len(d2screen), NDIM,)).astype("int")
  
  else:
      screen = np.zeros((len(d2screen), NDIM,))
      screenalt = np.zeros((len(d2screen), NDIM,))
  




  imgplot = plt.imshow(1 - Flatten(grid)).set_cmap('hot')  # show 1st scenario
  plt.axis('off')
  plt.show()
  
  
  #import pdb; pdb.set_trace()
  for iscr in range(screen.shape[0]):
      screenplot = plt.plot(np.arange(0,NDIM), screen[iscr,:])  # show 1st scenario
  plt.show()
  
  
  

  #########################################
  #########################################
  ##
  ##  INITIAL CONDITIONS (end)
  ##
  #########################################
  #########################################

  thisparity = 0
  
  
  TMAX = max([1200000, nplotmod])
 
  minchunk = []
  maxchunk = []
  minmaxrefresh = 1000000 #1000
  minthischunk = 0
  maxthischunk = 0
 
  #i_scenario = 0
  
  print2("sum " + str(np.sum(grid)))
  

            
  t_lastdump = 0
  #screenpty = ((d2split + d2screen + 0) % 2) 
  slotApty = ((dA + 0) % 2) 
  slotBpty = ((NDIM-(dA)) % 2) 
  
  if slotApty != slotBpty:
      import pdb; pdb.set_trace()
      print("ERROR: emission points should have same parity ")
 
  NPer = pulsetime
  refreshtime = 10*pulsetime # wipe the screen each time
  
  
  bNoWipe = True # disables any grid erasures -- test for amplitube preservation
  
  #import pdb; pdb.set_trace()
  
  for t in range(TMAX):
  
  
    if t % refreshtime == 0 and not(bNoWipe):
        #import pdb; pdb.set_trace()
        grid = grid * 0


    if not(bNoWipe):
        #grid[:, d2split, : ,:] = 0 #  
        grid[:, d2split-1, : ,:] = 0 #  
        grid[:, d2split-2, : ,:] = 0 #  
    
    

    #import pdb; pdb.set_trace()
    if True: # t % 10 == 0:  
      if True: # t % pulsetime == 0:          
          if bDiscrete:
              for i in range(grid.shape[0]):
                  #import pdb; pdb.set_trace()
                  
                  if t == 0 or not(bNoWipe):
                      import pdb; pdb.set_trace()
                      #grid[i, d2split, dA ,:] = randamp(grid[i, d2split, dA ,:].shape, np.cos(2 * np.pi * (t % NPer) / float(NPer))) # 1 #  
                      #grid[i, d2split, -dA ,:] = randamp(grid[i, d2split, -dA ,:].shape, -np.cos(2 * np.pi * (t % NPer) / float(NPer)))# -1 #  

                      for iscr in range(grid.shape[0]):
                          # first boundary condition
                          grid[iscr, d2split, dA ,0] = 1 #randamp(grid[:, d2split, dA ,0].shape, np.cos(2 * np.pi * (t % NPer) / float(NPer))) # 1 #  
                          
                          
                          # second opposing boundary condition
                          #grid[iscr, d2split, -dA ,0] = -1 #randamp(grid[:, d2split, -dA ,0].shape, -np.cos(2 * np.pi * (t % NPer) / float(NPer)))# -1 #  
                      
                      import pdb; pdb.set_trace()
                      
                      
                      if iscr == 0:
                          print("here is the sum of amplitudes in LHS of grid (scenario 0): " + str(np.sum(grid[0, 0, 0:8])))
                          print("here is the sum of amplitudes in RHS of grid (scenario 0): " + str(np.sum(grid[0, 0, 9:])))
                  
                  
                  
          else:
              #import pdb; pdb.set_trace()
              grid[:, d2split, dA ,:] = np.cos(2 * np.pi * (t % NPer) / float(NPer)) # 1 #  
              grid[:, d2split, -dA ,:] = -np.cos(2 * np.pi * (t % NPer) / float(NPer))# -1 #  
              
              

    if bDiscrete:
        #propagate_rotate90oncollision_ONESWEEP(gridlatticegass)
        #classicalwave(grid, gridlatticegass)
        #bCplusplus = False
        if bCplusplus:
            NDimX = NDIM
            NDimY = NDIM
            NDir = 4
            randlong = int(np.random.rand() * 1000 * 1000 * 1000)
            #grid = statgrid2dclassical.cythonclassicalwave(grid,  np.zeros(grid.shape).astype("int"), NScenarios, NDimX, NDimY, NDir, randlong)
        else:
          #bUseLatticeGasForRand = False 
          #import pdb; pdb.set_trace()
          if bUseLatticeGasForRand:
              propagate_rotate90oncollision_ONESWEEP(gridlatticegass)
              
          #import pdb; pdb.set_trace()
          #classicalwave(grid, gridlatticegass) # discrete version
          
          #import pdb; pdb.set_trace()
          wavegridupdate(t, grid, #fermpenalty, n2penalty, bFermion = False, 
                              4)
          #import pdb; pdb.set_trace()
            
    else:

        classicalwavefloatingpt(grid)
    
    if bEnableSingleParticleTracking:
      tracker.postupdate(grid)
      
      
    if t < 100 and bDiscrete:
        #import pdb; pdb.set_trace()
        grid[1,:,:,:] = grid[1,:,:,:] * 0
        #import pdb; pdb.set_trace()
        flatten2dgridandcsv(grid, True, t, "flattenedgrid.csv")
      
    #zero_out_parity(grid, t % 2)   
    mingrid = np.min(grid)
    maxgrid = np.max(grid)
    
    if t % minmaxrefresh == 0:
        minchunk.append(copy.copy(minthischunk))
        maxchunk.append(copy.copy(maxthischunk))
        minthischunk = 0
        maxthischunk = 0
    
    if mingrid < minthischunk:
        minthischunk = copy.copy(mingrid)
    if maxgrid > maxthischunk:
        maxthischunk = copy.copy(maxgrid)
        
        
    print2(str(t))
    print2("sum     " + str(np.sum(grid)))
    print2("sumXmom " + str(np.sum(grid[:,:,:,Comp000] - grid[:,:,:,Comp180])))
    #print2("sumYmom " + str(np.sum(grid[:,:,:,Comp090] - grid[:,:,:,Comp270])) +  ' max ' + str(np.max(grid)) + ' min ' + str(np.min(grid)))
    print2("sumYmom " + str(np.sum(grid[:,:,:,Comp090] - grid[:,:,:,Comp270])) +  ' max ' + str(np.max(grid)) + ' min ' + str(np.min(grid)) +  ' 90pct ' + str(np.round(np.percentile(grid, 90))) + ' 10pct ' + str(np.round(np.percentile(grid, 10)))                    )
    

    
    if True: # t % 10 == 0:

      
      bFlatten = True
      
      if bFlatten:
        if t % nplotmod == 0:
            imggrid = 1 - Flatten(grid)
        
        #import pdb; pdb.set_trace()
        
        if bEnableSingleParticleTracking:
            tracker.addhistory2grid(imggrid, history_xval)
        
        if t % nplotmod == 0:  
            imgplot = plt.imshow(imggrid).set_cmap('hot')
      
      else:
        if t % nplotmod == 0:
            imggrid0 = Flatten(grid)
            imggrid = 1 - imggrid0
        
        #import pdb; pdb.set_trace()
        
        if bEnableSingleParticleTracking:
            tracker.addhistory2grid(imggrid, history_xval)

 
        if t % nplotmod == 0:
            imgplot = plt.contour( imggrid0)
      
      
      
      bSkipPlot = bDiscrete
      if (t % nplotmod == 0) and not(bSkipPlot):
          
          plt.axis('off')
          plt.show()
      
          #import pdb; pdb.set_trace()
          for iscr in range(screen.shape[0]):
              screenplot = plt.plot(np.arange(0,NDIM), screen[iscr,:])  # show 1st scenario
          plt.show()

          
          for iscr in range(screenalt.shape[0]):
              screenplot = plt.plot(np.arange(0,NDIM), screenalt[iscr,:])  # show 1st scenario
          screenplot = plt.plot(screenalt)  # show 1st scenario
          #import pdb; pdb.set_trace()
          plt.show()
          

  
          plt.plot(np.arange(screen.shape[0]), np.sum(screen, 1))
          
          plt.plot(np.arange(screenalt.shape[0]), np.sum(screenalt, 1))
          plt.show()
          

      if bDiscrete:
          subscreen = np.zeros(screen.shape).astype("int")
          subscreenalt = np.zeros(screen.shape).astype("int")
      else:
          subscreen = np.zeros(screen.shape)
          subscreenalt = np.zeros(screen.shape)
        


      
      if bDiscrete:
          for iscen in range(0, grid.shape[0], 2):
              for j in range(grid.shape[3]):
                  for iscr in range(subscreen.shape[0]):
                      subscreen[iscr, :] += grid[iscen, d2split + d2screen[iscr], : , j] * grid[min([iscen+1, NScenarios-1]), d2split + d2screen[iscr], : , j]
                      subscreenalt[iscr, :] += grid[iscen, d2split + d2screen[iscr] + 1, : , j] * grid[min([iscen+1, NScenarios-1]), d2split + d2screen[iscr], : , j]
      
      else:
          for iscen in range(0, grid.shape[0]):
              for j in range(grid.shape[3]):
                  for iscr in range(subscreen.shape[0]):
                      subscreen[iscr, :] += grid[iscen, d2split + d2screen[iscr], : , j] * grid[min([iscen, NScenarios-1]), d2split + d2screen[iscr], : , j]
                      subscreenalt[iscr, :] += grid[iscen, d2split + d2screen[iscr] + 1, : , j] * grid[min([iscen, NScenarios-1]), d2split + d2screen[iscr], : , j]
          
      screen += subscreen
      screenalt += subscreenalt
      
      #import pdb;pdb.set_trace()

      if False and t > 0 and (t % nplotmod) == 0: # 5000
          #import pdb; pdb.set_trace()
          subscreenfname = 'c:\\quant\\' + opts.picklefile + '_' + str(opts.seed) + '.pkl'
          subscreenfnamealt = 'c:\\quant\\' + opts.picklefile + '_' + str(opts.seed) + 'alt.pkl'
          if os.path.exists(subscreenfname):
              f = open(subscreenfname, 'rb')
              narr, oldarr, oldminchunk, oldmaxchunk = pickle.load(f)
              
              f.close()
              
              if oldarr.shape != subscreen.shape:
                  print("cannot update screen -- wrong shape.")
                  import pdb; pdb.set_trace()
              narr  += t - t_lastdump
          else:
              oldarr = np.zeros(subscreen.shape).astype("int")
              narr = 0
          oldarr = oldarr + subscreen
          
          f = open(subscreenfname, 'wb')
          pickle.dump((narr, oldarr, minchunk, maxchunk), f)
          f.close()  
          
          
          if os.path.exists(subscreenfnamealt):
              f = open(subscreenfnamealt, 'rb')
              narr, oldarr, oldminchunk, oldmaxchunk = pickle.load(f)
              #import pdb;pdb.set_trace()
              f.close()
              
              if oldarr.shape != subscreenalt.shape:
                  print("cannot update screen -- wrong shape.")
                  import pdb; pdb.set_trace()
              narr  += t - t_lastdump
          else:
              oldarr = np.zeros(subscreenalt.shape).astype("int")
              narr = 0
          oldarr = oldarr + subscreenalt
          
          f = open(subscreenfnamealt, 'wb')
          pickle.dump((narr, oldarr, minchunk, maxchunk), f)
          f.close()  
          #import pdb; pdb.set_trace()        
      #grid[:, 0, NDIM//2, :] = 0
      
      if True: # t % pulsetime == 0:          

          pass
          
      if not(bNoWipe):
          grid[:, -1:, :, :] = 0
          # now absorb along the other axis too
          #grid[:, :,  :1,  :] = 0
          grid[:, :, 0, :] = 0
          grid[:, :, 1, :] = 0
          grid[:, :, -1, :] = 0
          grid[:, :, -2, :] = 0
  
      
      



if __name__ == '__main__':   
  ProcessFile()     






"""


# use this script to look at some of the resultant distributions

import matplotlib as plt
import pickle
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import matplotlib.pyplot as plt
import matplotlib.image as mpimg

f = open('c:\\quant\\contagain_0.pkl', 'rb')
narr, ngrid, nmin, nmax = pickle.load(f)
f.close()


plt.plot(np.arange(ngrid.shape[1]), ngrid.T)
plt.show()
np.sum(ngrid)



plt.plot(np.arange(ngrid.shape[1]), np.sum(ngrid,0)) # totally not legit, but whatever...
plt.show()

b = []
a = np.sum(ngrid,0)
for i in range(1, len(a)-1):
    b.append(np.mean([a[i-1],a[i],a[i+1]]))
plt.plot(b)
plt.show()


plt.plot(np.arange(ngrid.shape[1]), ngrid[-1,:])

"""




"""

# use this script for comparing "continuous" distributions to discrete ones, scaling them so that the overall area under the curves are equal

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import matplotlib.pyplot as plt
import matplotlib.image as mpimg

f = open('c:\\quant\\fixparity_357883476.pkl', 'rb')
narr, ngrid, nmin, nmax = pickle.load(f)
f.close()






f = open('c:\\quant\\contagain_0.pkl', 'rb')
narrcont, ngridcont, nmincont, nmaxcont = pickle.load(f)
f.close()



for i in range(ngrid.shape[0]):
    
    plt.plot(np.arange(ngrid.shape[1]), ngridcont[i,:])
    plt.show()
    



for i in range(ngrid.shape[0]):
    normvol = np.sum(ngridcont[i,:])
    normgrid = ngrid[i,:] * (normvol / np.sum(ngrid[i,:]) )
    y = np.zeros((2, normgrid.shape[0]))
    y[0, :] = normgrid
    y[1, :] = ngridcont[i,:]
    
    plt.plot(np.arange(ngrid.shape[1]), y.T)
    plt.show()
    
 


"""







#     python latticegas_youngsdoubleslitexp_penalty.py --discrete --seed 222 --picklefile trash 
