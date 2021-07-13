#!/cygdrive/c/Python27/python

import sys
import os
import optparse
import copy

import numpy as np

import os.path
import datetime, time

import glob

import pandas as pn    

import matplotlib.pyplot as plt      
from matplotlib.backends.backend_pdf import PdfPages

import matplotlib.pyplot as plt
import matplotlib.image as mpimg








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


  p.add_option("--minnonzero", default=20,
                  action="store", dest='minnonzero', type='int',
                  help="minimum nonzero days -- added this because DBEF only had 2 non zero trades in 2015 ")


  
                
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
    
    
    









def  ProcessFile():
  
  
  


  # Deprecated: it's possible to process entire i,j portion of grid at once, which is much faster
  def propagate(grid,  i, j):
    grid[:, i,j,Comp000] = grid[:, i-1,j,Comp000]
    grid[:, i,j,Comp090] = grid[:, i,j-1,Comp090]
    grid[:, i,j,Comp180] = grid[:, (i+1) % grid.shape[1],j,Comp180]
    grid[:, i,j,Comp270] = grid[:, i,(j+1) % grid.shape[2],Comp270]
  
  
  
  
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
        print("mismatch")
        
        
      grid[:, i,j,Comp000] = new000.astype(int)
      grid[:, i,j,Comp090] = new090.astype(int)
      grid[:, i,j,Comp180] = new180.astype(int)
      grid[:, i,j,Comp270] = new270.astype(int)








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
        print("mismatch")
        
         
      grid[:, :, :, Comp000] = new000.astype(int)
      grid[:, :, :, Comp090] = new090.astype(int)
      grid[:, :, :, Comp180] = new180.astype(int)
      grid[:, :, :, Comp270] = new270.astype(int)











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
      




  #INCOMPLETE
  def classicalwave(grid,  gridlatticegas, N_Per_Arc=1):
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
      new000 = (new000 + (new000 * gridlatticegas[:,:,:,Comp000] - np.trunc(new000))).astype("int")
      new090 = (new090 + (new090 * gridlatticegas[:,:,:,Comp090] - np.trunc(new090))).astype("int")
      new180 = (new180 + (new180 * gridlatticegas[:,:,:,Comp180] - np.trunc(new180))).astype("int")
      new270 = (new270 + (new270 * gridlatticegas[:,:,:,Comp270] - np.trunc(new270))).astype("int")
      #new000 = np.round(new000 + 0.01*(np.random.random_sample(new000.shape())-0.5)).astype("int")  #      (new000 * gridlatticegas[:,:,:,Comp000] - np.trunc(new000))).astype("int")
      #new090 = np.round(new090 + 0.01*(np.random.random_sample(new090.shape())-0.5)).astype("int")  #      (new090 + (new090 * gridlatticegas[:,:,:,Comp090] - np.trunc(new090))).astype("int")
      #new180 = np.round(new180 + 0.01*(np.random.random_sample(new180.shape())-0.5)).astype("int")  #      (new180 + (new180 * gridlatticegas[:,:,:,Comp180] - np.trunc(new180))).astype("int")
      #new270 = np.round(new270 + 0.01*(np.random.random_sample(new270.shape())-0.5)).astype("int")  #      (new270 + (new270 * gridlatticegas[:,:,:,Comp270] - np.trunc(new270))).astype("int")
      
      
      
      
        
      grid[:, :, :, Comp000] = new000  #.astype(int) * randgrid_shuffleprob + xup * (1-randgrid_shuffleprob)
      grid[:, :, :, Comp090] = new090  #.astype(int) * randgrid_shuffleprob + yup * (1-randgrid_shuffleprob)
      grid[:, :, :, Comp180] = new180  #.astype(int) * randgrid_shuffleprob + xdn * (1-randgrid_shuffleprob)
      grid[:, :, :, Comp270] = new270  #.astype(int) * randgrid_shuffleprob + ydn * (1-randgrid_shuffleprob)
      
      
      
      
      
      
      
      










      
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

	
	# declare the 2-d lattice (that has 2*2 outgoing flows at each point.
	
  NDIM = 100 # must be at least 10 or 20
  history_xval = 0  # this is used to track 1 specific particle
  
  NScenarios = 64  # we get smoother distributions if se generate several different gases and average the amplitube over them, but feel free to set it to 1
  halfband = 3


  grid = np.zeros((NScenarios, NDIM, NDIM, 4))
  gridlatticegas = np.zeros((NScenarios, NDIM, NDIM, 4)) # this will be used to feed the random numbers
  
  gridlatticegass = np.random.random_sample(grid[:,:,:,:].shape) < 0.5
  gridlatticegass = gridlatticegass.astype("int")
  #grid2 = np.zeros((2*NDIM, 2*NDIM))
  

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


  grid[:,:,:,:] = np.random.random_sample(grid[:,:,:,:].shape) <= 0.0
  grid[:, (NDIM//2 - halfband):(NDIM//2 + halfband), (NDIM//2 - halfband):(NDIM//2 + halfband),:] = (np.random.random_sample(grid[:, (NDIM//2 - halfband):(NDIM//2 + halfband), (NDIM//2 - halfband):(NDIM//2 + halfband),:].shape) < 0.9) * 512
  
  #grid[:,:,:,:] = np.random.random_sample(grid[:,:,:,:].shape) < 0.0
  #grid[:, (NDIM//2 - halfband):(NDIM//2 + halfband), (NDIM//2 - halfband):(NDIM//2 + halfband),:] = np.random.random_sample(grid[:, (NDIM//2 - halfband):(NDIM//2 + halfband), (NDIM//2 - halfband):(NDIM//2 + halfband),:].shape) < 0.9
  
  
  #import pdb; pdb.set_trace()
  for i_scenario in range(NScenarios):
    zero_out_parity(grid, 1)
  


  
  #grid = np.zeros(grid.shape)
  #grid[0,3,3,0] = 1
	



  imgplot = plt.imshow(1 - Flatten(grid)).set_cmap('hot')  # show 1st scenario
  plt.axis('off')
  plt.show()
  
  
  

  #########################################
  #########################################
  ##
  ##  INITIAL CONDITIONS (end)
  ##
  #########################################
  #########################################

  thisparity = 0
  
  
  TMAX = 500
 
  #i_scenario = 0
  
  print("sum " + str(np.sum(grid)))
  

  
  #for t in range(TMAX):
  #  for i in range(grid.shape[1]):
  #    for j in range(grid.shape[2]):
  #      if (i + j + t) % 2 != thisparity:
  #        propagate_rotate90oncollision_multi(grid, i, j)
            
  
  for t in range(TMAX):
  
    #brownian(grid)
    #import pdb; pdb.set_trace()
    # this will be used as the random number generator in the statistical-wave grid
    
    
    
    #poissonbrownian(grid, 0.000)
    
    
    
    #propagate_rotate90oncollision_ONESWEEP(gridlatticegass)
    #classicalwave(grid, gridlatticegass)


    classicalwavefloatingpt(grid)
    
    if bEnableSingleParticleTracking:
      tracker.postupdate(grid)
      
      
    zero_out_parity(grid, t % 2)   
    
    print(str(t))
    print("sum     " + str(np.sum(grid)))
    print("sumXmom " + str(np.sum(grid[:,:,:,Comp000] - grid[:,:,:,Comp180])))
    print("sumYmom " + str(np.sum(grid[:,:,:,Comp090] - grid[:,:,:,Comp270])))
    
    
    
      #print "sum " + str(np.sum(grid)) + " " + str(t) + " " + str(thisparity) + " min " + str(np.min(grid)) + " max " + str(np.max(grid))
      #imgplot = plt.imshow(1 - Flatten(grid, grid2, i_scenario)).set_cmap('hot') 
    
    
    
    if True: # t % 10 == 0:
      #import pdb; pdb.set_trace()
      
      
      
      bFlatten = True
      
      if bFlatten:
       
        imggrid = 1 - Flatten(grid)
        
        #import pdb; pdb.set_trace()
        
        if bEnableSingleParticleTracking:
            tracker.addhistory2grid(imggrid, history_xval)
          
          
        imgplot = plt.imshow(imggrid).set_cmap('hot')
      
      else:

        imggrid0 = Flatten(grid)
        imggrid = 1 - imggrid0
        
        #import pdb; pdb.set_trace()
        
        if bEnableSingleParticleTracking:
            tracker.addhistory2grid(imggrid, history_xval)
          
        #x, y =   
        imgplot = plt.contour( imggrid0)
      
      
      
      

      plt.axis('off')
      plt.show()
    #import pdb ;pdb.set_trace()
    #if t % 20 == 0:
    #  import pdb; pdb.set_trace()




if __name__ == '__main__':   
  ProcessFile()     



# python latticegas_2Dsystolic.py



# python e:\\quant\\latticegas_2Dsystolic_multi.py



# python /Volumes/NONAME/quant/latticegas_2Dsystolic_multi.py


# python latticegas_2Dsystolic_classicalwave.py



                                                                                                        