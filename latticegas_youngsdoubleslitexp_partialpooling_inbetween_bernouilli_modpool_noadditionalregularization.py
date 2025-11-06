#!/cygdrive/c/Python27/python


print("""


Sep-2024

We introduce another form of reducing sum(abs(grid)); not by de-whorling, but
rather, by spewing out rangomly walking particle-antiparticle pairs, with
amplitudes constructed so that they maximally reduce the abs(grid) at whatever
node this pair is emitted. Eventually (and that time may be as long as infinity)
the random-walk particles will meander and annihilate with each other (or
else with  opposite particles emitted at some other node). In either case, the
resultant closed-circuit or loop in space time serves the same purpose as the
whorls that are used to offset the particle count in the previous iteration.

Note that because the distribution of these random walking particles is Gaussian,
then in the continuum limit, the contribution will cancel everywhere (even if
the amplitudes they offset are such that  what gets emitted at one arc is always
positive, whereas what gets emitted at the other arc is negative)


Result: Seems to work, but it took a while, so maybe there are errors. For now,
The Fermi percentage (i.e. number of arcks which have abs amplitudes less than 1)
is close to 85%, same as with the earlier method of de-whorling by plaquettes.

The adjunct "heat bath" of random walking particles (that represent the endpoints
of the homogenous circuits we keep forcing onto the space in order to make the
absolute number of particles bound) has a density that I first feared might
also grow without bound, but the absolute square of that amplitude likewise 
levels off to a value proportional to the size of the grid (doubling a
2-d space increases the equilibrium value of that sum by 4, etc.)

"""
    
    
    
"""
Sep-2024

This is same as latticegas_youngsdoubleslitexp_partialpooling.py (we turn off bDeWhorl)
except we add in an alpha -- if alpha is 1.0, there is "mod-pooling" and growth is sqrt(T)

If we set alpha to zero, then the growth is ~T (i.e. pure Brownian-Huygens)

For alpha in between those extremes, we take twice the incoming sum of arcs at any node and
then divide into two chunks -- sized alpha and 1-alpha -- and mod-pool the former, while
doing a pure Brownian-Huygets for the latter.

I expect that the growth will eventually always be linear if alpha > 0, but depending on how
small it is, growth will be sart(T) until the linear term overwhelms the growth


Result: hunch is correct, but given the way we set it up, you need to set the alpha close to 1.0 -- like .99 to see a gradual
transition from sqrt(T) to T as T increases; for significantly lower alpha, the growth quickly becomes like T

"""

"""

Jan-2024
In this version, the variance reduction is much like before,
except that instead of creating just a single particle-antiparticle
pair to offset any instance of a higher-than-desired
coincidence of positive and negative particles along to outgoing
arcs, we will instead create a LOOP or WHORL. (Note
that if the two arcs are oppositely aligned, the whorl will
span at least two "squares" and can connect the two
outgoing arcs in two ways.)

RESULT: This actually seems successful. Here are some things
I learned.

If you set the threshold ("compresss threshold") higher,
(say 20) then the fermi fraction (i.e. the portion of arcs
in which the absolute amplitude is 1 or 0) is
low; if you set the compress threshold to 1 then
it's around 85% REGARDLESS OF THE START DENSITY
(assuming it's sparse -- obviously this won't be true
if  the fermi faction starts out tiny).

The fraction doesn't depend on the size of the lattice either, which
is even more suprirising.


It doesn't matter if, for the 0-2 or 1-3
offset (for which the offsetting whirl spans two
adjacent plaquettes) one chooses the left/right or
else top/bottom two plaquettes. The results are the same.



""")
      

print("""

Jan-2024
In this version, the ONLY variance reduction is that the
Brownian particles produced at any node are annihilated
(so that the outgoing Brownian particles are all the same sign).
The method we've used previously not only pools the
outgoing Brownian particles but also smoothly spreads them out so
that if the pooled number is a multiple of N (i.e. number of arcs)
then the output for that node is deterministic. My hunch is that
the growth in absolute numer of particles is sqrt(T) for either
of these cases, since smoothing out the outward flow
along all the available arcs doesn't affect the total
number of particles produced there


RESULT: My hunch was incorrect. The partial pooling is pretty close
to 1 -- about 95% (but maybe as high as 100%).
So the evening out of  the arcs is also vital in reducing the growth.




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
#import pandas as pn    

import matplotlib.pyplot as plt      
from matplotlib.backends.backend_pdf import PdfPages

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import random as rn
from bisect import bisect_left

#import statgrid2dclassical # don't need this cython routine after all




tHisTArr = []
absTotArr = []
l2Arr = []
maxarr = []
ferminoic = []
add2HistoryIncrement = 100
pklFileDumpIncrement = 1000


global CUTOFF
global FAC
global EQWT


CUTOFF = 100
FAC = 0.1 # FAC is only used in varreductionC -- if, after modding and such, a random particle must be generated, the prob
            # of choosing an arc X whose amplitude (prior to putting on the random particle) is of the opposite sign will be
            #
            #  Prob_x  is_proportional_to  1 + X_amp_signFac * (abs(ampX)-cutoff)
            # 
            # where X_amp_signFac is 1 (or 0) depending on whether the amplitude at X is of the opposite (or same) sign
            # the overall normalizaton const is then the Prob_X for all the outgoing arcs

            # by making CUTOFF very large or FAC very small, we can approach linearity
            #
            # we implicitly make the RANDIMIZATION ASSUMPTION -- i.e. the boost in probability averages out to zero over all the
            # scenarios in which probability density would be measured.
CUTOFF = 100
FAC = 0.0

EQWT = [0.25, 0.50, 0.75]



global NHistogram
global HistPartChange    # total number of particles going out minus tot number coming in (ignore no-particle)
global HistSpanChange    # max(inamp) - min(inamp) - max(outamp) - min(outamp) where the max/min refers to all the arcs converging on  a point
        



NHistogram = 20
HistPartChange = np.zeros((NHistogram,))     # total number of particles going out minus tot number coming in (ignore no-particle)
HistSpanChange = np.zeros((NHistogram,))      # max(inamp) - min(inamp) - max(outamp) - min(outamp) where the max/min refers to all the arcs converging on  a point
       











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


  p.add_option("--alpha", default=1.0,
                  action="store", dest='alpha', type='float',
                  help="if 1.0, the process is pure 'modulo/pooled'; if 0, it is pure bernouilli (i.e. 'true Brownian-Huygens'); the alpha multiplies twice times the sum of the incoming flow, etc. ")


  
                
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


def print2grid(grid, t):
    def chstr(xstr, pos, c):
        mylist = list(xstr)
        mylist[pos] = c
        return ''.join(mylist)
    spc = 6
    vrt = 4

    flines = [] 
    #              sum(0,0)                x+                  sum(1)                    x-       
    flines.append('_' * spc + '=' * spc + '_' * spc + '=' * spc + '_' * spc + '=' * spc + '_' * spc + '=' * spc)
    for i in range(vrt):
        flines.append('  ||  ' + ' ' * 3 * spc +            '  ||  '    + ' ' * spc + ' ' * spc + ' ' * spc)
    flines.append('_' * spc + ' ' * spc + ' ' * spc  + ' ' * spc + '_' * spc + ' ' * spc) # y- and y+
    for i in range(vrt):
        flines.append('  ||  ' + ' ' * 3 * spc +            '  ||  '    + ' ' * spc + ' ' * spc + ' ' * spc)
    flines.append('_' * spc + '=' * spc + '_' * spc + '=' * spc + '_' * spc + '=' * spc + '_' * spc + '=' * spc) # sum(0,1) and sum(1,1)
    for i in range(vrt):
        flines.append('  ||  ' + ' ' * 3 * spc  +          '  ||  '    + ' ' * spc + ' ' * spc + ' ' * spc)
    flines.append('_' * spc + ' ' * spc + ' ' * spc  + ' ' * spc + '_' * spc + ' ' * spc) # y- and y+
    for i in range(vrt):
        flines.append('  ||  ' + ' ' * 3 * spc  +          '  ||  '    + ' ' * spc + ' ' * spc + ' ' * spc)

    
    if t % 2 == 0:
        flines[0] = chstr(flines[0], 4*spc-1, '>')
        flines[0] = chstr(flines[0], 5*spc, '<')
        flines[1] = chstr(flines[1], 4*spc+2, '/')
        flines[1] = chstr(flines[1], 4*spc+3, '\\')
        flines[2*vrt+1] = chstr(flines[2*vrt+1], 2, '\\')
        flines[2*vrt+1] = chstr(flines[2*vrt+1], 3, '/')
        flines[2*vrt+3] = chstr(flines[2*vrt+3], 2, '/')
        flines[2*vrt+3] = chstr(flines[2*vrt+3], 3, '\\')
        flines[2*vrt+2] = chstr(flines[2*vrt+2], spc, '<')
        flines[2*vrt+2] = chstr(flines[2*vrt+2], 8*spc-1, '>')
        flines[-1] = chstr(flines[-1], 4*spc+2, '\\')
        flines[-1] = chstr(flines[-1], 4*spc+3, '/')
        a = list(flines[0])
        import pdb; pdb.set_trace()
        b = '{:^6d}'.format(np.sum(grid[0,0,0,:])) 
        a[0:spc] = list( b )
        b = '{:^6d}'.format(grid[0,0,0,0])
        a[2*spc:3*spc] = list( b)
        #b = '{:^6d}'.format(np.sum(grid[0,0,1,:])) 
        #a[4*spc:5*spc] = list( b )
        b = '{:^6d}'.format(grid[0,0,0,2])
        a[6*spc:7*spc] = list(b)
        flines[0] = ''.join(a)

        a = list(flines[vrt+1])
        #import pdb; pdb.set_trace()
        b = '{:^6d}'.format(grid[0,0,0,1]) 
        a[0:spc] = list( b )
        b = '{:^6d}'.format(grid[0,1,1,3]) 
        a[4*spc:5*spc] = list( b)
        flines[vrt+1] = ''.join(a)


        a = list(flines[2*vrt+2])
        #b = '{:^6d}'.format(np.sum(grid[0,0,1,:])) 
        #a[0:spc] = list( b )
        b = '{:^6d}'.format(grid[0,1,1,2])
        a[2*spc:3*spc] = list( b)
        b = '{:^6d}'.format(np.sum(grid[0,1,1,:])) 
        a[4*spc:5*spc] = list( b )
        b = '{:^6d}'.format(grid[0,1,1,0])
        a[6*spc:7*spc] = list(b)
        flines[2*vrt+2] = ''.join(a)

        a = list(flines[3*vrt+3])
        #import pdb; pdb.set_trace()
        b = '{:^6d}'.format(grid[0,0,0,3]) 
        a[0:spc] = list( b )
        b = '{:^6d}'.format(grid[0,1,1,1]) 
        a[4*spc:5*spc] = list( b)
        flines[3*vrt+3] = ''.join(a)

    if t % 2 == 1:
        #import pdb; pdb.set_trace()
        flines[0] = chstr(flines[0], spc, '<')
        flines[0] = chstr(flines[0], -1, '>')
        flines[1] = chstr(flines[1], 2, '/')
        flines[1] = chstr(flines[1], 3, '\\')
        
        flines[2*vrt+1] = chstr(flines[2*vrt+1], 4*spc + 2, '\\')
        flines[2*vrt+1] = chstr(flines[2*vrt+1], 4*spc + 3, '/')
        flines[2*vrt+3] = chstr(flines[2*vrt+3], 4*spc + 2, '/')
        flines[2*vrt+3] = chstr(flines[2*vrt+3], 4*spc + 3, '\\')
        flines[2*vrt+2] = chstr(flines[2*vrt+2], 4*spc-1, '>')
        flines[2*vrt+2] = chstr(flines[2*vrt+2], 5*spc, '<')
        flines[-1] = chstr(flines[-1], 2, '\\')
        flines[-1] = chstr(flines[-1], 3, '/')
        a = list(flines[0])
        #import pdb; pdb.set_trace()
        #b = '{:^6d}'.format(np.sum(grid[0,0,0,:])) 
        #a[0:spc] = list( b )
        b = '{:^6d}'.format(grid[0,1,0,2])
        a[2*spc:3*spc] = list( b)
        b = '{:^6d}'.format(np.sum(grid[0,1,0,:])) 
        a[4*spc:5*spc] = list( b )
        b = '{:^6d}'.format(grid[0,1,0,0])
        a[6*spc:7*spc] = list(b)
        flines[0] = ''.join(a)

        a = list(flines[vrt+1])
        #import pdb; pdb.set_trace()
        b = '{:^6d}'.format(grid[0,0,1,3]) 
        a[0:spc] = list( b )
        b = '{:^6d}'.format(grid[0,1,0,1]) 
        a[4*spc:5*spc] = list( b)
        flines[vrt+1] = ''.join(a)


        a = list(flines[2*vrt+2])
        b = '{:^6d}'.format(np.sum(grid[0,0,1,:])) 
        a[0:spc] = list( b )
        b = '{:^6d}'.format(grid[0,0,1,0])
        a[2*spc:3*spc] = list( b)
        #b = '{:^6d}'.format(np.sum(grid[0,1,1,:])) 
        #a[4*spc:5*spc] = list( b )
        b = '{:^6d}'.format(grid[0,0,1,2])
        a[6*spc:7*spc] = list(b)
        flines[2*vrt+2] = ''.join(a)

        a = list(flines[3*vrt+3])
        #import pdb; pdb.set_trace()
        b = '{:^6d}'.format(grid[0,0,1,1]) 
        a[0:spc] = list( b )
        b = '{:^6d}'.format(grid[0,1,0,3]) 
        a[4*spc:5*spc] = list( b)
        flines[3*vrt+3] = ''.join(a)


    for il in flines:
       print(il)
    import pdb; pdb.set_trace()


def print4grid(grid, t):
    def chstr(xstr, pos, c):
        mylist = list(xstr)
        mylist[pos] = c
        return ''.join(mylist)
    spc = 6
    vrt = 4

    #import pdb; pdb.set_trace()
    flines = [] 
    #              sum(0,0)                x+                  sum(0,1)                    x-       
    flines.append(('_' * spc + '=' * spc + '_' * spc + '=' * spc + '_' * spc + '=' * spc + '_' * spc + '=' * spc)*2)
    for i in range(vrt):
        flines.append(('  ||  ' + ' ' * 3 * spc +            '  ||  '    + ' ' * spc + ' ' * spc + ' ' * spc)*2)
    flines.append(('_' * spc + ' ' * spc + ' ' * spc  + ' ' * spc + '_' * spc + ' ' * spc * 3 )*2) # y- and y+
    for i in range(vrt):
        flines.append(('  ||  ' + ' ' * 3 * spc +            '  ||  '    + ' ' * spc + ' ' * spc + ' ' * spc) *2)
    flines.append(('_' * spc + '=' * spc + '_' * spc + '=' * spc + '_' * spc + '=' * spc + '_' * spc + '=' * spc)*2) # sum(0,1) and sum(1,1)
    for i in range(vrt):
        flines.append(('  ||  ' + ' ' * 3 * spc  +          '  ||  '    + ' ' * spc + ' ' * spc + ' ' * spc)*2)
    flines.append(('_' * spc + ' ' * spc + ' ' * spc  + ' ' * spc + '_' * spc + ' ' * spc * 3 )*2) #  # y- and y+
    for i in range(vrt):
        flines.append(('  ||  ' + ' ' * 3 * spc  +          '  ||  '    + ' ' * spc + ' ' * spc + ' ' * spc)*2)
    
    flines = flines + copy.copy(flines)

    
    
    if t % 2 == 0:
        flines[0] = chstr(flines[0], 4*spc-1, '>')
        flines[0] = chstr(flines[0], 5*spc, '<')
        flines[0] = chstr(flines[0], 12*spc-1, '>')
        flines[0] = chstr(flines[0], 13*spc, '<')



        flines[1] = chstr(flines[1], 4*spc+2, '/')
        flines[1] = chstr(flines[1], 4*spc+3, '\\')        
        flines[1] = chstr(flines[1], 12*spc+2, '/')
        flines[1] = chstr(flines[1], 12*spc+3, '\\')

        flines[2*vrt+1] = chstr(flines[2*vrt+1], 2, '\\')
        flines[2*vrt+1] = chstr(flines[2*vrt+1], 3, '/')
        flines[2*vrt+1] = chstr(flines[2*vrt+1], 8*spc+2, '\\')
        flines[2*vrt+1] = chstr(flines[2*vrt+1], 8*spc+3, '/')



        flines[2*vrt+3] = chstr(flines[2*vrt+3], 2, '/')
        flines[2*vrt+3] = chstr(flines[2*vrt+3], 3, '\\')
        flines[2*vrt+3] = chstr(flines[2*vrt+3], 8*spc+2, '/')
        flines[2*vrt+3] = chstr(flines[2*vrt+3], 8*spc+3, '\\')



        #import pdb; pdb.set_trace()
        flines[2*vrt+2] = chstr(flines[2*vrt+2], spc, '<')
        flines[2*vrt+2] = chstr(flines[2*vrt+2], -1, '>')
        flines[2*vrt+2] = chstr(flines[2*vrt+2], 8*spc-1, '>')
        flines[2*vrt+2] = chstr(flines[2*vrt+2], 9*spc, '<')

        flines[4*vrt+4] = chstr(flines[4*vrt+4], 4*spc-1, '>')
        flines[4*vrt+4] = chstr(flines[4*vrt+4], 5*spc, '<')
        flines[4*vrt+4] = chstr(flines[4*vrt+4], 12*spc-1, '>')
        flines[4*vrt+4] = chstr(flines[4*vrt+4], 13*spc, '<')

        
        flines[6*vrt+6] = chstr(flines[6*vrt+6], spc, '<')
        flines[6*vrt+6] = chstr(flines[6*vrt+6], -1, '>')
        flines[6*vrt+6] = chstr(flines[6*vrt+6], 8*spc-1, '>')
        flines[6*vrt+6] = chstr(flines[6*vrt+6], 8*spc+spc, '<')
      
        flines[6*vrt+5] = chstr(flines[6*vrt+5], 2, '\\')
        flines[6*vrt+5] = chstr(flines[6*vrt+5], 3, '/')
        flines[6*vrt+5] = chstr(flines[6*vrt+5], 8*spc+2, '\\')
        flines[6*vrt+5] = chstr(flines[6*vrt+5], 8*spc+3, '/')

        flines[6*vrt+7] = chstr(flines[6*vrt+7], 2, '/')
        flines[6*vrt+7] = chstr(flines[6*vrt+7], 3, '\\')
        flines[6*vrt+7] = chstr(flines[6*vrt+7], 8*spc+2, '/')
        flines[6*vrt+7] = chstr(flines[6*vrt+7], 8*spc+3, '\\')




        hlfwy = len(flines) // 2
        flines[hlfwy+1] = chstr(flines[hlfwy+1], 4*spc+2, '/')
        flines[hlfwy+1] = chstr(flines[hlfwy+1], 4*spc+3, '\\')        
        flines[hlfwy+1] = chstr(flines[hlfwy+1], 12*spc+2, '/')
        flines[hlfwy+1] = chstr(flines[hlfwy+1], 12*spc+3, '\\')

        flines[hlfwy-1] = chstr(flines[hlfwy-1], 4*spc+2, '\\')
        flines[hlfwy-1] = chstr(flines[hlfwy-1], 4*spc+3, '/')
        flines[hlfwy-1] = chstr(flines[hlfwy-1], 12*spc+2, '\\')
        flines[hlfwy-1] = chstr(flines[hlfwy-1], 12*spc+3, '/')




        flines[-1] = chstr(flines[-1], 4*spc+2, '\\')
        flines[-1] = chstr(flines[-1], 4*spc+3, '/')
        flines[-1] = chstr(flines[-1], 12*spc+2, '\\')
        flines[-1] = chstr(flines[-1], 12*spc+3, '/')

        #import pdb; pdb.set_trace()
        a = list(flines[0])
        #import pdb; pdb.set_trace()
        b = '{:^6d}'.format(np.sum(grid[0,0,0,:])) 
        a[0:spc] = list( b )
        b = '{:^6d}'.format(grid[0,0,0,0])
        a[2*spc:3*spc] = list( b)
        #b = '{:^6d}'.format(np.sum(grid[0,0,1,:])) 
        #a[4*spc:5*spc] = list( b )
        b = '{:^6d}'.format(grid[0,2,0,2])
        a[6*spc:7*spc] = list(b)
        #quadrant2,0
        b = '{:^6d}'.format(np.sum(grid[0,2,0,:])) 
        a[8*spc:9*spc] = list( b )
        b = '{:^6d}'.format(grid[0,2,0,0])
        a[10*spc:11*spc] = list( b)
        #b = '{:^6d}'.format(np.sum(grid[0,0,1,:])) 
        #a[4*spc:5*spc] = list( b )
        b = '{:^6d}'.format(grid[0,0,0,2])
        a[14*spc:15*spc] = list(b)
        flines[0] = ''.join(a)

        a = list(flines[vrt+1])
        #import pdb; pdb.set_trace()
        b = '{:^6d}'.format(grid[0,0,0,1]) 
        a[0:spc] = list( b )
        b = '{:^6d}'.format(grid[0,1,1,3]) 
        a[4*spc:5*spc] = list( b)
        #quadrant2,1
        b = '{:^6d}'.format(grid[0,2,0,1]) 
        a[8*spc:9*spc] = list( b )
        b = '{:^6d}'.format(grid[0,3,1,3]) 
        a[12*spc:13*spc] = list( b)
        flines[vrt+1] = ''.join(a)



        a = list(flines[2*vrt+2])
        #b = '{:^6d}'.format(np.sum(grid[0,0,1,:])) 
        #a[0:spc] = list( b )
        b = '{:^6d}'.format(grid[0,1,1,2])
        a[2*spc:3*spc] = list( b)
        b = '{:^6d}'.format(np.sum(grid[0,1,1,:])) 
        a[4*spc:5*spc] = list( b )
        b = '{:^6d}'.format(grid[0,1,1,0])
        a[6*spc:7*spc] = list(b)
        #quadrant2,2
        #b = '{:^6d}'.format(np.sum(grid[0,0,1,:])) 
        #a[0:spc] = list( b )
        b = '{:^6d}'.format(grid[0,3,1,2])
        a[10*spc:11*spc] = list( b)
        b = '{:^6d}'.format(np.sum(grid[0,3,1,:])) 
        a[12*spc:13*spc] = list( b )
        b = '{:^6d}'.format(grid[0,3,1,0])
        a[14*spc:15*spc] = list(b)
        flines[2*vrt+2] = ''.join(a)

        
        a = list(flines[3*vrt+3])
        b = '{:^6d}'.format(grid[0,0,2,3]) 
        a[0:spc] = list( b )
        b = '{:^6d}'.format(grid[0,1,1,1]) 
        a[4*spc:5*spc] = list( b)
        #quadrant2,3
        b = '{:^6d}'.format(grid[0,2,2,3]) 
        a[8*spc:9*spc] = list( b )
        b = '{:^6d}'.format(grid[0,3,1,1]) 
        a[12*spc:13*spc] = list( b)
        flines[3*vrt+3] = ''.join(a)
        #import pdb; pdb.set_trace()


        # 2nd half
        hlfwy = len(flines) // 2
    
        a = list(flines[hlfwy+0])
        #import pdb; pdb.set_trace()
        b = '{:^6d}'.format(np.sum(grid[0,0,2,:])) 
        a[0:spc] = list( b )
        b = '{:^6d}'.format(grid[0,0,2,0])
        a[2*spc:3*spc] = list( b)
        #b = '{:^6d}'.format(np.sum(grid[0,0,1,:])) 
        #a[4*spc:5*spc] = list( b )
        b = '{:^6d}'.format(grid[0,2,2,2])
        a[6*spc:7*spc] = list(b)
        #quadrant2,0
        b = '{:^6d}'.format(np.sum(grid[0,2,2,:])) 
        a[8*spc:9*spc] = list( b )
        b = '{:^6d}'.format(grid[0,2,2,0])
        a[10*spc:11*spc] = list( b)
        #b = '{:^6d}'.format(np.sum(grid[0,0,1,:])) 
        #a[4*spc:5*spc] = list( b )
        b = '{:^6d}'.format(grid[0,0,2,2])
        a[14*spc:15*spc] = list(b)
        flines[hlfwy+0] = ''.join(a)

        a = list(flines[hlfwy+vrt+1])
        #import pdb; pdb.set_trace()
        b = '{:^6d}'.format(grid[0,0,2,1]) 
        a[0:spc] = list( b )
        b = '{:^6d}'.format(grid[0,1,3,3]) 
        a[4*spc:5*spc] = list( b)
        #quadrant2,1
        b = '{:^6d}'.format(grid[0,2,2,1]) 
        a[8*spc:9*spc] = list( b )
        b = '{:^6d}'.format(grid[0,3,3,3]) 
        a[12*spc:13*spc] = list( b)
        flines[hlfwy+vrt+1] = ''.join(a)



        a = list(flines[hlfwy+2*vrt+2])
        #b = '{:^6d}'.format(np.sum(grid[0,0,1,:])) 
        #a[0:spc] = list( b )
        b = '{:^6d}'.format(grid[0,1,3,2])
        a[2*spc:3*spc] = list( b)
        b = '{:^6d}'.format(np.sum(grid[0,1,3,:])) 
        a[4*spc:5*spc] = list( b )
        b = '{:^6d}'.format(grid[0,1,3,0])
        a[6*spc:7*spc] = list(b)
        #quadrant2,2
        #b = '{:^6d}'.format(np.sum(grid[0,0,1,:])) 
        #a[0:spc] = list( b )
        b = '{:^6d}'.format(grid[0,3,3,2])
        a[10*spc:11*spc] = list( b)
        b = '{:^6d}'.format(np.sum(grid[0,3,3,:])) 
        a[12*spc:13*spc] = list( b )
        b = '{:^6d}'.format(grid[0,3,3,0])
        a[14*spc:15*spc] = list(b)
        flines[hlfwy+2*vrt+2] = ''.join(a)


        a = list(flines[hlfwy+3*vrt+3])
        b = '{:^6d}'.format(grid[0,0,0,3]) 
        a[0:spc] = list( b )
        b = '{:^6d}'.format(grid[0,1,3,1]) 
        a[4*spc:5*spc] = list( b)
        #quadrant2,3
        b = '{:^6d}'.format(grid[0,2,0,3]) 
        a[8*spc:9*spc] = list( b )
        b = '{:^6d}'.format(grid[0,3,3,1]) 
        a[12*spc:13*spc] = list( b)
        flines[hlfwy+3*vrt+3] = ''.join(a)




    if t % 2 == 1:
        #import pdb; pdb.set_trace()

        flines[0] = chstr(flines[0], spc, '<')
        flines[0] = chstr(flines[0], -1, '>')
        flines[0] = chstr(flines[0], 8*spc-1, '>')
        flines[0] = chstr(flines[0], 9*spc, '<')


        flines[1] = chstr(flines[1], 2, '/')
        flines[1] = chstr(flines[1], 3, '\\')
        flines[1] = chstr(flines[1], 8*spc+2, '/')
        flines[1] = chstr(flines[1], 8*spc+3, '\\')
        
        flines[2*vrt+1] = chstr(flines[2*vrt+1], 4*spc + 2, '\\')
        flines[2*vrt+1] = chstr(flines[2*vrt+1], 4*spc + 3, '/')
        flines[2*vrt+1] = chstr(flines[2*vrt+1], 12*spc + 2, '\\')
        flines[2*vrt+1] = chstr(flines[2*vrt+1], 12*spc + 3, '/')

        flines[2*vrt+3] = chstr(flines[2*vrt+3], 4*spc + 2, '/')
        flines[2*vrt+3] = chstr(flines[2*vrt+3], 4*spc + 3, '\\')
        flines[2*vrt+3] = chstr(flines[2*vrt+3], 12*spc + 2, '/')
        flines[2*vrt+3] = chstr(flines[2*vrt+3], 12*spc + 3, '\\')


        flines[2*vrt+2] = chstr(flines[2*vrt+2], 4*spc-1, '>')
        flines[2*vrt+2] = chstr(flines[2*vrt+2], 5*spc, '<')
        flines[2*vrt+2] = chstr(flines[2*vrt+2], 12*spc-1, '>')
        flines[2*vrt+2] = chstr(flines[2*vrt+2], 13*spc, '<')


        #???
        hlfwy = len(flines)//2

        flines[hlfwy+0] = chstr(flines[hlfwy+0], spc, '<')
        flines[hlfwy+0] = chstr(flines[hlfwy+0], -1, '>')
        flines[hlfwy+0] = chstr(flines[hlfwy+0], 8*spc-1, '>')
        flines[hlfwy+0] = chstr(flines[hlfwy+0], 9*spc, '<')


        flines[hlfwy+1] = chstr(flines[hlfwy+1], 2, '/')
        flines[hlfwy+1] = chstr(flines[hlfwy+1], 3, '\\')
        flines[hlfwy+1] = chstr(flines[hlfwy+1], 8*spc+2, '/')
        flines[hlfwy+1] = chstr(flines[hlfwy+1], 8*spc+3, '\\')
        
        flines[hlfwy+2*vrt+1] = chstr(flines[hlfwy+2*vrt+1], 4*spc + 2, '\\')
        flines[hlfwy+2*vrt+1] = chstr(flines[hlfwy+2*vrt+1], 4*spc + 3, '/')
        flines[hlfwy+2*vrt+1] = chstr(flines[hlfwy+2*vrt+1], 12*spc + 2, '\\')
        flines[hlfwy+2*vrt+1] = chstr(flines[hlfwy+2*vrt+1], 12*spc + 3, '/')

        flines[hlfwy+2*vrt+3] = chstr(flines[hlfwy+2*vrt+3], 4*spc + 2, '/')
        flines[hlfwy+2*vrt+3] = chstr(flines[hlfwy+2*vrt+3], 4*spc + 3, '\\')
        flines[hlfwy+2*vrt+3] = chstr(flines[hlfwy+2*vrt+3], 12*spc + 2, '/')
        flines[hlfwy+2*vrt+3] = chstr(flines[hlfwy+2*vrt+3], 12*spc + 3, '\\')


        flines[hlfwy+2*vrt+2] = chstr(flines[hlfwy+2*vrt+2], 4*spc-1, '>')
        flines[hlfwy+2*vrt+2] = chstr(flines[hlfwy+2*vrt+2], 5*spc, '<')
        flines[hlfwy+2*vrt+2] = chstr(flines[hlfwy+2*vrt+2], 12*spc-1, '>')
        flines[hlfwy+2*vrt+2] = chstr(flines[hlfwy+2*vrt+2], 13*spc, '<')



        flines[hlfwy-1] = chstr(flines[hlfwy-1], 2, '\\')
        flines[hlfwy-1] = chstr(flines[hlfwy-1], 3, '/')
        flines[hlfwy-1] = chstr(flines[hlfwy-1], 8*spc+2, '\\')
        flines[hlfwy-1] = chstr(flines[hlfwy-1], 8*spc+3, '/')


        flines[-1] = chstr(flines[-1], 2, '\\')
        flines[-1] = chstr(flines[-1], 3, '/')
        flines[-1] = chstr(flines[-1], 8*spc+2, '\\')
        flines[-1] = chstr(flines[-1], 8*spc+3, '/')




        a = list(flines[0])
        #import pdb; pdb.set_trace()
        #b = '{:^6d}'.format(np.sum(grid[0,0,0,:])) 
        #a[0:spc] = list( b )
        b = '{:^6d}'.format(grid[0,1,0,2])
        a[2*spc:3*spc] = list( b)
        b = '{:^6d}'.format(np.sum(grid[0,1,0,:])) 
        a[4*spc:5*spc] = list( b )
        b = '{:^6d}'.format(grid[0,1,0,0])
        a[6*spc:7*spc] = list(b)

        b = '{:^6d}'.format(grid[0,3,0,2])
        a[10*spc:11*spc] = list( b)
        b = '{:^6d}'.format(np.sum(grid[0,3,0,:])) 
        a[12*spc:13*spc] = list( b )
        b = '{:^6d}'.format(grid[0,3,0,0])


        a[14*spc:15*spc] = list(b)

        flines[0] = ''.join(a)


        #import pdb; pdb.set_trace()
        a = list(flines[vrt+1])
        #import pdb; pdb.set_trace()
        b = '{:^6d}'.format(grid[0,0,1,3]) 
        a[0:spc] = list( b )
        b = '{:^6d}'.format(grid[0,1,0,1]) 
        a[4*spc:5*spc] = list( b)

        #import pdb; pdb.set_trace()
        b = '{:^6d}'.format(grid[0,2,1,3]) 
        a[(8*spc+0):(8*spc+spc)] = list( b )
        b = '{:^6d}'.format(grid[0,3,0,1]) 
        a[(8*spc+4*spc):(8*spc+5*spc)] = list( b)


        flines[vrt+1] = ''.join(a)

        

        a = list(flines[2*vrt+2])
        b = '{:^6d}'.format(np.sum(grid[0,0,1,:])) 
        a[0:spc] = list( b )
        b = '{:^6d}'.format(grid[0,0,1,0])
        a[2*spc:3*spc] = list( b)
        #b = '{:^6d}'.format(np.sum(grid[0,1,1,:])) 
        #a[4*spc:5*spc] = list( b )
        b = '{:^6d}'.format(grid[0,0,1,2])
        a[14*spc:15*spc] = list(b)

        b = '{:^6d}'.format(np.sum(grid[0,2,1,:]))
        a[(8*spc+0):(8*spc+spc)] = list( b )
        b = '{:^6d}'.format(grid[0,2,1,0])
        a[(8*spc + 2*spc):(8*spc + 3*spc)] = list( b)
        #b = '{:^6d}'.format(np.sum(grid[0,1,1,:])) 
        #a[4*spc:5*spc] = list( b )
        b = '{:^6d}'.format(grid[0,2,1,2])
        a[6*spc:7*spc] = list(b)

        flines[2*vrt+2] = ''.join(a)



        a = list(flines[3*vrt+3])
        #import pdb; pdb.set_trace()
        b = '{:^6d}'.format(grid[0,0,1,1]) 
        a[0:spc] = list( b )
        b = '{:^6d}'.format(grid[0,1,2,3]) 
        a[4*spc:5*spc] = list( b)

        b = '{:^6d}'.format(grid[0,2,1,1]) 
        a[(8*spc + 0*spc):(8*spc + 1*spc)] = list( b )
        b = '{:^6d}'.format(grid[0,3,2,3]) 
        a[(8*spc + 4*spc):(8*spc + 5*spc)] = list( b)


        flines[3*vrt+3] = ''.join(a)




        # lower half


        a = list(flines[hlfwy+0])
        #import pdb; pdb.set_trace()
        #b = '{:^6d}'.format(np.sum(grid[0,0,0,:])) 
        #a[0:spc] = list( b )
        b = '{:^6d}'.format(grid[0,1,2,2])
        a[2*spc:3*spc] = list( b)
        b = '{:^6d}'.format(np.sum(grid[0,1,2,:])) 
        a[4*spc:5*spc] = list( b )
        b = '{:^6d}'.format(grid[0,1,2,0])
        a[6*spc:7*spc] = list(b)

        b = '{:^6d}'.format(grid[0,3,2,2])
        a[10*spc:11*spc] = list( b)
        b = '{:^6d}'.format(np.sum(grid[0,3,2,:])) 
        a[12*spc:13*spc] = list( b )
        b = '{:^6d}'.format(grid[0,3,2,0])


        a[14*spc:15*spc] = list(b)

        flines[hlfwy+0] = ''.join(a)


        #import pdb; pdb.set_trace()
        a = list(flines[hlfwy+vrt+1])
        #import pdb; pdb.set_trace()
        b = '{:^6d}'.format(grid[0,0,3,3]) 
        a[0:spc] = list( b )
        b = '{:^6d}'.format(grid[0,1,2,1]) 
        a[4*spc:5*spc] = list( b)

        #import pdb; pdb.set_trace()
        b = '{:^6d}'.format(grid[0,2,3,3]) 
        a[(8*spc+0):(8*spc+spc)] = list( b )
        b = '{:^6d}'.format(grid[0,3,2,1]) 
        a[(8*spc+4*spc):(8*spc+5*spc)] = list( b)


        flines[hlfwy+vrt+1] = ''.join(a)

        

        a = list(flines[hlfwy+2*vrt+2])
        b = '{:^6d}'.format(np.sum(grid[0,0,3,:])) 
        a[0:spc] = list( b )
        b = '{:^6d}'.format(grid[0,0,3,0])
        a[2*spc:3*spc] = list( b)
        #b = '{:^6d}'.format(np.sum(grid[0,1,1,:])) 
        #a[4*spc:5*spc] = list( b )
        b = '{:^6d}'.format(grid[0,2,3,2])
        a[6*spc:7*spc] = list(b)

        b = '{:^6d}'.format(np.sum(grid[0,2,3,:]))
        a[(8*spc+0):(8*spc+spc)] = list( b )
        b = '{:^6d}'.format(grid[0,2,3,0])
        a[(8*spc + 2*spc):(8*spc + 3*spc)] = list( b)
        #b = '{:^6d}'.format(np.sum(grid[0,1,1,:])) 
        #a[4*spc:5*spc] = list( b )
        b = '{:^6d}'.format(grid[0,0,3,2])
        a[14*spc:15*spc] = list(b)

        flines[hlfwy+2*vrt+2] = ''.join(a)



        a = list(flines[hlfwy+3*vrt+3])
        #import pdb; pdb.set_trace()
        b = '{:^6d}'.format(grid[0,0,3,1]) 
        a[0:spc] = list( b )
        b = '{:^6d}'.format(grid[0,1,0,3]) 
        a[4*spc:5*spc] = list( b)

        b = '{:^6d}'.format(grid[0,2,3,1]) 
        a[(8*spc + 0*spc):(8*spc + 1*spc)] = list( b )
        b = '{:^6d}'.format(grid[0,3,0,3]) 
        a[(8*spc + 4*spc):(8*spc + 5*spc)] = list( b)


        flines[hlfwy+3*vrt+3] = ''.join(a)



    for il in flines:
       print(il)
    #import pdb; pdb.set_trace()
      
def GetRandomLocation(x, val):
    ans = []
    for i, ix in enumerate(x):
      if ix == val:
         ans.append(i)
    return rn.choice(ans)


   # unlike the index function, this gets a random

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
          sgnivec = np.sign(ivec)
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


def printgrid(grid, iscen=0):
    print(grid[0,:,:,0])
    print(grid[0,:,:,1])    
    print(grid[0,:,:,2])
    print(grid[0,:,:,3])

# slower than classicalwave() but it preserves amplitube at each node   
def wavegridupdate(t, grid, grid_br): # fermpenalty, n2penalty, bFermion = False, Narc=4):
    global alpha
    """
    assume that even parity cells will be updated
    """
    global abscount
    global reduction_count
    global inc_count
    global reduction_amt
    global inc_amt
    global CUTOFF
    global FAC
    global EQWT
    global modpoolvec

    global NHistogram
    global HistPartChange    # total number of particles going out minus tot number coming in (ignore no-particle)
    global HistSpanChange    # max(inamp) - min(inamp) - max(outamp) - min(outamp) where the max/min refers to all the arcs converging on  a point
           

    N_2zero = 0
    N_all = 0
    N_allzero = 0
    #import pdb; pdb.set_trace()

    # in this version, iscen tells you which edge in the multigraph is being updated

    # remember, grid = np.zeros((NScenarios, NDIM, NDIM, 4)).astype("int"), so that means the number of
    
    def GetLookup(Xamp, sgn):
        if int(np.sign(Xamp)) == int(sgn):
           return 1.0
        absXamp = np.abs(Xamp)
        if absXamp <= CUTOFF:
            return 1.0
        return 1.0 + (absXamp - CUTOFF) * FAC
    
    def GetAllInd(mylist, x):
        myretval = []
        for i, il in enumerate(mylist):
            if x == il:
                myretval.append(i)
        return myretval

    def Shuffle_NoUnmoved(N): # shuffle the indices 0...N-1 such that no index is left unshifted
        randomized_list = list(range(N)) 
        while True:
            rn.shuffle(randomized_list)
            for ii, i in enumerate(randomized_list):
                if ii == i:
                    break
            else:
                return randomized_list

    def incoming_indices(NDIM, ix, iy, i):
        if i == 0:
            return (ix-1, iy)
        elif i == 1:
            return (ix, iy-1)
        elif i == 2:
            return ((ix+1)%NDIM, iy)
        elif i == 3:
            return (ix, (iy+1)%NDIM)

    
    def NetOut(inlist):
        N = len(inlist)
        posportion = 0
        negportion = 0
        poslist = []
        neglist = []
        for ii, iflow in enumerate(inlist):
            if iflow > 0:
                poslist.append(iflow)
                neglist.append(0)
                posportion += iflow
            elif iflow < 0:
                neglist.append(iflow)
                poslist.append(0)
                negportion += iflow
            else:
                poslist.append(0)
                neglist.append(0)

        
        if posportion == -negportion:
            return [0] * N
        elif posportion > -negportion:
            fraclist = np.array(poslist)/posportion
            try:
                shave = np.random.multinomial(-negportion, fraclist)    
            except:
                print("pos")        
            return [poslist[i] - shave[i] for i in range(N)]
        elif -negportion > posportion:
            fraclist = -np.array(neglist)/-negportion # all non-netative
            try:
                shave = np.random.multinomial(posportion, fraclist)
            except:
                import pdb; pdb.set_trace()
            return [neglist[i] + shave[i] for i in range(N)]



    N_pointswithzerobefore = 0
    N_pointswithzeroafter = 0


    bAnnihilateAcrossScenarios = False # changing this to true doesn't change the qrt(T) grown ot the absolute particle count
    #if bAnnihilateAcrossScenarios:
    #   modpoolvec = np.zeros(grid.shape[:-1])
















    bBrownianRegularization = False








    #if t == 2:
    #    import pdb; pdb.set_trace()

    #import pdb; pdb.set_trace()         
    N_arc = grid.shape[3]
    NDIM = grid.shape[1]

    for iscen in range(grid.shape[0]):
        for ix in range(grid.shape[1]):
            for iy in range(grid.shape[2]):
                if (ix + iy + t) % 2 != 1:
                    continue # wrong parity  
                
                #if (ix == 1 and iy == 1 and t == 1):
                #    print("for -- should add 2")
                #    import pdb; pdb.set_trace()

                #if t==1 and np.sum(np.abs([grid_br[iscen, ix-1, iy, 0], grid_br[iscen, ix, iy-1, 1], 
                #     grid_br[iscen, (ix+1) % grid.shape[1], iy, 2], grid_br[iscen, ix, (iy+1) % grid.shape[2], 3]])) != 0:
                #     import pdb; pdb.set_trace()

                
                
                #x = np.array([grid_br[iscen, ix-1, iy, 0], grid_br[iscen, ix, iy-1, 1], grid_br[iscen, (ix+1) % grid.shape[1], iy, 2], grid_br[iscen, ix, (iy+1) % grid.shape[2], 3]])

                if bBrownianRegularization: # and np.sum(np.abs(x)) != 0:

                    x = np.array([grid_br[iscen, ix-1, iy, 0], grid_br[iscen, ix, iy-1, 1], grid_br[iscen, (ix+1) % grid.shape[1], iy, 2], grid_br[iscen, ix, (iy+1) % grid.shape[2], 3]])


                    #if np.sum(np.abs(x)) > 0:
                    #    import pdb; pdb.set_trace()
                    # first we net out whatever can be netted, so that all the incoming flows on the brownian grid are of one sign (or zero)
                    shaved = NetOut([grid_br[iscen, ix-1, iy, 0], grid_br[iscen, ix, iy-1, 1], grid_br[iscen, (ix+1) % grid.shape[1], iy, 2], grid_br[iscen, ix, (iy+1) % grid.shape[2], 3]])
                    
                    # add it to the regular grid
                    grid[iscen, ix-1, iy, 0] += shaved[0]
                    grid[iscen, ix, iy-1, 1] += shaved[1]
                    grid[iscen, (ix+1) % grid.shape[1], iy, 2] += shaved[2]
                    grid[iscen, ix, (iy+1) % grid.shape[2], 3] += shaved[3]

                    grid_br[iscen, ix-1, iy, 0] = 0
                    grid_br[iscen, ix, iy-1, 1] = 0
                    grid_br[iscen, (ix+1) % grid.shape[1], iy, 2] = 0
                    grid_br[iscen, ix, (iy+1) % grid.shape[2], 3] = 0



                    #incoming_along_xpos_br = shaved[0]  # note the origin of this flow is in the xneg direction, i.e. ix-1
                    #incoming_along_ypos_br = shaved[1]
                    #incoming_along_xneg_br = shaved[2]
                    #incoming_along_yneg_br = shaved[3]

                    #incoming_br_list = [incoming_along_xpos_br, incoming_along_ypos_br, incoming_along_xneg_br, incoming_along_yneg_br]

                    incoming_br_list = copy.copy(shaved)
                    incoming_br_sum = np.sum(shaved)
                    incoming_br_sgn = np.sign(incoming_br_sum)
                    incoming_br_abs = np.abs(incoming_br_sum)
                    incoming_br_div = incoming_br_abs // N_arc

                    outgoing_br_mod = incoming_br_abs % N_arc # note we take this on the absolute value of the sum, so we'll multiply it by a sign factor
                    outgoing_br_rand = np.random.multinomial(outgoing_br_mod, [1/float(N_arc)]*N_arc)
                    outgoing_br = incoming_br_div  + np.array(outgoing_br_rand) # up to an overall sign factor, this is the "flat" portion plus the randomized (modded) portion
                    outgoing_br = -incoming_br_sgn * outgoing_br # the aforementioned sign factor

                    grid_br[iscen, ix, iy, :] = grid_br[iscen, ix, iy, :] + outgoing_br
                    # add incoming flow (but negated) from some other (non-identical) arc
                    # this is the effective next link of the random walk

                    #OppFlows = Shuffle_NoUnmoved(NDIM) # early attempt in which we tried to prevent random walks with reversals; upon reflection, that was seen to be unnecessary
                    #for ii, i in enumerate(OppFlows):
                    #    ix_other, iy_other = incoming_indices(NDIM, ix, iy, i)
                     #   grid_br[iscen, ix_other, iy_other, i] -= incoming_br_list[ii]

                    
                inspan = max(grid[iscen,ix,iy,:]) - min(grid[iscen,ix,iy,:])
                
                incoming_along_xpos = copy.copy(grid[iscen, ix-1, iy, 0])  # note the origin of this flow is in the xneg direction, i.e. ix-1
                incoming_along_ypos = copy.copy(grid[iscen, ix, iy-1, 1])
                incoming_along_xneg = copy.copy(grid[iscen, (ix+1) % grid.shape[1], iy, 2])
                incoming_along_yneg = copy.copy(grid[iscen, ix, (iy+1) % grid.shape[2], 3])


                    
                    



                if np.max(np.abs([incoming_along_xpos, incoming_along_ypos, incoming_along_xneg, incoming_along_yneg])) == 0:
                    N_allzero += 1
                    N_all += 1
                    if not(bBrownianRegularization):
                        continue



                """ we'll do this below, using pythonic looping"""
                #grid[iscen, ix-1, iy, 0] = 0 # zero out the prior step's amplitudes
                #grid[iscen, ix, iy-1, 1] = 0
                #grid[iscen, (ix+1) % grid.shape[1], iy, 2] = 0
                #grid[iscen, ix, (iy+1) % grid.shape[2], 3] = 0

                # huygens backprop contribution
                outgoing_along_xpos = -incoming_along_xneg
                outgoing_along_ypos = -incoming_along_yneg
                outgoing_along_xneg = -incoming_along_xpos
                outgoing_along_yneg = -incoming_along_ypos

                brownian_pool = 2 * (incoming_along_xpos + incoming_along_ypos + incoming_along_xneg + incoming_along_yneg)

                

                # (incoming_along_xpos, incoming_along_ypos, incoming_along_xneg, incoming_along_yneg)



                if brownian_pool != 0:

                    #modpool = abspool % N_arc

                    #import pdb; pdb.set_trace()
                    

                    if alpha != 0.0 and alpha != 1.0: # the general "mixed case", handles every case of alpha,
                                                       # but we'll go ahead and do the pure B-H and pure "modpooling" separately

                        pospool = np.max([incoming_along_xneg,0]) + np.max([incoming_along_yneg,0]) + np.max([incoming_along_xpos,0]) + np.max([incoming_along_ypos,0])
                        negpool = np.min([incoming_along_xneg,0]) + np.min([incoming_along_yneg,0]) + np.min([incoming_along_xpos,0]) + np.min([incoming_along_ypos,0])

                        pospool *= 2
                        negpool *= 2

                        pospool = int(np.floor((1.0-alpha) * pospool))
                        negpool = int(np.floor((1.0-alpha) * negpool))

                        brownian_pool -= pospool + negpool

                        abspool = np.abs(brownian_pool)
                        sgnpool = np.sign(brownian_pool)
                        modpool = abspool % N_arc
                        divpool = (abspool // N_arc) 

                        # ad any modulo remainder bits to either the positive-sign B-H output, or the negative-sign B-H output
                        if sgnpool > 0:
                            pospool += modpool
                        elif sgnpool < 0:
                            negpool -= modpool

                        posdeviate = np.random.multinomial(pospool, [1/float(N_arc)]*N_arc)
                        negdeviate = np.random.multinomial(-negpool, [1/float(N_arc)]*N_arc) 

                        outdeviatemod = posdeviate - negdeviate


                        divpool *= sgnpool
                        outgoing_along_xpos += outdeviatemod[0] + divpool
                        outgoing_along_ypos += outdeviatemod[1] + divpool
                        outgoing_along_xneg += outdeviatemod[2] + divpool
                        outgoing_along_yneg += outdeviatemod[3] + divpool
                    
                    elif alpha == 0: # the "pure" Brownian-Huygens case, for which sum(abs(grid)) grows like t

                        # note that for this case we did not need the above calculation of brownianpool, since that is unused here

                        pospool = np.max([incoming_along_xneg,0]) + np.max([incoming_along_yneg,0]) + np.max([incoming_along_xpos,0]) + np.max([incoming_along_ypos,0])
                        negpool = np.min([incoming_along_xneg,0]) + np.min([incoming_along_yneg,0]) + np.min([incoming_along_xpos,0]) + np.min([incoming_along_ypos,0])

                        pospool *= 2
                        negpool *= 2

                        #pospool = int(np.floor((1.0-alpha) * pospool))
                        #negpool = int(np.floor((1.0-alpha) * negpool))

                        #brownian_pool -= pospool + negpool

                        #abspool = np.abs(brownian_pool)
                        #sgnpool = np.sign(brownian_pool)
                        #modpool = abspool % N_arc
                        #divpool = (abspool // N_arc) 

                        # ad any modulo remainder bits to either the positive-sign B-H output, or the negative-sign B-H output
                        #if sgnpool > 0:
                        #    pospool += modpool
                        #elif sgnpool < 0:
                        #    negpool -= modpool

                        posdeviate = np.random.multinomial(pospool, [1/float(N_arc)]*N_arc)
                        negdeviate = np.random.multinomial(-negpool, [1/float(N_arc)]*N_arc) 

                        outdeviatemod = posdeviate - negdeviate


                        #divpool *= sgnpool
                        outgoing_along_xpos += outdeviatemod[0] # + divpool
                        outgoing_along_ypos += outdeviatemod[1] # + divpool
                        outgoing_along_xneg += outdeviatemod[2] # + divpool
                        outgoing_along_yneg += outdeviatemod[3] # + divpool
                    
                    else: # alpha == 1 (the "modpool" case for which sum(abs(grid)) grows like sqrt(t) )

                        # in this case, we don't need a separate pospool and negpool, we just mod out the brownian (i.e. sum of neg and pos) pool
                        #pospool = np.max([incoming_along_xneg,0]) + np.max([incoming_along_yneg,0]) + np.max([incoming_along_xpos,0]) + np.max([incoming_along_ypos,0])
                        #negpool = np.min([incoming_along_xneg,0]) + np.min([incoming_along_yneg,0]) + np.min([incoming_along_xpos,0]) + np.min([incoming_along_ypos,0])

                        #pospool *= 2
                        #negpool *= 2

                        #pospool = int(np.floor((1.0-alpha) * pospool))
                        #negpool = int(np.floor((1.0-alpha) * negpool))

                        #brownian_pool -= pospool + negpool

                        abspool = np.abs(brownian_pool)
                        sgnpool = np.sign(brownian_pool)
                        modpool = abspool % N_arc
                        divpool = (abspool // N_arc) # this (times sgnpool) is a scalar value that gets added to every output arc

                        # this is the randomized "leftover" part
                        outdeviatemod = sgnpool * np.random.multinomial(modpool, [1/float(N_arc)]*N_arc)

                        divpool *= sgnpool
                        outgoing_along_xpos += outdeviatemod[0] + divpool
                        outgoing_along_ypos += outdeviatemod[1] + divpool
                        outgoing_along_xneg += outdeviatemod[2] + divpool
                        outgoing_along_yneg += outdeviatemod[3] + divpool


                inparticlecount = abs(incoming_along_xpos) + abs(incoming_along_ypos) + abs(incoming_along_xneg) + abs(incoming_along_yneg)

                


                grid[iscen, ix, iy, 0] = copy.copy(outgoing_along_xpos)
                grid[iscen, ix, iy, 1] = copy.copy(outgoing_along_ypos)
                grid[iscen, ix, iy, 2] = copy.copy(outgoing_along_xneg)
                grid[iscen, ix, iy, 3] = copy.copy(outgoing_along_yneg)
                outspan = max(grid[iscen,ix,iy,:]) - min(grid[iscen,ix,iy,:])



                outparticlecount = abs(outgoing_along_xpos) + abs(outgoing_along_ypos) + abs(outgoing_along_xneg) + abs(outgoing_along_yneg)

                croppedtotparticlecount = np.min([np.max([outparticlecount - inparticlecount, -NHistogram//2]), NHistogram//2-1]) 
                croppedspancount = np.min([np.max([outspan - inspan, -NHistogram//2]), NHistogram//2-1]) 
                HistPartChange[croppedtotparticlecount//2 + NHistogram//2] += 1
                HistSpanChange[croppedspancount//2 + NHistogram//2] += 1








    bDewhorl = False # i.e. remove "whorls"  of the form
    #              a
    #          --------->
    #         |          ^
    #      b  |          | d
    #         v          |
    #           <--------
    #              c
    #  where a*b < 0 and c*d < 0 and b*c < 0 (implying a**d < 0)
    

    if bDewhorl:


        starttot = 1

        endtot = 0

        #while starttot > endtot:
        if True:
            
            starttot = np.sum(np.abs(grid))
            for iscen in range(grid.shape[0]): 
                for ix in range(grid.shape[1]):
                    for iy in range(grid.shape[2]):
                        #if (ix + iy + t) % 2 != 0:
                        #    continue
                        #if (iscen,ix,iy,t) == (0,3,0,46):
                        #    import pdb; pdb.set_trace()
                        ixp1 = (ix + 1) % grid.shape[1]
                        iyp1 = (iy + 1) % grid.shape[2]
                        a = grid[iscen, ix, iy, 0]
                        b = grid[iscen, ix, iy, 1]
                        c = grid[iscen, ixp1, iyp1, 2]
                        d = grid[iscen, ixp1, iyp1, 3]


                            

                        totpartcountbefore = np.sum(np.abs(np.array([a, b, c, d])))
                        

                        
                        absa = np.abs(a)
                        absb = np.abs(b)
                        absc = np.abs(c)
                        absd = np.abs(d)


                        endresa = (np.abs(a - a) + np.abs(b + a) + np.abs(c - a) + np.abs(d + a))
                        endresb = (np.abs(a - b) + np.abs(b + b) + np.abs(c - b) + np.abs(d + b))
                        endresc = (np.abs(a - c) + np.abs(b + c) + np.abs(c - c) + np.abs(d + c))
                        endresd = (np.abs(a - d) + np.abs(b + d) + np.abs(c - d) + np.abs(d + d))

                        endresA = (np.abs(a - -a) + np.abs(b + -a) + np.abs(c - -a) + np.abs(d + -a))
                        endresB = (np.abs(a - -b) + np.abs(b + -b) + np.abs(c - -b) + np.abs(d + -b))
                        endresC = (np.abs(a - -c) + np.abs(b + -c) + np.abs(c - -c) + np.abs(d + -c))
                        endresD = (np.abs(a - -d) + np.abs(b + -d) + np.abs(c - -d) + np.abs(d + -d))



                        
                        ind = [endresa, endresb, endresc, endresd, endresA, endresB, endresC, endresD].index(np.min([endresa, endresb, endresc, endresd, endresA, endresB, endresC, endresD]))


                        #if ind > 3:
                        #    import pdb; pdb.set_trace()
                        if ind == 0: # a is minimum
                            aa = copy.copy(a)
                            if (absa + absb + absc + absd) <= (np.abs(a - aa) + np.abs(b + aa) + np.abs(c - aa) + np.abs(d + aa)):
                                #print("wouldn't change: ", a, b, c, d, " start ", (absa + absb + absc + absd), "end", (np.abs(a - aa) + np.abs(b + aa) + np.abs(c - aa) + np.abs(d + aa)))
                                #import pdb; pdb.set_trace()
                                continue


                            grid[iscen, ix, iy, 0] -= aa
                            grid[iscen, ix, iy, 1] += aa
                            grid[iscen, ixp1, iyp1, 2] -= aa
                            grid[iscen, ixp1, iyp1, 3] += aa
                        elif ind == 1:
                            bb = copy.copy(b)
                            if (absa + absb + absc + absd) <= (np.abs(a - bb) + np.abs(b + bb) + np.abs(c - bb) + np.abs(d + bb)):
                                #print("wouldn't change: ", a, b, c, d, " start ", (absa + absb + absc + absd), "end", (np.abs(a - bb) + np.abs(b + bb) + np.abs(c - bb) + np.abs(d + bb)))
                                #import pdb; pdb.set_trace()
                                continue
                            
                            grid[iscen, ix, iy, 0] -= bb
                            grid[iscen, ix, iy, 1] += bb
                            grid[iscen, ixp1, iyp1, 2] -= bb
                            grid[iscen, ixp1, iyp1, 3] += bb
                        elif ind == 2:
                            cc = copy.copy(c)
                            if (absa + absb + absc + absd) <= (np.abs(a - cc) + np.abs(b + cc) + np.abs(c - cc) + np.abs(d + cc)):
                                #print("wouldn't change: ", a, b, c, d, " start ", (absa + absb + absc + absd), "end", (np.abs(a - cc) + np.abs(b + cc) + np.abs(c - cc) + np.abs(d + cc)))
                                #import pdb; pdb.set_trace()
                                continue
                                                    
                            grid[iscen, ix, iy, 0] -= cc
                            grid[iscen, ix, iy, 1] += cc
                            grid[iscen, ixp1, iyp1, 2] -= cc
                            grid[iscen, ixp1, iyp1, 3] += cc  
                        elif ind == 3:               
                            dd = copy.copy(d)
                            if (absa + absb + absc + absd) <= (np.abs(a - dd) + np.abs(b + dd) + np.abs(c - dd) + np.abs(d + dd)):
                                #print("wouldn't change: ", a, b, c, d, " start ", (absa + absb + absc + absd), "end", (np.abs(a - dd) + np.abs(b + dd) + np.abs(c - dd) + np.abs(d + dd)))
                                #import pdb; pdb.set_trace()
                                continue
                                                                
                            grid[iscen, ix, iy, 0] -= dd
                            grid[iscen, ix, iy, 1] += dd
                            grid[iscen, ixp1, iyp1, 2] -= dd
                            grid[iscen, ixp1, iyp1, 3] += dd

                        if ind == 4: # a is minimum
                            aa = -copy.copy(a)
                            if (absa + absb + absc + absd) <= (np.abs(a - aa) + np.abs(b + aa) + np.abs(c - aa) + np.abs(d + aa)):
                                #print("wouldn't change: ", a, b, c, d, " start ", (absa + absb + absc + absd), "end", (np.abs(a - aa) + np.abs(b + aa) + np.abs(c - aa) + np.abs(d + aa)))
                                #import pdb; pdb.set_trace()
                                continue


                            grid[iscen, ix, iy, 0] -= aa
                            grid[iscen, ix, iy, 1] += aa
                            grid[iscen, ixp1, iyp1, 2] -= aa
                            grid[iscen, ixp1, iyp1, 3] += aa
                        elif ind == 5:
                            bb = -copy.copy(b)
                            if (absa + absb + absc + absd) <= (np.abs(a - bb) + np.abs(b + bb) + np.abs(c - bb) + np.abs(d + bb)):
                                #print("wouldn't change: ", a, b, c, d, " start ", (absa + absb + absc + absd), "end", (np.abs(a - bb) + np.abs(b + bb) + np.abs(c - bb) + np.abs(d + bb)))
                                #import pdb; pdb.set_trace()
                                continue
                            
                            grid[iscen, ix, iy, 0] -= bb
                            grid[iscen, ix, iy, 1] += bb
                            grid[iscen, ixp1, iyp1, 2] -= bb
                            grid[iscen, ixp1, iyp1, 3] += bb
                        elif ind == 6:
                            cc = -copy.copy(c)
                            if (absa + absb + absc + absd) <= (np.abs(a - cc) + np.abs(b + cc) + np.abs(c - cc) + np.abs(d + cc)):
                                #print("wouldn't change: ", a, b, c, d, " start ", (absa + absb + absc + absd), "end", (np.abs(a - cc) + np.abs(b + cc) + np.abs(c - cc) + np.abs(d + cc)))
                                #import pdb; pdb.set_trace()
                                continue
                                                    
                            grid[iscen, ix, iy, 0] -= cc
                            grid[iscen, ix, iy, 1] += cc
                            grid[iscen, ixp1, iyp1, 2] -= cc
                            grid[iscen, ixp1, iyp1, 3] += cc  
                        elif ind == 7:               
                            dd = -copy.copy(d)
                            if (absa + absb + absc + absd) <= (np.abs(a - dd) + np.abs(b + dd) + np.abs(c - dd) + np.abs(d + dd)):
                                #print("wouldn't change: ", a, b, c, d, " start ", (absa + absb + absc + absd), "end", (np.abs(a - dd) + np.abs(b + dd) + np.abs(c - dd) + np.abs(d + dd)))
                                #import pdb; pdb.set_trace()
                                continue
                                                                
                            grid[iscen, ix, iy, 0] -= dd
                            grid[iscen, ix, iy, 1] += dd
                            grid[iscen, ixp1, iyp1, 2] -= dd
                            grid[iscen, ixp1, iyp1, 3] += dd
                        


                        #import pdb; pdb.set_trace() # (a,b,c,d,t,iscen,ix,iy) # grid[iscen, ix, iy, 0], grid[iscen, ix, iy, 1], grid[iscen, ixp1, iyp1, 2], grid[iscen, ixp1, iyp1, 3]


                        totpartcountafter = np.sum(np.abs(np.array([grid[iscen, ix, iy, 0], grid[iscen, ix, iy, 1], 
                                        grid[iscen, ixp1, iyp1, 2], grid[iscen, ixp1, iyp1, 3]])))


                        if totpartcountafter >= totpartcountbefore:
                            #print("didn't work")
                            import pdb; pdb.set_trace()
                        else:
                            pass
                            #print("went down from ", str(totpartcountbefore), " to ", totpartcountafter, starttot,np.sum(np.abs(grid)))
            endtot = np.sum(np.abs(grid))

            #if starttot > endtot:
            print("startendPOST", starttot, endtot, N_allzero, N_2zero, N_all)
            #    import pdb; pdb.set_trace()
















    bUseWhorlsToReduce = False
    STARTCOMPRESS = 1 # 20 # 1
    if bUseWhorlsToReduce:


        for iscen in range(grid.shape[0]):
            for ix in range(grid.shape[1]):
                for iy in range(grid.shape[2]):
                    if (ix + iy + t) % 2 != 1:
                        continue # wrong parity  
                    
                    iloop = 0
                    while True:

                        maxout = max(grid[iscen,ix,iy,:])
                        minout = min(grid[iscen,ix,iy,:])
                        outspan = maxout - minout

                        minspan = min([maxout, -minout]) # = 500
                        
                        if outspan < 2 * STARTCOMPRESS:
                            break
                        if iloop >= 3: # consider output arcs of the form [2, 2, 2, -6]; if startcompress is 1, then we can dewhorl t his 3 times.
                            break
                        iloop += 1

                        
                        outgoing_along_xpos = copy.copy(grid[iscen, ix, iy, 0])  
                        outgoing_along_ypos = copy.copy(grid[iscen, ix, iy, 1])
                        outgoing_along_xneg = copy.copy(grid[iscen, ix, iy, 2])
                        outgoing_along_yneg = copy.copy(grid[iscen, ix, iy, 3])

                        outgoing_along_xneg = copy.copy(grid[iscen, ix, iy, 2])
                        outarcs = [outgoing_along_xpos, outgoing_along_ypos, outgoing_along_xneg, outgoing_along_yneg]
                        maxind = rn.choice(GetAllInd(outarcs, maxout))
                        minind = rn.choice(GetAllInd(outarcs, minout))

                        #import pdb; pdb.set_trace() # maxind,minind,outarcs, ix, iy
                        if (maxind - minind) % 2 == 0:
                            

                            # in this case, the two indices are oppositely oriented 
                            #bWhich = rn.choice([True, False])
                            if maxind == 0 or maxind == 2: # implies minind == 2
                                if maxind == 0:
                                    iminspan = minspan * +1
                                elif maxind == 2:
                                    iminspan = minspan * -1

                                post_outarcs_top_sumpart = abs(grid[iscen, ix, iy, 0] - iminspan) + abs(grid[iscen, ix, iy, 2] + iminspan) + abs(grid[iscen, (ix + 1) % grid.shape[1], (iy + 1) % grid.shape[2], 3] + iminspan) + abs(grid[iscen, (ix + 1) % grid.shape[1], (iy + 1) % grid.shape[2], 2] - iminspan) + abs(grid[iscen, (ix - 1) % grid.shape[1], (iy + 1) % grid.shape[2], 0] + iminspan) + abs(grid[iscen, (ix - 1) % grid.shape[1], (iy + 1) % grid.shape[2], 3] + -iminspan)
                                pre_outarcs_top_sumpart = abs(grid[iscen, ix, iy, 0] - 0*iminspan) + abs(grid[iscen, ix, iy, 2] + 0*iminspan) + abs(grid[iscen, (ix + 1) % grid.shape[1], (iy + 1) % grid.shape[2], 3] + 0*iminspan) + abs(grid[iscen, (ix + 1) % grid.shape[1], (iy + 1) % grid.shape[2], 2] - 0*iminspan) + abs(grid[iscen, (ix - 1) % grid.shape[1], (iy + 1) % grid.shape[2], 0] + 0*iminspan) + abs(grid[iscen, (ix - 1) % grid.shape[1], (iy + 1) % grid.shape[2], 3] - 0*iminspan)

                                post_outarcs_bot_sumpart = abs(grid[iscen, ix, iy, 0] - iminspan) + abs(grid[iscen, ix, iy, 2] + iminspan) + abs(grid[iscen, (ix + 1) % grid.shape[1], iy - 1, 1] + iminspan) + abs(grid[iscen, (ix + 1) % grid.shape[1], iy - 1, 2] + -iminspan) + abs(grid[iscen, (ix - 1) % grid.shape[1], iy - 1, 0] + iminspan) + abs(grid[iscen, (ix - 1) % grid.shape[1], iy - 1, 1] + -iminspan)
                                pre_outarcs_bot_sumpart = abs(grid[iscen, ix, iy, 0] - 0*iminspan) + abs(grid[iscen, ix, iy, 2] + 0*iminspan) + abs(grid[iscen, (ix + 1) % grid.shape[1], iy - 1, 1] + 0*iminspan) + abs(grid[iscen, (ix + 1) % grid.shape[1], iy - 1, 2] - 0*iminspan) + abs(grid[iscen, (ix - 1) % grid.shape[1], iy - 1, 0] + 0*iminspan) + abs(grid[iscen, (ix - 1) % grid.shape[1], iy - 1, 1] - 0*iminspan)



                                grid[iscen, ix, iy, 0] += -iminspan
                                grid[iscen, ix, iy, 2] += iminspan
                                if (pre_outarcs_top_sumpart - post_outarcs_top_sumpart) > (post_outarcs_bot_sumpart - pre_outarcs_bot_sumpart) and (pre_outarcs_top_sumpart - post_outarcs_top_sumpart) > 0:
                                    grid[iscen, (ix + 1) % grid.shape[1], (iy + 1) % grid.shape[2], 3] += iminspan
                                    grid[iscen, (ix + 1) % grid.shape[1], (iy + 1) % grid.shape[2], 2] += -iminspan
                                    grid[iscen, (ix - 1) % grid.shape[1], (iy + 1) % grid.shape[2], 0] += iminspan
                                    grid[iscen, (ix - 1) % grid.shape[1], (iy + 1) % grid.shape[2], 3] += -iminspan


                                else:
                                    grid[iscen, (ix + 1) % grid.shape[1], iy - 1, 1] += iminspan
                                    grid[iscen, (ix + 1) % grid.shape[1], iy - 1, 2] += -iminspan
                                    grid[iscen, (ix - 1) % grid.shape[1], iy - 1, 0] += iminspan
                                    grid[iscen, (ix - 1) % grid.shape[1], iy - 1, 1] += -iminspan
                                    
                            elif maxind == 1 or maxind == 3: # implies minind == 2
                                if maxind == 1:
                                    iminspan = minspan * +1
                                elif maxind == 3:
                                    iminspan = minspan * -1



                                post_outarcs_top_sumpart = abs(grid[iscen, ix, iy, 1] - iminspan) + abs(grid[iscen, ix, iy, 3] + iminspan) + abs(grid[iscen, (ix + 1) % grid.shape[1], iy - 1, 1] - iminspan) + abs(grid[iscen, (ix + 1) % grid.shape[1], iy - 1, 2] + iminspan) + abs(grid[iscen, (ix + 1) % grid.shape[1], (iy + 1) % grid.shape[2], 2] + iminspan) + abs(grid[iscen, (ix + 1) % grid.shape[1], (iy + 1) % grid.shape[2], 3] - iminspan)
                                pre_outarcs_top_sumpart = abs(grid[iscen, ix, iy, 1] - 0*iminspan) + abs(grid[iscen, ix, iy, 3] + 0*iminspan) + abs(grid[iscen, (ix + 1) % grid.shape[1], iy - 1, 1] - 0*iminspan) + abs(grid[iscen, (ix + 1) % grid.shape[1], iy - 1, 2] + 0*iminspan) + abs(grid[iscen, (ix + 1) % grid.shape[1], (iy + 1) % grid.shape[2], 2] + 0*iminspan) + abs(grid[iscen, (ix + 1) % grid.shape[1], (iy + 1) % grid.shape[2], 3] - 0*iminspan)

                                post_outarcs_bot_sumpart = abs(grid[iscen, ix, iy, 1] - iminspan) + abs(grid[iscen, ix, iy, 3] + iminspan) + abs(grid[iscen, (ix - 1), iy - 1, 1] - iminspan) + abs(grid[iscen, (ix - 1), iy - 1, 0] + iminspan) + abs(grid[iscen, (ix - 1), (iy + 1) % grid.shape[2], 0] + iminspan) + abs(grid[iscen, (ix - 1), (iy + 1) % grid.shape[2], 3] - iminspan)
                                pre_outarcs_bot_sumpart = abs(grid[iscen, ix, iy, 1] - 0*iminspan) + abs(grid[iscen, ix, iy, 3] + 0*iminspan) + abs(grid[iscen, (ix - 1), iy - 1, 1] - 0*iminspan) + abs(grid[iscen, (ix - 1), iy - 1, 0] + 0*iminspan) + abs(grid[iscen, (ix - 1), (iy + 1) % grid.shape[2], 0] + 0*iminspan) + abs(grid[iscen, (ix - 1), (iy + 1) % grid.shape[2], 3] - 0*iminspan)





                                grid[iscen, ix, iy, 1] += -iminspan
                                grid[iscen, ix, iy, 3] += iminspan
                                if (pre_outarcs_top_sumpart - post_outarcs_top_sumpart) > (post_outarcs_bot_sumpart - pre_outarcs_bot_sumpart) and (pre_outarcs_top_sumpart - post_outarcs_top_sumpart) > 0:

                                    grid[iscen, (ix + 1) % grid.shape[1], iy - 1, 1] += -iminspan
                                    grid[iscen, (ix + 1) % grid.shape[1], iy - 1, 2] += iminspan
                                    grid[iscen, (ix + 1) % grid.shape[1], (iy + 1) % grid.shape[2], 2] += iminspan
                                    grid[iscen, (ix + 1) % grid.shape[1], (iy + 1) % grid.shape[2], 3] += -iminspan
                                else:
                                    grid[iscen, (ix - 1), iy - 1, 1] += -iminspan
                                    grid[iscen, (ix - 1), iy - 1, 0] += iminspan
                                    grid[iscen, (ix - 1), (iy + 1) % grid.shape[2], 0] += iminspan
                                    grid[iscen, (ix - 1), (iy + 1) % grid.shape[2], 3] += -iminspan

                        else:
                            # in this case, the two arcs can be handled by a single whirl
                            if (maxind == 0 and minind == 3) or (maxind == 3 and minind == 0):
                                #checked
                                if (maxind == 0 and minind == 3):
                                    iminspan = minspan * +1
                                elif (maxind == 3 and minind == 0):
                                    iminspan = minspan * -1

                                post = abs(grid[iscen, ix, iy, 0] - iminspan) + abs(grid[iscen, ix, iy, 3] + iminspan) + abs(grid[iscen, (ix + 1) % grid.shape[1], iy - 1, 1] + iminspan) + abs(grid[iscen, (ix + 1) % grid.shape[1], iy - 1, 2] - iminspan)
                                pre = abs(grid[iscen, ix, iy, 0]) + abs(grid[iscen, ix, iy, 3]) + abs(grid[iscen, (ix + 1) % grid.shape[1], iy - 1, 1]) + abs(grid[iscen, (ix + 1) % grid.shape[1], iy - 1, 2])
                                if post < pre:
                                    grid[iscen, ix, iy, 0] += -iminspan
                                    grid[iscen, ix, iy, 3] += iminspan
                                    grid[iscen, (ix + 1) % grid.shape[1], iy - 1, 1] += iminspan
                                    grid[iscen, (ix + 1) % grid.shape[1], iy - 1, 2] += -iminspan

                            elif (maxind == 0 and minind == 1) or (maxind == 1 and minind == 0):
                                # checked
                                if (maxind == 0 and minind == 1):
                                    iminspan = minspan * +1
                                elif (maxind == 1 and minind == 0):
                                    iminspan = minspan * -1

                                post = abs(grid[iscen, ix, iy, 0] - iminspan) + abs(grid[iscen, ix, iy, 1] + iminspan) + abs(grid[iscen, (ix + 1) % grid.shape[1], (iy + 1) % grid.shape[2], 3] + iminspan) + abs(grid[iscen, (ix + 1) % grid.shape[1], (iy + 1) % grid.shape[2], 2] - iminspan)    
                                pre = abs(grid[iscen, ix, iy, 0]) + abs(grid[iscen, ix, iy, 1]) + abs(grid[iscen, (ix + 1) % grid.shape[1], (iy + 1) % grid.shape[2], 3]) + abs(grid[iscen, (ix + 1) % grid.shape[1], (iy + 1) % grid.shape[2], 2])   
                                if post < pre:
                                    grid[iscen, ix, iy, 0] += -iminspan
                                    grid[iscen, ix, iy, 1] += iminspan
                                    grid[iscen, (ix + 1) % grid.shape[1], (iy + 1) % grid.shape[2], 3] += iminspan
                                    grid[iscen, (ix + 1) % grid.shape[1], (iy + 1) % grid.shape[2], 2] += -iminspan

                            elif (maxind == 2 and minind == 3) or (maxind == 3 and minind == 2):
                                if (maxind == 2 and minind == 3):
                                    iminspan = minspan * +1
                                elif (maxind == 3 and minind == 2):
                                    iminspan = minspan * -1

                                post = abs(grid[iscen, ix, iy, 2] - iminspan) + abs(grid[iscen, ix, iy, 3] + iminspan) + abs(grid[iscen, ix - 1, iy - 1, 1] - iminspan) + abs(grid[iscen, ix - 1, iy - 1, 0] + iminspan)
                                pre = abs(grid[iscen, ix, iy, 2]) + abs(grid[iscen, ix, iy, 3]) + abs(grid[iscen, ix - 1, iy - 1, 1]) + abs(grid[iscen, ix - 1, iy - 1, 0])

                                if post < pre:
                                    grid[iscen, ix, iy, 2] += -iminspan
                                    grid[iscen, ix, iy, 3] += iminspan
                                    grid[iscen, ix - 1, iy - 1, 1] += -iminspan
                                    grid[iscen, ix - 1, iy - 1, 0] += iminspan

                            elif (maxind == 1 and minind == 2) or (maxind == 2 and minind == 1):
                                # checked
                                if (maxind == 1 and minind == 2):
                                    iminspan = minspan * -1
                                elif (maxind == 2 and minind == 1):
                                    iminspan = minspan * +1

                                post = abs(grid[iscen, ix, iy, 2] - iminspan) + abs(grid[iscen, ix, iy, 1] + iminspan) + abs(grid[iscen, ix - 1, (iy + 1) % grid.shape[2], 3] + iminspan) + abs(grid[iscen, ix - 1, (iy + 1) % grid.shape[2], 0] - iminspan)
                                pre = abs(grid[iscen, ix, iy, 2] - iminspan) + abs(grid[iscen, ix, iy, 1]) + abs(grid[iscen, ix - 1, (iy + 1) % grid.shape[2], 3]) + abs(grid[iscen, ix - 1, (iy + 1) % grid.shape[2], 0])

                                if post < pre:
                                    grid[iscen, ix, iy, 2] += -iminspan
                                    grid[iscen, ix, iy, 1] += iminspan
                                    grid[iscen, ix - 1, (iy + 1) % grid.shape[2], 3] += iminspan
                                    grid[iscen, ix - 1, (iy + 1) % grid.shape[2], 0] += -iminspan

                                        


                    
                    #incoming_along_xpos = copy.copy(grid[iscen, ix-1, iy, 0])  # note the origin of this flow is in the xneg direction, i.e. ix-1
                    #incoming_along_ypos = copy.copy(grid[iscen, ix, iy-1, 1])
                    #incoming_along_xneg = copy.copy(grid[iscen, (ix+1) % grid.shape[1], iy, 2])
                    #incoming_along_yneg = copy.copy(grid[iscen, ix, (iy+1) % grid.shape[2], 3])






    # now, wipe the previous time step's particles (though we already dit it)
    if t % 2 == 0:
        grid[:, 0::2, 0::2, :] = 0
        grid[:, 1::2, 1::2, :] = 0
    else:
        grid[:, 1::2, 0::2, :] = 0
        grid[:, 0::2, 1::2, :] = 0
    
    modpoolvec = 0 * modpoolvec

    #import pdb; pdb.set_trace()
    if bBrownianRegularization:
        if t % 2 == 0:
            grid_br[:, 0::2, 0::2, :] = 0
            grid_br[:, 1::2, 1::2, :] = 0
        else:
            grid_br[:, 1::2, 0::2, :] = 0
            grid_br[:, 0::2, 1::2, :] = 0



        for iscen in range(grid.shape[0]):
            for ix in range(grid.shape[1]):
                for iy in range(grid.shape[2]):
                    if (ix + iy + t) % 2 != 1:
                        continue # wrong parity  
                    
                    #iloop = 0
                    #if  (t,ix,iy) == (0, 3, 2):
                    #    print("IT")
                    #    import pdb; pdb.set_trace()

                    origarcs = copy.copy(grid[iscen, ix, iy, :])

                    minarcs = np.min(origarcs)
                    maxarcs = np.max(origarcs)

                    #import pdb; pdb.set_trace()
                   
                    while minarcs * maxarcs < 0 and np.min([maxarcs, -minarcs]) > 1:
                        #import pdb; pdb.set_trace()
                        minind = GetRandomLocation(origarcs, minarcs)
                        maxind = GetRandomLocation(origarcs, maxarcs)

                        fillamt = np.min([maxarcs, -minarcs])

                        #the brownian grid will already be added to the reg grid so do NOT add to reg grid here, but do update the copy of the grid
                        #import pdb; pdb.set_trace()
                        grid_br[iscen, ix, iy, maxind] -= fillamt
                        grid_br[iscen, ix, iy, minind] += fillamt
                        origarcs[maxind] -= fillamt
                        origarcs[minind] += fillamt
                        #import pdb; pdb.set_trace()

                        """
                        if maxind == 0:
                            grid_br[iscen, (ix+1)%NDIM, iy, maxind] -= fillamt
                        elif maxind == 1:
                            grid_br[iscen, ix, (iy+1)%NDIM, maxind] -= fillamt
                        elif maxind == 2:
                            grid_br[iscen, (ix-1)%NDIM, iy, maxind] -= fillamt
                        elif maxind == 3:
                            grid_br[iscen, ix, (iy-1)%NDIM, maxind] -= fillamt
                        
                        if minind == 0:
                            grid_br[iscen, (ix+1)%NDIM, iy, minind] += fillamt
                        elif minind == 1:
                            grid_br[iscen, ix, (iy+1)%NDIM, minind] += fillamt
                        elif minind == 2:
                            grid_br[iscen, (ix-1)%NDIM, iy, minind] += fillamt
                        elif minind == 3:
                            grid_br[iscen, ix, (iy-1)%NDIM, minind] += fillamt
                        """


                        minarcs = np.min(origarcs)
                        maxarcs = np.max(origarcs)






def ProcessFile():
  
  global bUseLatticeGasForRand
  global bCplusplus
  global modpoolvec
  global alpha
  
  opts, args = ReadParams()

  alpha = opts.alpha
  #import pdb; pdb.set_trace()
  seed = abs( opts.seed )
  
  rn.seed( seed  )
  bDiscrete = opts.bDiscrete
  
  nplotmod = 1200000 # 100 #20000 #1200000 
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
  NDIM = 8 #64 # 32 # 64 # 4 # 64 # 100 # must be at least 10 or 20 
  
  if NDIM % 2 == 1:
      print("Make the grid length even parity to make it truly systolic")
      import pdb; pdb.set_trace()
  
  
  history_xval = 0  # this is used to track 1 specific particle
  
  NScenarios = 2 # we get smoother distributions if se generate several different gases and average the amplitube over them
                             # in this version, the different scenarios are just the different edges connecting each pair of neighbors
  
  
  if not(bDiscrete):
      NScenarios = 1
  
  halfband = 3
  
  if bDiscrete:
      grid = np.zeros((NScenarios, NDIM, NDIM, 4)).astype("int")
      grid_brown = np.zeros((NScenarios, NDIM, NDIM, 4)).astype("int")
      #probably don't need any of this
      gridlatticegas = np.zeros((NScenarios, NDIM, NDIM, 4)).astype("int") # this will be used to feed the random numbers
      modpoolvec = np.zeros(grid.shape[:-1]).astype("int")
  else:
      grid = np.zeros((NScenarios, NDIM, NDIM, 4))
      grid_brown = np.zeros((NScenarios, NDIM, NDIM, 4)).astype("int")
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
  dellist = []
  for i in range(len(d2screen)):
    if d2screen[i] >= grid.shape[1]: 
        dellist.append(i)
  dellist.reverse()
  for i in dellist:
    del d2screen[i]
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
  



  """
  imgplot = plt.imshow(1 - Flatten(grid)).set_cmap('hot')  # show 1st scenario
  plt.axis('off')
  plt.show()
  
  
  #import pdb; pdb.set_trace()
  for iscr in range(screen.shape[0]):
      screenplot = plt.plot(np.arange(0,NDIM), screen[iscr,:])  # show 1st scenario
  plt.show()
  """

  global abscount
  global reduction_count
  global inc_count
  global reduction_amt
  global inc_amt


  abscount = 0
  reduction_count = 0
  inc_count = 0
  reduction_amt = 0
  inc_amt = 0
 

  #########################################
  #########################################
  ##
  ##  INITIAL CONDITIONS (end)
  ##
  #########################################
  #########################################

  thisparity = 0
  
  
  TMAX = max([50 * 1000 * 1000, nplotmod])
 
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
  
    if False: #t >= 1: # t % 1000 == 0: #  in [0,1] and t > 1:     t > 0: #
       print4grid(grid, t)
       #import pdb; pdb.set_trace()

    if t % refreshtime == 0 and not(bNoWipe):
        #import pdb; pdb.set_trace()
        grid = grid * 0


    if not(bNoWipe):
        #grid[:, d2split, : ,:] = 0 #  
        grid[:, d2split-1, : ,:] = 0 #  
        grid[:, d2split-2, : ,:] = 0 #  
    
    iscr = 0

    #import pdb; pdb.set_trace()
    if True: # t % 10 == 0:  
      if True: # t % pulsetime == 0:          
          if bDiscrete:
            for i in range(grid.shape[0]):
                #import pdb; pdb.set_trace()
                
                if t == 0 or not(bNoWipe):
                    
                    #grid[i, d2split, dA ,:] = randamp(grid[i, d2split, dA ,:].shape, np.cos(2 * np.pi * (t % NPer) / float(NPer))) # 1 #  
                    #grid[i, d2split, -dA ,:] = randamp(grid[i, d2split, -dA ,:].shape, -np.cos(2 * np.pi * (t % NPer) / float(NPer)))# -1 #  


                    #import pdb; pdb.set_trace()
                    #grid[:, d2split, dA ,:] = 1   # 1 #  
                    #grid[:, d2split, dA-2 ,:] = -1 # -1 #  
                    grid = np.zeros(grid.shape).astype("int") # 10000 * np.ones(grid.shape).astype("int")


                    # 1-square
                    #grid[:, 1, 1, 0] = 1
                    #grid[:, 1, 1, 1] = -1
                    #grid[:, 2, 2, 3] = -1
                    #grid[:, 2, 2, 2] = 1

                    # outward 
                    grid[:, 3, 3, 2] = -1
                    grid[:, 3, 3, 3] = 1


                    # 3x3square ping-pong (in the mod version)
                    """
                    grid[:, 0, 0, 0] = 1 
                    grid[:, 0, 0, 1] = -1
                    grid[:, 0, 2, 3] = 5
                    grid[:, 0, 2, 1] = -5
                    grid[:, 1, 3, 2] = 1
                    grid[:, 1, 3, 0] = -1
                    grid[:, 3, 3, 2] = 1
                    grid[:, 3, 3, 3] = -1
                    grid[:, 3, 1, 1] = 5
                    grid[:, 3, 1, 3] = -5
                    grid[:, 2, 0, 0] = 5
                    grid[:, 2, 0, 2] = -5
                    """

                    
                    #import pdb; pdb.set_trace()
                    grid[:, 0, 0, 0] = 1 
                    grid[:, 1, 1, 1] = -1
                    grid[:, 0, 2, 3] = 5
                    grid[:, 0, 2, 1] = -5
                    grid[:, 1, 1, 2] = 1
                    grid[:, 1, 3, 0] = -1
                    grid[:, 3, 1, 2] = 1
                    grid[:, 1, 1, 3] = -1
                    grid[:, 3, 1, 1] = 5
                    grid[:, 3, 1, 3] = -5
                    grid[:, 2, 0, 0] = 5
                    grid[:, 2, 0, 2] = -5
                    


                  
                    """
                    grid_brown[:, 0, 0, 0] = 1 
                    grid_brown[:, 1, 1, 1] = -1
                    grid_brown[:, 0, 2, 3] = 5
                    grid_brown[:, 0, 2, 1] = -5
                    grid_brown[:, 1, 1, 2] = 1
                    grid_brown[:, 1, 3, 0] = -1
                    grid_brown[:, 3, 1, 2] = 1
                    grid_brown[:, 1, 1, 3] = -1
                    grid_brown[:, 3, 1, 1] = 5
                    grid_brown[:, 3, 1, 3] = -5
                    grid_brown[:, 2, 0, 0] = 5
                    grid_brown[:, 2, 0, 2] = -5
                    """
                    
                    #import pdb; pdb.set_trace()
                    #grid_brown = 0 * grid_brown
                    #grid_brown[0, 1, 1, 1] = 1

                    # 3x2-square
                    #grid[:, 0, 0, 0] = 1 
                    #grid[:, 0, 0, 1] = -1
                    #grid[:, 0, 2, 3] = 1
                    #grid[:, 0, 2, 0] = -1
                    #grid[:, 2, 2, 2] = 1
                    #grid[:, 2, 2, 0] = -1
                    #grid[:, 3, 1, 1] = 1
                    #grid[:, 3, 1, 3] = -1
                    #grid[:, 2, 0, 0] = 1
                    #grid[:, 2, 0, 2] = -1



                    #grid[:, 1, 3, 1] = -1
                    #grid[:, 1, 3, 0] = 1
                    #grid[:, 2, 2, 1] = -1
                    #grid[:, 2, 2, 3] = 1
                    #grid[:, 2, 0, 1] = -1
                    #grid[:, 2, 0, 2] = 1
                    #grid[:, 0, 0, 0] = -1


                    #grid[:, 2, 2, 2] += 1
                    #grid[:, 2, 2, 3] += -1
                    
                    #another curve
                    #grid[:, 0, 0 , 1] = 1 
                    #grid[:, 0, 2, 3] = -1
                    #grid[:, 0, 2, 0] = 1
                    #grid[:, 1, 3, 3] = -1
                    #grid[:, 1, 3, 0] = 1
                    #grid[:, 2, 2, 1] = -1
                    #grid[:, 2, 2, 3] = 1
                    #grid[:, 2, 0, 1] = -1
                    #grid[:, 2, 0, 2] = 1
                    #grid[:, 0, 0, 0] = -1

                    
                    #grid[:, NDIM//2:, NDIM//2:, :] += 1
                    #grid[:, NDIM//2:, NDIM//2:, :] += -1

                    print4grid(grid, t)
                    #import pdb; pdb.set_trace()
                    # I will generate PAIRS of particles according to a cos distribution; I suspect doing by pairs is cheating but 
                    # so be it -- this is after all little more than a gedankenexperiment with code
                    """
                    preR1 = np.cos(2 * np.pi * (t % NPer) / float(NPer))
                    r1 = 0
                    if np.random.rand() < np.abs(preR1):
                        r1 = np.sign(preR1)
                    #import pdb; pdb.set_trace()
                    if r1 != 0: #for iscr in range(grid.shape[0]):
                        # first boundary condition
                        #import pdb; pdb.set_trace()
                        grid[0, d2split, dA ,0] = r1 #randamp(grid[:, d2split, dA ,0].shape, np.cos(2 * np.pi * (t % NPer) / float(NPer))) # 1 #  
                        #grid[0, d2split, dA+2 ,2] = -1 #randamp(grid[:, d2split, dA ,0].shape, np.cos(2 * np.pi * (t % NPer) / float(NPer))) # 1 #  
                        grid[0, d2split, -dA ,0] = -r1 #randamp(grid[:, d2split, dA ,0].shape, np.cos(2 * np.pi * (t % NPer) / float(NPer))) # 1 #  
                        
                        
                        # second opposing boundary condition
                        #grid[iscr, d2split, dA-2 ,0] = -1 #grid[iscr, d2split, -dA ,0] = -1 #randamp(grid[:, d2split, -dA ,0].shape, -np.cos(2 * np.pi * (t % NPer) / float(NPer)))# -1 #  
                    
                    
                    
                    
                    if iscr == 0:
                    
                        if grid.shape[1] >= 8:
                            try:
                                print("here is the sum of amplitudes in LHS of grid (scenario 0): " + str(np.sum(grid[0, 0, 0:NDIM//2])))
                                print("here is the sum of amplitudes in RHS of grid (scenario 0): " + str(np.sum(grid[0, 0, 9:])))
                            except:
                                import pdb; pdb.set_trace()
                    
                    """
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
          wavegridupdate(t, grid, grid_brown) #, #fermpenalty, n2penalty, bFermion = False, 4)
          #import pdb; pdb.set_trace()
            
    else:

        classicalwavefloatingpt(grid)
    
    if bEnableSingleParticleTracking:
        tracker.postupdate(grid)
      
      
    #if t < 100 and bDiscrete and grid.shape[0] > 1:
    #    #import pdb; pdb.set_trace()
    #    grid[1,:,:,:] = grid[1,:,:,:] * 0
    #    flatten2dgridandcsv(grid, True, t, "flattenedgrid.csv")
      
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
        
       
        
    if t % add2HistoryIncrement == 0:
        #import pdb; pdb.set_trace()
        tHisTArr.append(t)
        absTotArr.append(np.sum(np.abs(grid[:,:,:,:])))
        l2Arr.append(np.sum(grid[:,:,:,:]*grid[:,:,:,:]))
        maxarr.append(np.max(np.abs(grid)))
        gridct = float(np.prod(list(grid.shape)))/2.0
        fermipct = ((np.sum(np.abs(grid) <= 1) - gridct)) / float(gridct)
        ferminoic.append(fermipct)
        if len(maxarr) > 1:
            maxarr[-1] = copy.copy(maxarr[-1] - maxarr[0])
            #import pdb; pdb.set_trace()
        if t % add2HistoryIncrement == 0:
           import pickle
           fname = opts.picklefile
           if len(fname) > 4 and fname[-4:] != '.pkl':
              fname = fname + '.pkl'
           f = open(fname, 'wb')

           pickle.dump((tHisTArr, absTotArr, l2Arr, maxarr, ferminoic), f)
           f.close()



    sumgrid = np.sum(grid)
    maxgrid = np.max(grid)
    absgrid = np.sum(np.abs(grid))



    sumgrid2 = np.sum(grid[0,:,:,:] * grid[0,:,:,:])
    sumbraket = np.sum(grid[0,:,:,:] * grid[1,:,:,:])
    sumamp2 =  np.sum(np.sum(grid[:, :, :, :], 3) * np.sum(grid[:, :, :, :], 3))



    sumgrid_brown = np.sum(grid_brown)
    maxgrid_brown = np.max(grid_brown)
    absgrid_brown = np.sum(np.abs(grid_brown))    

    sumgrid_brown2 = np.sum(grid_brown[0,:,:,:] * grid_brown[0,:,:,:])
    if sumgrid_brown2 < 0:
        import pdb; pdb.set_trace()





    
    if True: #t == 0:
        emamax =  maxgrid
        emaabs = absgrid

        beta = 0 # 0.01
        betaalt = 1 - beta
    else:
        emamax = beta * maxgrid + betaalt * emamax
        emaabs = beta * absgrid + betaalt * emaabs

    minabs = 1000 * 1000 * 1000 * 1000 * np.ones((NScenarios, )).astype("int")
    #if t >= 1000:
    #    import pdb; pdb.set_trace()
    for iscen in range(grid.shape[0]):
        for ix in range(grid.shape[1]):
            for iy in range(grid.shape[2]):
                if (ix + iy) % 2 == 0:
                    continue

                minabs[iscen] = int(np.min([minabs[iscen], np.min(np.abs(grid[iscen, ix, iy, :]))]))
                #if t >= 1000:
                #   import pdb; pdb.set_trace()

                
    gridct = float(np.prod(list(grid.shape)))/2.0
    fermipct = ((np.sum(np.abs(grid) <= 1) - gridct)) / float(gridct)
    print2(str(t))

    #import pdb; pdb.set_trace()
    print2("sum     " + str(sumgrid) + " " + "%5.1f" % (np.sum(np.abs(grid)),) + " " + "%5.1f" % (emaabs/gridct,) + " %5.1f" % (sumgrid2/gridct,) + " %5.1f" % (sumbraket/gridct,) + " fermipct %5.4f" % (np.sqrt(fermipct),)) # + " " + str(np.sum(grid[:,0::2,0::2,:])) + " " + str(np.sum(grid[:,0::2,1::2,:]))  )
    print2("sumXmom " + str(np.sum(grid[:,:,:,Comp000] - grid[:,:,:,Comp180])) + " sumYmom " + str(np.sum(grid[:,:,:,Comp090] - grid[:,:,:,Comp270])) +  ' max ' + str(maxgrid) + ' min ' + str(np.min(grid)) +  ' 90pct ' + str(np.round(np.percentile(grid, 90))) + ' 10pct ' + str(np.round(np.percentile(grid, 10))) + ' mnab ' + str(minabs[0]) + ' ' + str(minabs[1]) )
 
    print2("sumBR     " + str(sumgrid_brown) + " " + str(absgrid_brown) + " " + " %5.1f" % (sumgrid_brown2/gridct,)) # + " " + str(np.sum(grid[:,0::2,0::2,:])) + " " + str(np.sum(grid[:,0::2,1::2,:]))  )
    #print2("sumXmomBR " + str(np.sum(grid_brown[:,:,:,Comp000] - grid[:,:,:,Comp180])) + " sumYmom " + str(np.sum(grid[:,:,:,Comp090] - grid[:,:,:,Comp270])) +  ' max ' + str(maxgrid) + ' min ' + str(np.min(grid)) +  ' 90pct ' + str(np.round(np.percentile(grid, 90))) + ' 10pct ' + str(np.round(np.percentile(grid, 10))) + ' mnab ' + str(minabs[0]) + ' ' + str(minabs[1]) )
 
    #print2("sumYmom " + str(np.sum(grid[:,:,:,Comp090] - grid[:,:,:,Comp270])) +  ' max ' + str(np.max(grid)) + ' min ' + str(np.min(grid)))
    #print2("sumYmom " + str(np.sum(grid[:,:,:,Comp090] - grid[:,:,:,Comp270])) +  ' max ' + str(np.max(grid)) + ' min ' + str(np.min(grid)) +  ' 90pct ' + str(np.round(np.percentile(grid, 90))) + ' 10pct ' + str(np.round(np.percentile(grid, 10)))                    )
    #if t >= 4000: # t % 1000 in [0,1] and t > 1:  
    #   print2grid(grid, t)

    print("inc_amb/incct/redamt/redct", inc_amt, inc_count, reduction_amt, reduction_count)


    table_data = [ [ (' %5d' % (x,)) for x in list(np.arange(-NHistogram, NHistogram, 2))], 
                  [ (' %5d' % (x,)) for x in list(HistPartChange)] ]
    for row in table_data:
        print(("{: >6}" * NHistogram).format(*row))
    #print("partchg", ('%3d '*NHistogram) % tuple(np.arange(-NHistogram, NHistogram, 2)) )
    #print("partchg", ('%3d '*NHistogram) % tuple(HistPartChange) )
    #print("spanchg", ('%3d '*NHistogram) % tuple(HistSpanChange) )

    if t == 0:
        emamaxarr = []
        emaabsarr = []
        rmsarr = []
        #ferminoic = [] # gives percentage of nodes where abs(amplitude) < 2
        tarr = []
        gridsnapshot = []
        tarr2 = []
        
    if t % 100 == 0 and t > 0:
        filename1 = "blah.csv"
        #import pdb; pdb.set_trace()
        emamaxarr.append(emamax)
        emaabsarr.append(emaabs)
        rmsarr.append(np.sqrt(sumamp2)) # (np.sqrt(np.abs(np.sum(grid[:,:,:,:]*grid[:,:,:,:]))))
        #gridct = float(np.prod(list(grid.shape)))/2.0
        #fermipct = ((np.sum(np.abs(grid) <= 1) - gridct)) / float(gridct)
        #ferminoic.append(fermipct)
        tarr.append(t)

        emaarrnp = np.vstack([np.array(tarr), np.array(emaabsarr), np.array(emamaxarr), rmsarr])

        

        np.savetxt(filename1, emaarrnp.T, fmt='%5.1f', delimiter=',')
    if t % 10000 == 0 and t > 0:
        filename2 = "blah2.pkl"
        #import pickle
        if t == 10000:
            
            fp = open(filename2, "wb")
            gridsnapshot.append( [t, grid[0,:,:,0], grid[0,:,:,1], grid[0,:,:,2], grid[0,:,:,3]] )
            pickle.dump(gridsnapshot, fp)
            fp.close()
        else:
            
            fp = open(filename2, "rb")
            oldgridsnapshot = pickle.load(fp)
            fp.close()
            fp = open(filename2, "wb")
            oldgridsnapshot.append( [t, grid[0,:,:,0], grid[0,:,:,1], grid[0,:,:,2], grid[0,:,:,3]]) 
            pickle.dump(oldgridsnapshot, fp)
            fp.close()
            


    
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
      
      
      
      bSkipPlot = True # bDiscrete
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
        


      if grid.shape[1] > 4:
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








"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

xarr = np.genfromtxt('blah.csv', delimiter=',')
plt.plot(xarr[:,0], xarr[:,3])
plt.show()


import pickle
filename2 = "blah2.pkl"


fp = open(filename2, "rb")
oldgridsnapshot = pickle.load(fp)
fp.close()

whichshot = -100
ogs = oldgridsnapshot
grid = np.zeros((1,4,4,4)).astype("int")
for idir in range(4):
    grid[0,:,:,idir] = ogs[whichshot][idir+1]

print(grid[0,:,:,0])
print(" ")
print(grid[0,:,:,1])
print(" ")
print(grid[0,:,:,2])
print(" ")
print(grid[0,:,:,3])





"""



"""


import pickle
import numpy as np
import matplotlib.pyplot as plt
f = open("trash.pkl", 'rb') # open("polygon/bhpropagation_varreduction_recordY.pkl", 'rb')
outtup = pickle.load(f)
f.close()
len(outtup)
outtup[0]


plt.plot(outtup[0], outtup[1])  # used to be 1 instead of -1
plt.show()


plt.plot(outtup[0], outtup[1]) # used to be 2 instead of -1
plt.show()


plt.plot(outtup[0], np.sqrt(np.array(outtup[1]))) # used to be 2 instead of -1
plt.show()


plt.plot(outtup[0], outtup[1]) # used to be 2 instead of -1
plt.plot(outtup[0], np.sqrt(np.array(outtup[1]))) # used to be 2 instead of -1
plt.show()


a = np.log(np.array(outtup[0][1:])) 
b = np.log(np.array(outtup[1][1:])) # used to be 1 instead of -1

plt.plot(a, b)
plt.legend()
plt.show()

if True:
    x = a[:10000]
    y = b[:10000]
    slope, intercept = np.polyfit(x, y, 1)
    y_fit = slope * x + intercept
    print(f"Slope: {slope:.2f}, Intercept: {intercept:.2f}")
    plt.scatter(x, y) #, label='Original Data')
    #plt.plot(x, y_fit, color='red', label=f'Best Fit Line (y = {slope:.2f}x + {intercept:.2f})')
    plt.plot(x, y_fit, color='red', label=f'Best fit line (y ~  {slope:.2f}x)')
    plt.xlabel('log(Timesteps)')
    plt.ylabel('log(N_particles + N_antiparticles)')
    plt.title("B-H Propagation with Variance Reduction")
    plt.legend()
    plt.show()

from scipy.stats import linregress
linregress(a, b)


a = np.array(outtup[0][1:])
b = np.array(outtup[1][1:])

plt.plot(a, b)
plt.show()

IF SLOPE IS 1 (which it is in case of pure Brownian Huygens) then the growth is LINEAR








SOME COmPARISON DATASETS WITH A SPECIFIC EXPONENTIAL GROWTN

import copy
import numpy as np
N = 200
lam = 0.2
x = np.array(np.arange(1,N+1)).astype("float") 
y = copy.copy(x)
y[0] = 1.0


#for i in range(1,N):
#    y[i] = float(np.sqrt(float(x[i]*x[i])))

for i in range(N):
    y[i] = float(np.power(float(np.log(x[i])), 5))

y[:7]

plt.plot(np.log(x), np.log(y))  # used to be 1 instead of -1
plt.show()


plt.plot(x, np.log(y))  # used to be 1 instead of -1
plt.show()


a = np.log(np.array(x))
b = np.log(np.array(y))



from scipy.stats import linregress
linregress(a[1:], b[1:])









# python  ~/quant/latticegas_youngsdoubleslitexp_partialpooling_inbetween_bernouilli_modpool.py   --discrete --seed 24 --picklefile trash.pkl --alpha 1.0

"""