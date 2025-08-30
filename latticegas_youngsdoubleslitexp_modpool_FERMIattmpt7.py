#!/cygdrive/c/Python27/python


print("""


Oct-2024
Instead of minimizing the number of particles, we instead maximize
the percentage of edges where Pauli Exclusion Principle is satisfied.

Attempt 1(ended in failure): don't stop generating particle/antiparticle
pairs so as to make all outgoing particles
of one sign (or zero). Keep doing it until all outgoing edges satisfy
Pauli Exclusion principle.

Again, this fails -- for instance, we encounter the outgoing arcs
(2, -2, 3, 3). No matter how many particle/antiparticle pairs you
use to annihilate stuff, what's left will not be able to satisfy PEP.

Ergo, move on to 

Attempt 2: allow  multiple particles to be in same edge so long as they
are directly related (i.e. parent and child, and so forth) -- in that case
the multiple particles are just a weight of a given path, more so than
multiple fermions in a given state. 

In this approach, the randomness is just a permutation 
of incoming arcs into outgoing arcs
with the 2 brownian particles of a given arc being sent to the arc 
dictated by the permutation. There will be a problem if the arc and the
permuted arc have particles that are the opposite sign, because in that 
case, the Brownian particles will be piling onto the arc that also has
Huygens particles of the same sign -- violating PEP. 

But luckily for me, I've already put together the Brownian whorl 
enhancement, which conveniently annihilates until all arcs coming into
a node are of the same sign (or zero). 

We also must, for the same reason, be wary about brownian paths -- they
have to pick a step that lands them into an oppositely signed (or zero)
arc. If no such arc exists in a given  node that they find themselves in,
then their oppositely signed subsequent step must be in the opposite 
direction it came (i.e. the loop stays more or less stationary for a step 
or more, until it can find an oppositely signed path to move into).

It seems sketchy, but it's doable, and it's all I got.


UPDATE: Cannot do Brownian Whorls as planned -- the "closing" of any loop epends on the
relative amplitude distributions, I think. If there's a large positive amplitude, then it
will preferantly result in a loop closing somewhere where things are negative. However,
that won't work. The loop closings have to be random



""")

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


PROBLEM: for a large amplitude with net value of zero, the wave fizzles
away to nothingness (i.e. no particles anywhere), which is clearly wrong.

I think the Brownian whorl approach unavoidably mixes positive amplitude
areas with negative amplitude areas to result in zero


Jul-2025

Found the problem with this Fermi approach. Because it "shaves" (i.e. uses Brownian-heat-bath particles) to rectify INCOMING arcs 
as well as outgoing, this thing is just wrong. E.g. consider the incoming case [1,-1,1,-1]; in that case, the shaving
process changes the arcs into the zero vector (but also creates 2 incoming pairs of Brownian particle/antiparticle pairs). However,
in creating output flows, thes Brownian particle/antiparticle pairs just get wiped.

I.e. you've taken a perfectly legit incoming flow and basically just annihilated it, and there's no way to recover from that.

So the updated approach is to leave the incoming particles alone. In the case where incoming arcs are of opposite sign, you still
have to shave the Brownian back-propagatign particle (we choose to shave those instead of the back-propagating Huygens particle
because it leads to a greater reduction in total number of particles), and we also might need to shave further if the
total number of particles keeps growing, but both those shaves will only affect the outgoing arcs. We will henceforth
give up on the notion of trying to adjust or manipulate incoming arcs.

""")



import sys
import os
import optparse
from copy import copy

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
    
    flines = flines + copy(flines)

    
    
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
  
    old_i = copy(self.i)
    old_j = copy(self.j)
    old_dir = copy(self.dir)
  
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
      self.i = copy(old_i)
      self.j = copy(old_j)
      self.dir = copy(old_dir)
    

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
           
    def inflow_coords(dirn, i, j, Ndir):
        #Ndir = 4 
        if dirn == 0:
            return (i-1, j, 0)
        elif dirn == 1:
            return (i, j-1, 1)
        elif dirn == 2:
            return ((i + 1) % Ndir, j, 2)
        elif dirn == 3:
            return (i, (j + 1) % Ndir, 3)

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
            #try:
            shave = np.random.multinomial(-negportion, fraclist)    
            #except:
            #print("pos")        
            return [poslist[i] - shave[i] for i in range(N)]
        elif -negportion > posportion:
            fraclist = -np.array(neglist)/-negportion # all non-netative
            #try:
            shave = np.random.multinomial(posportion, fraclist)
            #except:
            #import pdb; pdb.set_trace()
            return [int(neglist[i] + shave[i]) for i in range(N)]

    
    def NetOutStrict(preinlist):
        """ Unlike NetOut(), this version
        always returns a list whose components are all the same sign (or else zero)
        """

        def intPosIsBigger(inlist):
            sumpos = 0
            sumneg = 0
            for i in inlist:
                if i > 0:
                    sumpos += i
                elif i < 0:
                    sumneg += i
            if sumpos == -sumneg:
                return 0, sumpos, sumneg
            if sumpos > -sumneg:
                return 1, sumpos, sumneg
            return -1, sumpos, sumneg


        def SelectIthBall(xarr, i):
            "gling from left to right, pick the ith particle in the array, assuming they're all the same sign"
            
            cumul = 0
            for iic, ic in enumerate(xarr):
                if cumul + ic > i:
                    return iic
                cumul += ic
            
            return len(xarr) + 1
        

        #if tuple(preinlist) == (0,  0, -2,  2):
        #    import pdb; pdb.set_trace()

        N = len(preinlist)
        signflag, sumposlist, sumneglist = intPosIsBigger(preinlist)

        if signflag == 0:
            return np.zeros((N,)).astype("int")

        inlist = copy(preinlist)
        if signflag == -1:
            for ii, i in enumerate(inlist):
                inlist[ii] = -i
            sumposlist, sumneglist = (-sumneglist, -sumposlist)

        posind = []
        posamt = []

        for ii,i in enumerate(inlist):
            if i > 0:
                posind.append(ii)
                posamt.append(i)
              
        if sumposlist > -sumneglist:
            sumneglistleft = copy(sumneglist)
            for i in range(-sumneglist):
                if len(posamt) == 1:
                    posamt[0] -= -sumneglistleft
                    break
                    
                irand = rn.randint(0, sumposlist-1)
                whichcol = SelectIthBall(posamt, irand)

                posamt[whichcol] -= 1
                sumposlist -= 1
                sumneglistleft += 1

                if posamt[whichcol] <= 0:
                    np.delete(posamt, whichcol)
                    np.delete(posind, whichcol)
        else:
            print("ERROR in NetOutStrict")
            import pdb; pdb.set_trace()

        myretlist = np.zeros((N,)).astype("int")
        for ic, iamt in zip(posind, posamt):
            myretlist[ic] = iamt

        if signflag == -1:
            myretlist = - myretlist

        return myretlist
        


  
  
    def LumpyNetOut(inlist):
        """
        all the balls in a given slot must be lumped together when netting out
        """
        N = len(inlist)
        posportion = 0
        negportion = 0
        poslist = []
        neglist = []

        pospositionlist = []
        negpositionlist = []

        for ii, iflow in enumerate(inlist):
            if iflow > 0:
                poslist.append(iflow)
                neglist.append(0)                
                posportion += iflow
                pospositionlist.append(ii)
            elif iflow < 0:
                neglist.append(iflow)
                poslist.append(0)
                negportion += iflow
                negpositionlist.append(ii)
            else:
                poslist.append(0)
                neglist.append(0)

        
        if posportion == -negportion:
            return [0] * N

        bDone = False
        iloop = -1
        while not(bDone):
            iloop += 1
            if iloop >= N:
                print("Loop is too long here!")
                import pdb; pdb.set_trace()
            copypospositionlist = copy(pospositionlist)
            copynegpositionlist = copy(negpositionlist)
            rn.shuffle(copypospositionlist)
            rn.shuffle(copynegpositionlist)
            if len(pospositionlist) <= len(negpositionlist):
                for i in range(len(copypospositionlist)):
                    ipos = pospositionlist[i]
                    ineg = negpositionlist[i]
                    if inlist[ipos] >= -inlist[ineg]:
                        inlist[ipos] += inlist[ineg]
                        inlist[ineg] = 0
                        negpositionlist.remove(ineg)
                    else:
                        inlist[ineg] += inlist[ipos]
                        inlist[ipos] = 0    
                        pospositionlist.remove(ipos)
            else:                
                for i in range(len(negpositionlist)):
                    ipos = pospositionlist[i]
                    ineg = negpositionlist[i]
                    if inlist[ipos] >= -inlist[ineg]:
                        inlist[ipos] += inlist[ineg]
                        inlist[ineg] = 0
                        negpositionlist.remove(ineg)
                    else:
                        inlist[ineg] += inlist[ipos]
                        inlist[ipos] = 0    
                        pospositionlist.remove(ipos)
            
            if np.min([len(pospositionlist), len(negpositionlist)]) == 0:
                return inlist


    def ModPool(inlist):

        brownian_pool = 2 * np.sum(inlist)
        alpha = 1

        # first, the Huygens portion 

        incoming_along_xpos = -inlist[0]
        incoming_along_ypos = -inlist[1]
        incoming_along_xneg = -inlist[2]
        incoming_along_yneg = -inlist[3]


        if brownian_pool != 0:

            if False and alpha != 0.0 and alpha != 1.0: # the general "mixed case", handles every case of alpha,
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

                posdeviate = np.random.multinomial(pospool, [1/float(N_arc)]*N_arc)
                negdeviate = np.random.multinomial(-negpool, [1/float(N_arc)]*N_arc) 

                outdeviatemod = posdeviate - negdeviate

                incoming_along_xpos += outdeviatemod[0] # + divpool
                incoming_along_ypos += outdeviatemod[1] # + divpool
                incoming_along_xneg += outdeviatemod[2] # + divpool
                incoming_along_yneg += outdeviatemod[3] # + divpool
            
            else: # alpha == 1 (the "modpool" case for which sum(abs(grid)) grows like sqrt(t) )

                # in this case, we don't need a separate pospool and negpool, we just mod out the brownian (i.e. sum of neg and pos) pool

                abspool = np.abs(brownian_pool)
                sgnpool = np.sign(brownian_pool)
                modpool = abspool % N_arc
                divpool = (abspool // N_arc) # this (times sgnpool) is a scalar value that gets added to every output arc


                # if modpool is > N_arc//2, then we should change things so we only have to generate N_arc - modpool random variables.
                NegateModpoolFactor = 1
                if modpool > N_arc // 2:
                    NegateModpoolFactor = -1
                    modpool = N_arc - modpool
                    divpool += 1

                # this is the randomized "leftover" part
                #outdeviatemod = sgnpool * np.random.multinomial(modpool, [1/float(N_arc)]*N_arc)

                # this will work better for fermions in that the modpool portion will be distributed to different arcs
                # i.e. they will be spread out and never lump up
                outdeviatemod = np.zeros((N_arc,)).astype("int")
                modpoolchoices = np.random.choice(N_arc, modpool, False, None)
                for i in modpoolchoices:
                    outdeviatemod[i] += sgnpool * NegateModpoolFactor

                divpool *= sgnpool
                incoming_along_xpos += outdeviatemod[0] + divpool
                incoming_along_ypos += outdeviatemod[1] + divpool
                incoming_along_xneg += outdeviatemod[2] + divpool
                incoming_along_yneg += outdeviatemod[3] + divpool




        return [incoming_along_xpos, incoming_along_ypos, incoming_along_xneg, incoming_along_yneg]


    def ModPoolJustBrownian(inlist):
        # THIS IS NOT BROWNIAN-HUYGENS. It just pools incoming amps, and (assuming the result is a multiple of NArcs)
        # spews the same amplitude in all outgoing arcs. If the modpoos isn't a multiple of NArcs, the remaining bits
        # get assigned to random output arcs

        brownian_pool = np.sum(inlist) # Note this is not multiplied by 2
        alpha = 1

        #Note there is no Huygens contribution
        # first, the Huygens portion 
        incoming_along_xpos = 0 # -inlist[0]
        incoming_along_ypos = 0 # -inlist[1]
        incoming_along_xneg = 0 # -inlist[2]
        incoming_along_yneg = 0 # -inlist[3]


        if brownian_pool != 0:

            if False and alpha != 0.0 and alpha != 1.0: # the general "mixed case", handles every case of alpha,
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

                posdeviate = np.random.multinomial(pospool, [1/float(N_arc)]*N_arc)
                negdeviate = np.random.multinomial(-negpool, [1/float(N_arc)]*N_arc) 

                outdeviatemod = posdeviate - negdeviate

                incoming_along_xpos += outdeviatemod[0] # + divpool
                incoming_along_ypos += outdeviatemod[1] # + divpool
                incoming_along_xneg += outdeviatemod[2] # + divpool
                incoming_along_yneg += outdeviatemod[3] # + divpool
            
            else: # alpha == 1 (the "modpool" case for which sum(abs(grid)) grows like sqrt(t) )

                # in this case, we don't need a separate pospool and negpool, we just mod out the brownian (i.e. sum of neg and pos) pool

                abspool = np.abs(brownian_pool)
                sgnpool = np.sign(brownian_pool)
                modpool = abspool % N_arc
                divpool = (abspool // N_arc) # this (times sgnpool) is a scalar value that gets added to every output arc

                # this is the randomized "leftover" part
                outdeviatemod = sgnpool * np.random.multinomial(modpool, [1/float(N_arc)]*N_arc)

                divpool *= sgnpool
                incoming_along_xpos += outdeviatemod[0] + divpool
                incoming_along_ypos += outdeviatemod[1] + divpool
                incoming_along_xneg += outdeviatemod[2] + divpool
                incoming_along_yneg += outdeviatemod[3] + divpool

        return [incoming_along_xpos, incoming_along_ypos, incoming_along_xneg, incoming_along_yneg]

    def ModPoolJustBrownian(inlist):


        N_arc = len(inlist)
        brownian_pool = np.sum(inlist)

        abspool = np.abs(brownian_pool)
        sgnpool = np.sign(brownian_pool)
        modpool = abspool % N_arc
        divpool = (abspool // N_arc) # this (times sgnpool) is a scalar value that gets added to every output arc

        outdeviatemod = sgnpool * np.random.multinomial(modpool, [1/float(N_arc)]*N_arc)

        # this is the randomized "leftover" part
        outdeviatemod = sgnpool * np.random.multinomial(modpool, [1/float(N_arc)]*N_arc)

        divpool *= sgnpool

        myres = outdeviatemod + divpool

        return myres





    def ShuffleInToOut(a):
        if len(a) == 4:
            b = (a[2], a[3], a[0], a[1])
        elif len(a) == 6:
            b = (a[3], a[4], a[5], a[0], a[1], a[2])
        
        return list(b)

    def PermuteAccordingToInputPerm(xarr, permut):
        outarr = np.zeros(xarr.shape).astype("int")
        for i in range(len(outarr)):
            outarr[permut[i]] = xarr[i]
        return outarr

    

    def FermiAddScalar(A, d, B):
        if A * d < 0:
            print("ERROR in FermiAddScalar: A and d must be of same sign")
            import pdb; pdb.set_trace()
        
        if A * B > 0:
            if np.abs(A + d + B) > np.abs(A):
                return A, d + B
            #elif (A >= 0 and A + d + B > 0) or (A < 0 and A + d + B < 0):
            #    return A + d + B, 0
            else: # in this case, there is a change of sign
                return A + d + B, 0
        else:
            if np.abs(B) > np.abs(d):
                return A + (B + d), 0
            else:
                return A, d + B

    def FermiAdd(A, displaced, B):
        myret = np.zeros(A.shape).astype("int")
        myretdisplaced = np.zeros(A.shape).astype("int")

        
        for i in range(len(A)):
            a, b = FermiAddScalar(A[i], displaced[i], B[i])
            myret[i] = a
            myretdisplaced[i] = b

        return myret, myretdisplaced





    def FermiOut(inputarcarr, permut):
        """
        return the rectified wave, the displacements (which need not sum up to zero)
        and the displacements plus offsets (i.e. the heat-bath brownian pair injection), which always sums up to zero
        since it is particle/antiparticle pairs
        """


        def PartnerArc(i_first, N):
            a = np.arange(N)
            a.remove(i_first)
            return rn.choice(a)


        N = len(inputarcarr)
        
        #def BHInToOut(xarr, permut):
        #    huyg = -ShuffleInToOut(xarr)
        #    brown = 2 * ShuffleInToOut(PermuteAccordingToInputPerm))

        #if np.sum(np.abs(inputarcarr)) > 0:
        #    print(np.sum(np.abs(inputarcarr)))

        unmatched = np.array(NetOut(inputarcarr)).astype("int")
        matched = inputarcarr - unmatched

        # the matched dynamics are trivial: input = -output (i.e. the incoming and outgoing "bounce" against particles with opposite sign)


        displaced = np.zeros(unmatched.shape).astype("int")
        # displaced will mean shunted into heath-bath mode because they convlict with PEP rules. Therefore, they must be offset
        # it is the sum of two shuntings (one that avoids violating PEP, and one to reduce the total number of particles)
        brownout = copy(unmatched)
        huygout = -unmatched # note that for rectangular lattices you will need an extra shuffle
        # to orient the arcs properly (e.g. the positive x direction for incoming is 0, but for outgoing the positive x direction is 2)
        # but we will do that right before we update the grid with what we're calculating here
        for i, ip in zip(np.arange(len(permut)), permut):
            # note that for any i such that i = permut[i], the amplitudes are automatically opposite sign, so that's fine
            brownout[ip] = unmatched[i] * 2
            if huygout[ip] * brownout[ip] > 0:
                displaced[ip] = -huygout[ip]
                huygout[ip] = 0

        if tuple(inputarcarr) == (1 , 0 ,-5, -1):
            import pdb; pdb.set_trace()
        if np.sum(inputarcarr) == 0 and np.sum(np.abs(inputarcarr)) != 0:
            import pdb; pdb.set_trace()
     

        simple_fermicorrected_offset_displaced = brownout + huygout # this is composed just of the unmatched amps


        
        bSimplify = True
        if bSimplify:
            simple_fermicorrected_offset_displaced = brownout + huygout # this is composed just of the unmatched amps
            simpledisplaced = np.copy(displaced)
            sumdisplaced = np.sum(displaced)

            forloopordering = np.arange(N).astype("int")
            rn.shuffle(forloopordering) 

            for iip, ip in enumerate(forloopordering):

                # if we make additional displacements (so as to make sumremainder==0) we must increase the amount
                # already displaced in one of the previously displaced arcs -- i.e. increase magnitude of the amount displaced --
                # as opposed to REDUCING the amount displaced in that arc, because the latter action would wind up undoing the
                # displacement that was required in order to allow the PEP to remain unviolated 

                bManagedToOffloadRemainder = False
                if simple_fermicorrected_offset_displaced[ip] * sumdisplaced <= 0:
                    simpledisplaced[ip] -= sumdisplaced
                    bManagedToOffloadRemainder = True
                    break

            rectified_simple_fermicorrected_offset_displaced = NetOut(simple_fermicorrected_offset_displaced)
            extra_displaced = simple_fermicorrected_offset_displaced - rectified_simple_fermicorrected_offset_displaced
            #return rectified_simple_fermicorrected_offset_displaced, extra_displaced

        # now you have PEP compatible output amplitudes (assuming you inject simple_fermicorrected_offset_displaced into heat-bath)


        if tuple(inputarcarr) == (1 , 0 ,-5, -1):
            import pdb; pdb.set_trace()       



        # note that brownout and huygout were calculated from unmatched, so that what follows keep track of everything.
        matchedplusbrown, matchedplusbrowndisplaced = FermiAdd(matched, displaced, brownout)
        fermicorrectedout, fermicorrecteddisplaced = FermiAdd(matchedplusbrown, matchedplusbrowndisplaced, huygout)

            


        # now add particles of opposite sign to what is already in fermicorrecteddisplaced (so as to make the sum of  fermicorrecteddisplaced zero)
        sumremainder = np.sum(fermicorrecteddisplaced)
        fermicorrected_offset_displaced = copy(fermicorrecteddisplaced)
        if sumremainder != 0:            
            # we must search the arcs randomly
            
            forloopordering = np.arange(len(fermicorrected_offset_displaced)).astype("int")
            rn.shuffle(forloopordering) 

            #if tuple(inputarcarr) == (5,  0, -7, -5):
            #    print("yy")
            #    import pdb; pdb.set_trace()
            
            for iip, ip in enumerate(forloopordering):

                # if we make additional displacements (so as to make sumremainder==0) we must increase the amount
                # already displaced in one of the previously displaced arcs -- i.e. increase magnitude of the amount displaced --
                # as opposed to REDUCING the amount displaced in that arc, because the latter action would wind up undoing the
                # displacement that was required in order to allow the PEP to remain unviolated 

                bManagedToOffloadRemainder = False
                if fermicorrected_offset_displaced[ip] * sumremainder <= 0:
                    fermicorrected_offset_displaced[ip] -= sumremainder
                    bManagedToOffloadRemainder = True
                    break


            if not(bManagedToOffloadRemainder):
                print("TROUBLE! must have at least one non-displaced arc for this to work")
                import pdb; pdb.set_trace()                   
            # now, having added in the offsetting displacements, the input arcs have been rendered Fermi-compatible

            # next, continue to rectify purely so as to reduce total number. Assume anything that goes into heat-bath will
            # get annihilated, so maximize for that
            #import pdb; pdb.set_trace()

            if tuple(inputarcarr) == (1 , 0 ,-5, -1):
                import pdb; pdb.set_trace()



            rectified_out = NetOut(fermicorrectedout)
            xtradisplaced = fermicorrectedout - rectified_out
            return rectified_out, fermicorrecteddisplaced + xtradisplaced

        else:
            # in this case the displacement to allow for PEP was already  balanced in terms of particles/antiparticles, so we just need to proceed to 
            # minimizing the total number of particles,, i.e.
            # we must do another rectify operation purely so as to reduce total number. Assume anything that goes into heat-bath will
            # get annihilated, so maximize for that

            if tuple(inputarcarr) == (1 , 0 ,-5, -1):
                import pdb; pdb.set_trace()


            rectified_out = NetOut(fermicorrectedout)
            xtradisplaced = fermicorrectedout - rectified_out
            
            return rectified_out, xtradisplaced + fermicorrected_offset_displaced

        """
        # next, continue to rectify purely so as to reduce total number. Assume anything that goes into heat-bath will
        # get annihilated, so maximize for that
        
        myressofar = huygout + brownout
        newdisplaced = np.zeros(displaced.shape).astype("int")
        if bBrownianRegularization:            
            if np.max(myressofar) * np.min(myressofar) < 0:
                finalres = NetOut(myressofar)
                newdisplaced = finalres - myressofar
                myressofar = finalres            

        bShuffle2D = False # we'll do this only to the final result when we update grid and grid_br
        if bShuffle2D:
            myressofar = ShuffleInToOut(myressofar)
            displaced = ShuffleInToOut(displaced)
            newdisplaced = ShuffleInToOut(newdisplaced)
            
        # 1  myressofar: PEP-satisfying output (and if bOverDrive is true, it is redndered as)
        # 2  displaced: the heat-bath adjustment needed to preserve PEP
        # 3  newdisplaced: the heat-bath adjustment resulting from 1)
        return np.array(myressofar).astype("int"), np.array(displaced).astype("int"), np.array(newdisplaced).astype("int")

        """


    def FermiOut2(inputarcarr, permut):
        """
        return the rectified wave, the displacements (which need not sum up to zero)
        and the displacements plus offsets (i.e. the heat-bath brownian pair injection), which always sums up to zero
        since it is particle/antiparticle pairs
        """



        N = len(inputarcarr)
        

        displaced = np.zeros(inputarcarr.shape).astype("int")
        # displaced will mean shunted into heath-bath mode because they convlict with PEP rules. Therefore, they must be offset
        # it is the sum of two shuntings (one that avoids violating PEP, and one to reduce the total number of particles)
        brownout = copy(inputarcarr)
        huygout = -inputarcarr # note that for rectangular lattices you will need an extra shuffle
        # to orient the arcs properly (e.g. the positive x direction for incoming is 0, but for outgoing the positive x direction is 2)
        # but we will do that right before we update the grid with what we're calculating here
        for i, ip in zip(np.arange(len(permut)), permut):
            # note that for any i such that i = permut[i], the amplitudes are automatically opposite sign, so that's fine
            brownout[ip] = inputarcarr[i] * 2
            if huygout[ip] * brownout[ip] > 0:
                displaced[ip] = -huygout[ip]
                huygout[ip] = 0



        fermicorrected_out = brownout + huygout
        if tuple(inputarcarr) == (1 , 0 ,-5, -1):
            import pdb; pdb.set_trace()
        if np.sum(inputarcarr) == 0 and np.sum(np.abs(inputarcarr)) != 0:
            import pdb; pdb.set_trace()
     

        # now add particles of opposite sign to what is already in displaced (so as to make the sum of  fermicorrected_offset_displaced zero)
        sumremainder = np.sum(displaced)
        fermicorrected_offset_displaced = copy(displaced)
        if sumremainder != 0:            
            # we must search the arcs randomly
            
            forloopordering = np.arange(len(fermicorrected_offset_displaced)).astype("int")
            rn.shuffle(forloopordering) 

            #if tuple(inputarcarr) == (5,  0, -7, -5):
            #    print("yy")
            #    import pdb; pdb.set_trace()
            
            for iip, ip in enumerate(forloopordering):

                # if we make additional displacements (so as to make sumremainder==0) we must increase the amount
                # already displaced in one of the previously displaced arcs -- i.e. increase magnitude of the amount displaced --
                # as opposed to REDUCING the amount displaced in that arc, because the latter action would wind up undoing the
                # displacement that was required in order to allow the PEP to remain unviolated 

                bManagedToOffloadRemainder = False
                if fermicorrected_offset_displaced[ip] * sumremainder <= 0:
                    fermicorrected_offset_displaced[ip] -= sumremainder
                    bManagedToOffloadRemainder = True
                    break


            if not(bManagedToOffloadRemainder):
                print("TROUBLE! must have at least one non-displaced arc for this to work")
                import pdb; pdb.set_trace()                   
            # now, having added in the offsetting displacements, the input arcs have been rendered Fermi-compatible

            # next, continue to rectify purely so as to reduce total number. Assume anything that goes into heat-bath will
            # get annihilated, so maximize for that
            #import pdb; pdb.set_trace()

            if tuple(inputarcarr) == (1 , 0 ,-5, -1):
                import pdb; pdb.set_trace()


        rectified_out = NetOut(fermicorrected_out)
        xtradisplaced = fermicorrected_out - rectified_out
        return rectified_out, fermicorrected_offset_displaced + xtradisplaced


        """
        # next, continue to rectify purely so as to reduce total number. Assume anything that goes into heat-bath will
        # get annihilated, so maximize for that
        
        myressofar = huygout + brownout
        newdisplaced = np.zeros(displaced.shape).astype("int")
        if bBrownianRegularization:            
            if np.max(myressofar) * np.min(myressofar) < 0:
                finalres = NetOut(myressofar)
                newdisplaced = finalres - myressofar
                myressofar = finalres            

        bShuffle2D = False # we'll do this only to the final result when we update grid and grid_br
        if bShuffle2D:
            myressofar = ShuffleInToOut(myressofar)
            displaced = ShuffleInToOut(displaced)
            newdisplaced = ShuffleInToOut(newdisplaced)
            
        # 1  myressofar: PEP-satisfying output (and if bOverDrive is true, it is redndered as)
        # 2  displaced: the heat-bath adjustment needed to preserve PEP
        # 3  newdisplaced: the heat-bath adjustment resulting from 1)
        return np.array(myressofar).astype("int"), np.array(displaced).astype("int"), np.array(newdisplaced).astype("int")

        """



    def FermiOut3(inputarcarr, permut):
        """
        In this version we have three contributions to the output amplitude as far as making the outputs PEP compliant
         (and then additional stuff if we want to rectify furhter):
        1) the "matched" incoming amplitudes (which by definition add up to zero) in which particles just bounce back in reverse direction
        2) the brownian particles (obtained by permuting the arcs of the input particles and multiplying by 2)
        3) the huygens particles (obtained by negating the input particles)

        (We won't call ShuffleInToOut() until just before we update grid[], so that we only have to do that once)

        Therefore, for any arc, we have 3 components. If some of them are different signs, then we can consolidate (and if the result leaves 2 consituents
        of differing signs, we can consolidate further); whatever is left over gets put into the "displaced" vector. When we're done adding those three up, we
        calculate sumdisplaced and then find a component of displaced that is either zero or the opposite sign of displaced (so that its absolute value INCREASES
        when sum_displaced is subtracted from it)
        """

        N = len(inputarcarr)

        unmatched = np.array(NetOut(inputarcarr)).astype("int")
        matched = inputarcarr - unmatched

        # here are out 3 constituents
        matchedout = -matched # these will simply "bounce" in reverse
        brownout = 2 * unmatched
        rn.shuffle(brownout)
        huygout = -unmatched # just as with 
        
        # and also the displaced, which holds any particles whose presence would mess up PEP
        displaced =  np.zeros(inputarcarr.shape).astype("int")
        outputarcarr = np.zeros(inputarcarr.shape).astype("int")

        for i in range(N):
            mo = matchedout[i]
            bo = brownout[i]
            ho = huygout[i]

            A = LumpyNetOut([mo, bo, ho])
            ordering = np.arange(len(A))
            rn.shuffle(ordering)


            # now, get PEP-compatible sum or A (along with displaced)
            thissign = 0
            idisplaced = 0
            sumiA = 0
            for j in ordering:
                if thissign == 0:
                    if A[j] > 0:
                        thissign = 1
                        sumiA += A[j]
                    elif A[j] < 0:
                        thissign = -1
                        sumiA += A[j]
                    else:
                        pass # A[j]==0 case
                else:
                    if np.sign(A[j]) == thissign: # if both are in same sign (state), then one must be shoved off into "displaced" category
                        idisplaced += A[j]
                    else:
                        sumiA += A[j]
            displaced[i] = idisplaced
            outputarcarr[i] = sumiA

        sumdisplaced = np.sum(displaced) 
        
        # as noted elsewhere, we have to adjust displaced so as to INCREASE the absolute value of whatever component is adjusted
        if sumdisplaced != 0:
            sgnsumdisplaced =  np.sign(sumdisplaced)
            bManagedToOffset = False
            ordering = np.arange(N)
            rn.shuffle(ordering)
            for i in ordering:
                if not(sgnsumdisplaced == -np.sign(outputarcarr[i])):

                    #if np.abs(displaced[i] - sumdisplaced) < np.abs(displaced[i]):
                    #    print("ERROR -- the absolute val of displaced[i] must INCREASE; get rid of this check once you see it works right")
                    outputarcarr[i] -= -sumdisplaced
                    displaced[i] += -sumdisplaced
                    bManagedToOffset = True
                    break
            if not(bManagedToOffset):
                print("DID NOT OFFSET -- HOW IS THIS POSSIBLE")
                import pdb; pdb.set_trace()
        
        # now, we still have to rectify the output 
        rectifiedoutputarcarr = NetOut(outputarcarr)
        xtradisplaced = outputarcarr - rectifiedoutputarcarr




        return rectifiedoutputarcarr, displaced + xtradisplaced
            


        """
        # next, continue to rectify purely so as to reduce total number. Assume anything that goes into heat-bath will
        # get annihilated, so maximize for that
        
        myressofar = huygout + brownout
        newdisplaced = np.zeros(displaced.shape).astype("int")
        if bBrownianRegularization:            
            if np.max(myressofar) * np.min(myressofar) < 0:
                finalres = NetOut(myressofar)
                newdisplaced = finalres - myressofar
                myressofar = finalres            

        bShuffle2D = False # we'll do this only to the final result when we update grid and grid_br
        if bShuffle2D:
            myressofar = ShuffleInToOut(myressofar)
            displaced = ShuffleInToOut(displaced)
            newdisplaced = ShuffleInToOut(newdisplaced)
            
        # 1  myressofar: PEP-satisfying output (and if bOverDrive is true, it is redndered as)
        # 2  displaced: the heat-bath adjustment needed to preserve PEP
        # 3  newdisplaced: the heat-bath adjustment resulting from 1)
        return np.array(myressofar).astype("int"), np.array(displaced).astype("int"), np.array(newdisplaced).astype("int")

        """


    def PosPart(xarr):
        myret = copy(xarr)
        for iic, ic in enumerate(myret):
            if ic < 0:
                myret[iic] = 0
        return myret
    
    def NegPart(xarr):
        myret = copy(xarr)
        for iic, ic in enumerate(myret):
            if ic > 0:
                myret[iic] = 0
        return myret
    
    def ReturnOneSign(xarr, PosNumberForPosPart_NegNumberForNegPart):
        if PosNumberForPosPart_NegNumberForNegPart > 0:
            return PosPart(xarr)
        if PosNumberForPosPart_NegNumberForNegPart < 0:
            return NegPart(xarr)
        print("Error in ReturnOneSign(); second arg must be either positive definite or negative definite; no zeros or indeterminates allowed.")
    
    def SplitPosAndNeg(xarr):
        myretpos = copy(xarr)
        myretneg = copy(xarr)
        for iic, ic in enumerate(xarr):
            if ic > 0:
                myretneg[iic] = 0
            elif ic < 0:
                myretpos[iic] = 0
            
        return myretpos, myretneg

    def UnloadToSavePEP(xarr, M):
        """
        Dump M particles into xarr provided that the arcs in which you dump are zero or of opposite sign.
        """

        if tuple(xarr) == (0, 0, 2, 2):
            import pdb; pdb.set_trace()
        myret = copy(xarr)
        indx = []
        amnt = []


        if np.sum(np.abs(amnt)) > 0:
            import pdb; pdb.set_trace()
        for iix,ix in enumerate(xarr):
            if ix * M <= 0:
                indx.append(iix)
                amnt.append(ix)
        
        amnt[0] -= M
        # now, smooth everything out

        amt = ModPoolJustBrownian(amnt)

        for i, iamt in zip(indx, amt):
            myret[i] += iamt

        return myret


    def Distribute(sumdisplaced, Narcs):
        """
        The displaced vector is one input, while the other is a binary (0/1) vector whose 1's indicate
        which arcs are able to receive contributions. Since these will all go into heat bath, we can
        just distribute the total amount displaced (sum(displaced)) among
        """

        if sumdisplaced == 0:
            return np.zeros((Narcs,)).astype("int")

        aa = [0] * (Narcs//2)
        aa[0] = sumdisplaced
        aa = ModPoolJustBrownian(aa)
        negaa = [-i for i in aa]
        aa = list(aa) + negaa
        rn.shuffle(list(aa))
        return np.array(aa).astype("int")




    def FermiOut6(inputarcarr, bModPoolTheIncoming=True):
        """
        In this version we will calculate the modpool output and then calc a deviation-from-PEP vector, and then
        a displacement vector (whose sum is zero) to make the deviation-from-PEP zero.
        """


        """
        def CapAbs1(xarr, cap=1):
            return( [ int(np.sgn(i) * np.min([np.abs(i), cap])) * np.sign(i)) for i in xarr] )

        rawout = ModPool(inutarcarr)
https://old.reddit.com/r/news/comments/1n1iiyo/minnesota_governor_says_a_shooting_has_occurred/
        capout = CapAbs1(rawout)

        displaced = rawout - capout

        """

        def WhichHaveMostOppositeComp(err, modpoolout, signsumerr):

            #import pdb; pdb.set_trace()
            capacity = np.zeros((len(err),)).astype("int")
            for i in range(len(err)):
                if err[i] * signsumerr <= 0:
                    capacity[i] = (modpoolout[i] * signsumerr - 1)


            minlevel = np.min(capacity)
            if minlevel == 0:
                return None
            if minlevel > 0:
                print("no more offloading possible:", inputarcarr)
                import pdb; pdb.set_trace()
            whichhavemin = []
            for jjerr, jerr in enumerate(capacity):
                if jerr == minlevel:
                    whichhavemin.append(jjerr)
            #import pdb; pdb.set_trace()
            return np.random.choice(whichhavemin)


        N = len(inputarcarr)

        if bModPoolTheIncoming:
            modpoolout = ModPool(inputarcarr)
        else:
            modpoolout = copy(inputarcarr)
        origmodpoolout = np.copy(modpoolout)
        if np.max(np.abs(modpoolout)) <= 1:
            return modpoolout,  np.zeros((N, )).astype("int")

        err = np.zeros((N, )).astype("int")
        for iim, im in  enumerate(modpoolout):
            if np.abs(im) > 1:
                err[iim] = im - np.sign(im)
                modpoolout[iim] = np.sign(im)



        sumerr = np.sum(err)
        signsumerr = int(np.sign(sumerr))

        if sumerr == 0:
            return modpoolout, err
        
        #if tuple(origmodpoolout) == (-1,  1,  2, -1):
        #    print("Ya")
        #    import pdb; pdb.set_trace()

        if False: # sumerr != 0:
        
            ordering = [i for i in range(N)]
            rn.shuffle(ordering)

            for i in ordering:
                if err[i] == 0 and np.sign(modpoolout[i]) != signsumerr:
                    if modpoolout[i] == 0 or np.abs(sumerr) == 1:
                        modpoolout[i] += signsumerr
                        err[i] -= signsumerr
                        sumerr -= signsumerr
                    elif modpoolout[i] == -signsumerr:
                        modpoolout[i] += signsumerr
                        err[i] -= signsumerr
                        sumerr -= signsumerr
                    else:
                        import pdb; pdb.set_trace()
                        modpoolout[i] = 2 * signsumerr
                        err[i] -= 2 * signsumerr
                        sumerr -= 2 * signsumerr
                    if sumerr == 0:
                        break     


        # we're just going to "drip" the error away, always filling the indices with the most capacity (i.e. whose values have the opposite sign of sumerr)
        #import pdb; pdb.set_trace()
        for i in range(np.abs(sumerr)):
            whichind = WhichHaveMostOppositeComp(err, modpoolout, signsumerr)
            if whichind is None:
                print("can't lay it all off AA")
                import pdb; pdb.set_trace()
            err[whichind] -= signsumerr
            modpoolout[whichind] += signsumerr

        if np.max(np.abs(modpoolout)) > 1:
            print("couldn't lay it all off")
            import pdb; pdb.set_trace()

        return modpoolout, err







    # not used!! -- discard
    def Unscramble(somethingshuffled, permut):
        N = len(somethingshuffled)
        displaced = np.zeros((N,)).astype("int")

        for iic, ic in enumerate(permut):
            displaced[ic] -= somethingshuffled[ic]
            displaced[iic] += somethingshuffled[permut[iic]]
        
        return displaced




    # not used!! -- discard
    def Unscramble(somethingshuffled, permut):
        N = len(somethingshuffled)
        displaced = np.zeros((N,)).astype("int")

        for iic, ic in enumerate(permut):
            displaced[ic] -= somethingshuffled[ic]
            displaced[iic] += somethingshuffled[permut[iic]]
        
        return displaced


    def FermiOut8(inputarcarr):
        """
        In this version we will calculate a PEP-compatible output (which, alas, blows up big-time), then
        calculate the heat-bath displacement needed to shift everything in this exploding PEP compatible output to have it instead
        be the same modpool output that we crunch for boson particles (which we know already do not blow up)
        """


        """
        def CapAbs1(xarr, cap=1):
            return( [ int(np.sgn(i) * np.min([np.abs(i), cap])) * np.sign(i)) for i in xarr] )

        rawout = ModPool(inutarcarr)

        capout = CapAbs1(rawout)

        displaced = rawout - capout

        """

        N = len(inputarcarr)

        testmodpool = ModPool(inputarcarr)


        #if tuple(inputarcarr) == (1,  0, -1, -1):
        #    import pdb; pdb.set_trace()
        #    import pdb; pdb.set_trace()



        #if np.sum(np.abs(inputarcarr)) > 3:
        #    print("bla")
        #    import pdb; pdb.set_trace()

        # we can net out the Brownian particles without violating PEP so we do that first
        brownpool = np.array(NetOutStrict(2*inputarcarr)).astype("int")
        sumbrownpool = np.sum(brownpool)

        shufflebrownpool = copy(brownpool)
        rn.shuffle(shufflebrownpool)
        shufflebrownpool = shufflebrownpool 

        displaced = shufflebrownpool - brownpool

        # now you have satisfied PEP (since, after the unscrambling, the brownian contributions will be aligned with the opposite-sign huygens contributions)

        huygout = -np.array(inputarcarr).astype("int")

        # check for compliance

        for ic_h, ic_b in zip(huygout, brownpool):
            if ic_h * ic_b > 0:
                print("SOMETHING WRONG -- not PEP compliant")
                import pdb; pdb.set_trace()

        pepcompliantout = huygout + brownpool # again, note that this is brownpool and NOT shuffledbrownpool

        bNetOutAgain = False
        if bNetOutAgain:
            nettedpepcompliant = NetOut(pepcompliantout)
            xtradisplaced = nettedpepcompliant - pepcompliantout

            return pepcompliantout, -(displaced + xtradisplaced) # we negate it because displaced is negated before being added to heat bath (for earlier reasons)
        else:
            return pepcompliantout, -displaced # we negate it because displaced is negated before being added to heat bath (for earlier reasons)









    def BoseOut(inputarcarr):

        unmatched = np.array(NetOut(inputarcarr)).astype("int")
        matched = inputarcarr - unmatched
        #if np.sum(np.abs(matched)) > 3:
        #    import pdb; pdb.set_trace()
        outputarcarr = matched + ModPool(unmatched)

        rectifiedoutputarcarr = NetOut(outputarcarr)
        displaced = outputarcarr - rectifiedoutputarcarr

        return rectifiedoutputarcarr, displaced
            




    def NetOutWOverdrive(inlist):
        """
        INCOMPLETE--DO NOT USE
        This differs from NetOut in that it goes beyond shaving until Brownian particles of only
        one sign remain. Rather, it keeps shaving until the remaining signed Brownian particles
        are as smoothly distributed as can be (even though that latter shifting around does not
        recude the total number of particles)

        """
        print("DO NOT USE -- DEPRECATED -- JUST BOOST THE ENTIRE SET OF ARCS AND USE NetOut()")
        N = len(inlist)
        posportion = 0
        negportion = 0
        poslist = []
        neglist = []
        suminlist = np.sum(inlist)
        sgninlist = int(np.sign(suminlist))

        posdivinlist = int(np.ceil(suminlist/float(N)))
        negdivinlist = int(np.floor(suminlist/float(N)))

        

        divinlist = sgninlist * (np.abs(suminlist) // len(inlist))
        modinlist = sgninlist * (np.abs(suminlist) % len(inlist))


        for ii, iflow in enumerate(inlist):
            if iflow > divinlist:
                poslist.append(iflow)
                neglist.append(0)
                posportion += iflow
            elif iflow < divinlist:
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




    bBrownianRegularization = True


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

                

                
                    
                inspan = max(grid[iscen,ix,iy,:]) - min(grid[iscen,ix,iy,:]) # not needed here
                
                incoming_along_xpos = copy(grid[iscen, ix-1, iy, 0])  # note the origin of this flow is in the xneg direction, i.e. ix-1
                incoming_along_ypos = copy(grid[iscen, ix, iy-1, 1])
                incoming_along_xneg = copy(grid[iscen, (ix+1) % grid.shape[1], iy, 2])
                incoming_along_yneg = copy(grid[iscen, ix, (iy+1) % grid.shape[2], 3])


                
                inparticlecount = np.abs(incoming_along_xpos) + np.abs(incoming_along_ypos) + np.abs(incoming_along_xneg) + np.abs(incoming_along_yneg)
                incarr = np.array([incoming_along_xpos, incoming_along_ypos, incoming_along_xneg, incoming_along_yneg]).astype("int")



                incoming_along_br_xpos = copy(grid_br[iscen, ix-1, iy, 0])  # note the origin of this flow is in the xneg direction, i.e. ix-1
                incoming_along_br_ypos = copy(grid_br[iscen, ix, iy-1, 1])
                incoming_along_br_xneg = copy(grid_br[iscen, (ix+1) % grid_br.shape[1], iy, 2])
                incoming_along_br_yneg = copy(grid_br[iscen, ix, (iy+1) % grid_br.shape[2], 3])

                inparticlecount_br = np.abs(incoming_along_br_xpos) + np.abs(incoming_along_br_ypos) + np.abs(incoming_along_br_xneg) + np.abs(incoming_along_br_yneg)
                incarr_br = np.array([incoming_along_br_xpos, incoming_along_br_ypos, incoming_along_br_xneg, incoming_along_br_yneg]).astype("int")


                #    print("in", np.sum(incarr), incarr, incarr_br)
                #    import pdb; pdb.set_trace()
                #    bPrintOutFlow = True

                # we have three things to do
                # 1) FermiOut(incarr) -- note this is on INCOMING
                # 2) consolidate any previous heat-bath particles (i.e. annihilate particle-antiparticles so that only one sign remains) 
                #  via the incoming of the surviving heat-bath  given by NetOut(incarr_br) and the outgoing
                # 3) supplement the remaining heat-bath particles with those produced by FermiOut (don't worry about PauliExPrinc since these are unobservable)
                #    AND GENERATE PARTICLE ANTIPARTICLE PAIRS IN THE MAIN GRID

                #outarr, displaced, newdisplaced = FermiOut(incarr)


                if inparticlecount > 2:
                    incoming_along_xpos_br = copy(grid_br[iscen, ix-1, iy, 0])  # note the origin of this flow is in the xneg direction, i.e. ix-1
                    incoming_along_ypos_br = copy(grid_br[iscen, ix, iy-1, 1])
                    incoming_along_xneg_br = copy(grid_br[iscen, (ix+1) % grid.shape[1], iy, 2])
                    incoming_along_yneg_br = copy(grid_br[iscen, ix, (iy+1) % grid.shape[2], 3])

                    incarr_br = np.array([incoming_along_xpos_br, incoming_along_ypos_br, incoming_along_xneg_br, incoming_along_yneg_br]).astype("int")
                    #if np.sum(np.abs(incarr_br)) > 2:
                    #    print("busy")      
                    #    import pdb; pdb.set_trace()
                        #brownout, displaced, leftover = FermiOut(incarr, thisperm)
                        #import pdb; pdb.set_trace()

                #if tuple(incarr) == (1,  0, -1, -1):
                #    import pdb; pdb.set_trace()

                #if tuple(incarr) == (19937,  -528,  -894, -5122):
                #    import pdb; pdb.set_trace()
                #if np.sum(np.abs(incarr)) > 100 and t > 10:
                #    import pdb; pdb.set_trace()


                bFermi = True
                if bFermi:
                    #thisperm = np.arange(grid.shape[3])
                    #rn.shuffle(thisperm)

                    if tuple(incarr) == ( 0, 0, 0, 1 ) and ix == 0 and iy == 2 and t == 1:
                        print("woo")
                        import pdb; pdb.set_trace()

                    #if t % 10 == 1:
                    #    import pdb; pdb.set_trace()

                    #if np.sum(np.abs(incarr)) > 0 and np.sum(incarr) == 0:
                    #    import pdb; pdb.set_trace()



                    # note the output here has been adjustd both to allow PEP to preserved, and then again to reduce overall number  of particles
                    fermibrownhuygout, displaced = FermiOut6(incarr)

                    #import pdb; pdb.set_trace()
                    #                        wavestuff        + forward/samesign heatbathspew     - reverse/oppositesignheatbathspew
                    #if np.sum(np.abs(incarr)) > 0 and np.max(np.abs(incarr)) == 1 and t > 12:
                    #    import pdb; pdb.set_trace()


                    #netted_heat_bath_predisplaced = NetOut(incarr_br) # 


                    # the reverse-propagation contribution cannot be netted out in any way -- assuming it came into the node, 
                    #then a reverse propagating outward contribution has to be emitted
                    # without any randomization
                    #reversed_heat_bath = ShuffleInToOut( np.array(incarr_br).astype("int") )

                    forward_emission_heat_bath_contribution = ModPoolJustBrownian(incarr_br + displaced)

                    reversed_heat_bath_contribution = -np.array(incarr_br).astype("int") 

                    netted_heat_bath_displaced = forward_emission_heat_bath_contribution + reversed_heat_bath_contribution

                    outgoing_inc_heat_bath = ShuffleInToOut(fermibrownhuygout + netted_heat_bath_displaced)


                    xtradisplaced = np.zeros((len(incarr),)).astype("int")
                    if np.max(np.abs(outgoing_inc_heat_bath)) > 1:
                        orig_outgoing_inc_heat_bath = copy(outgoing_inc_heat_bath)
                        #print("trouble ahead")
                        #import pdb; pdb.set_trace()
                        outgoing_inc_heat_bath, xtradisplaced = FermiOut6(outgoing_inc_heat_bath, False)
                    grid[iscen, ix, iy, :] = outgoing_inc_heat_bath
                    grid_br[iscen, ix, iy, :] = ShuffleInToOut(forward_emission_heat_bath_contribution + xtradisplaced)
                
                else:
                    bosebrownhuygout, displaced = BoseOut(incarr)

                    
                    netted_heat_bath = ModPoolJustBrownian(incarr_br) - displaced
                    # why do we subtract displaced? Because the amplitudes in displaced are what needs to be canceled out, so we add the negative of that to the heat bath

                    #for the heat-bath propagation, we have to deterministically spew (into grid) an opposite sign and opposite direction particle for every 
                    # existing particle at this arc of grid_br, and also a same sign particle in some random direction 
                    # (don't worry about PEP, because these are fundamentall unobservable)
                    
                    
                    # so here spew out the opposite-direction-opposite-sign versions -- this will wind up in grid, not grid_br -- note it will be negated
                    # when we include it into the heat bath
                    reverse_netted_heat_bath = np.array(netted_heat_bath).astype("int")
                    
                    # now do the same-sign ejection
                    randomized_netted_heat_bath = copy(netted_heat_bath)
                    rn.shuffle(randomized_netted_heat_bath) 

                    # arguably, we should run this through ShuffleInToOut(), but since we've randomized it already, that is unnucessary


                    outgoing_inc_heat_bath = ShuffleInToOut(bosebrownhuygout + randomized_netted_heat_bath - reverse_netted_heat_bath)

                    grid_br[iscen, ix, iy, :] = ShuffleInToOut(randomized_netted_heat_bath)


                    grid[iscen, ix, iy, :] = outgoing_inc_heat_bath


                outspan = max(grid[iscen,ix,iy,:]) - min(grid[iscen,ix,iy,:])



                #outparticlecount = abs(outgoing_along_xpos) + abs(outgoing_along_ypos) + abs(outgoing_along_xneg) + abs(outgoing_along_yneg)
                outparticlecount = np.sum(np.abs(outgoing_inc_heat_bath))

                if np.sum(incarr) != np.sum(grid[iscen,ix,iy,:]):
                    print("mismatch", incarr, grid_br[iscen,ix,iy,:])
                    import pdb; pdb.set_trace()

                #if(outparticlecount != 0):
                #    print(outparticlecount)

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
                            aa = copy(a)
                            if (absa + absb + absc + absd) <= (np.abs(a - aa) + np.abs(b + aa) + np.abs(c - aa) + np.abs(d + aa)):
                                #print("wouldn't change: ", a, b, c, d, " start ", (absa + absb + absc + absd), "end", (np.abs(a - aa) + np.abs(b + aa) + np.abs(c - aa) + np.abs(d + aa)))
                                #import pdb; pdb.set_trace()
                                continue


                            grid[iscen, ix, iy, 0] -= aa
                            grid[iscen, ix, iy, 1] += aa
                            grid[iscen, ixp1, iyp1, 2] -= aa
                            grid[iscen, ixp1, iyp1, 3] += aa
                        elif ind == 1:
                            bb = copy(b)
                            if (absa + absb + absc + absd) <= (np.abs(a - bb) + np.abs(b + bb) + np.abs(c - bb) + np.abs(d + bb)):
                                #print("wouldn't change: ", a, b, c, d, " start ", (absa + absb + absc + absd), "end", (np.abs(a - bb) + np.abs(b + bb) + np.abs(c - bb) + np.abs(d + bb)))
                                #import pdb; pdb.set_trace()
                                continue
                            
                            grid[iscen, ix, iy, 0] -= bb
                            grid[iscen, ix, iy, 1] += bb
                            grid[iscen, ixp1, iyp1, 2] -= bb
                            grid[iscen, ixp1, iyp1, 3] += bb
                        elif ind == 2:
                            cc = copy(c)
                            if (absa + absb + absc + absd) <= (np.abs(a - cc) + np.abs(b + cc) + np.abs(c - cc) + np.abs(d + cc)):
                                #print("wouldn't change: ", a, b, c, d, " start ", (absa + absb + absc + absd), "end", (np.abs(a - cc) + np.abs(b + cc) + np.abs(c - cc) + np.abs(d + cc)))
                                #import pdb; pdb.set_trace()
                                continue
                                                    
                            grid[iscen, ix, iy, 0] -= cc
                            grid[iscen, ix, iy, 1] += cc
                            grid[iscen, ixp1, iyp1, 2] -= cc
                            grid[iscen, ixp1, iyp1, 3] += cc  
                        elif ind == 3:               
                            dd = copy(d)
                            if (absa + absb + absc + absd) <= (np.abs(a - dd) + np.abs(b + dd) + np.abs(c - dd) + np.abs(d + dd)):
                                #print("wouldn't change: ", a, b, c, d, " start ", (absa + absb + absc + absd), "end", (np.abs(a - dd) + np.abs(b + dd) + np.abs(c - dd) + np.abs(d + dd)))
                                #import pdb; pdb.set_trace()
                                continue
                                                                
                            grid[iscen, ix, iy, 0] -= dd
                            grid[iscen, ix, iy, 1] += dd
                            grid[iscen, ixp1, iyp1, 2] -= dd
                            grid[iscen, ixp1, iyp1, 3] += dd

                        if ind == 4: # a is minimum
                            aa = -copy(a)
                            if (absa + absb + absc + absd) <= (np.abs(a - aa) + np.abs(b + aa) + np.abs(c - aa) + np.abs(d + aa)):
                                #print("wouldn't change: ", a, b, c, d, " start ", (absa + absb + absc + absd), "end", (np.abs(a - aa) + np.abs(b + aa) + np.abs(c - aa) + np.abs(d + aa)))
                                #import pdb; pdb.set_trace()
                                continue


                            grid[iscen, ix, iy, 0] -= aa
                            grid[iscen, ix, iy, 1] += aa
                            grid[iscen, ixp1, iyp1, 2] -= aa
                            grid[iscen, ixp1, iyp1, 3] += aa
                        elif ind == 5:
                            bb = -copy(b)
                            if (absa + absb + absc + absd) <= (np.abs(a - bb) + np.abs(b + bb) + np.abs(c - bb) + np.abs(d + bb)):
                                #print("wouldn't change: ", a, b, c, d, " start ", (absa + absb + absc + absd), "end", (np.abs(a - bb) + np.abs(b + bb) + np.abs(c - bb) + np.abs(d + bb)))
                                #import pdb; pdb.set_trace()
                                continue
                            
                            grid[iscen, ix, iy, 0] -= bb
                            grid[iscen, ix, iy, 1] += bb
                            grid[iscen, ixp1, iyp1, 2] -= bb
                            grid[iscen, ixp1, iyp1, 3] += bb
                        elif ind == 6:
                            cc = -copy(c)
                            if (absa + absb + absc + absd) <= (np.abs(a - cc) + np.abs(b + cc) + np.abs(c - cc) + np.abs(d + cc)):
                                #print("wouldn't change: ", a, b, c, d, " start ", (absa + absb + absc + absd), "end", (np.abs(a - cc) + np.abs(b + cc) + np.abs(c - cc) + np.abs(d + cc)))
                                #import pdb; pdb.set_trace()
                                continue
                                                    
                            grid[iscen, ix, iy, 0] -= cc
                            grid[iscen, ix, iy, 1] += cc
                            grid[iscen, ixp1, iyp1, 2] -= cc
                            grid[iscen, ixp1, iyp1, 3] += cc  
                        elif ind == 7:               
                            dd = -copy(d)
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

                        
                        outgoing_along_xpos = copy(grid[iscen, ix, iy, 0])  
                        outgoing_along_ypos = copy(grid[iscen, ix, iy, 1])
                        outgoing_along_xneg = copy(grid[iscen, ix, iy, 2])
                        outgoing_along_yneg = copy(grid[iscen, ix, iy, 3])

                        outgoing_along_xneg = copy(grid[iscen, ix, iy, 2])
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

                                        


                    
                    #incoming_along_xpos = copy(grid[iscen, ix-1, iy, 0])  # note the origin of this flow is in the xneg direction, i.e. ix-1
                    #incoming_along_ypos = copy(grid[iscen, ix, iy-1, 1])
                    #incoming_along_xneg = copy(grid[iscen, (ix+1) % grid.shape[1], iy, 2])
                    #incoming_along_yneg = copy(grid[iscen, ix, (iy+1) % grid.shape[2], 3])






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

                    #if ix == 0 and iy == 3 and t == 0:
                    #    import pdb; pdb.set_trace()

                    origarcs = copy(grid[iscen, ix, iy, :])

                    minarcs = np.min(origarcs)
                    maxarcs = np.max(origarcs)

                    #import pdb; pdb.set_trace()
                   
                    #while minarcs * maxarcs < 0 and np.min([maxarcs, -minarcs]) > 1:


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


                    #minarcs = np.min(origarcs)
                    #maxarcs = np.max(origarcs)

                        





def ProcessFile():
  
  global bUseLatticeGasForRand
  global bCplusplus
  global modpoolvec
  global alpha
  
  opts, args = ReadParams()

  alpha = opts.alpha
  #import pdb; pdb.set_trace()
  myseed = abs( opts.seed )
  
  # for reasons beyond me, both these seeds have to be set or results will not replicate:
  # see https://stackoverflow.com/questions/46661426/why-random-seed-does-not-make-results-constant-in-python
  # and also https://stackoverflow.com/questions/31057197/should-i-use-random-seed-or-numpy-random-seed-to-control-random-number-gener
  # also, try using
  rn.seed(myseed)
  np.random.seed( myseed  )


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
      
      grid[i_scenario, i,j,Comp000] = copy(new000)
      grid[i_scenario, i,j,Comp090] = copy(new090)
      grid[i_scenario, i,j,Comp180] = copy(new180)
      grid[i_scenario, i,j,Comp270] = copy(new270)
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
          incomingbraketlist = list(copy(grid.shape))
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
          outgoingbraketlist = list(copy(grid.shape))
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
                    grid[:, 0, 2, 3] = 1
                    grid[:, 0, 2, 1] = -1
                    grid[:, 1, 1, 2] = 1
                    grid[:, 1, 3, 0] = -1
                    grid[:, 3, 1, 2] = 1
                    grid[:, 1, 1, 3] = -1
                    grid[:, 3, 1, 1] = 1
                    grid[:, 3, 1, 3] = -1
                    grid[:, 2, 0, 0] = 1
                    grid[:, 2, 0, 2] = -1
                    
                    grid[:, 4, 0, 0] = 1
                    grid[:, 4, 0, 2] = -1

                    
                  
                    """
                    grid_brown[:, 0, 0, 0] = 1 
                    grid_brown[:, 1, 1, 1] = -1
                    grid_brown[:, 0, 2, 3] = 1
                    grid_brown[:, 0, 2, 1] = -1
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
        minchunk.append(copy(minthischunk))
        maxchunk.append(copy(maxthischunk))
        minthischunk = 0
        maxthischunk = 0
    
    if mingrid < minthischunk:
        minthischunk = copy(mingrid)
    if maxgrid > maxthischunk:
        maxthischunk = copy(maxgrid)
        
       
        
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
            maxarr[-1] = copy(maxarr[-1] - maxarr[0])
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
    print2("sum     " + str(sumgrid) + " " + str(absgrid) + " fermipct %5.4f" % (np.sqrt(fermipct),)) # + " " + str(np.sum(grid[:,0::2,0::2,:])) + " " + str(np.sum(grid[:,0::2,1::2,:]))  )
    print2("sumXmom " + str(np.sum(grid[:,:,:,Comp000] - grid[:,:,:,Comp180])) + " sumYmom " + str(np.sum(grid[:,:,:,Comp090] - grid[:,:,:,Comp270])) +  ' max ' + str(maxgrid) + ' min ' + str(np.min(grid)) +  ' 90pct ' + str(np.round(np.percentile(grid, 90))) + ' 10pct ' + str(np.round(np.percentile(grid, 10))) + ' mnab ' + str(minabs[0]) + ' ' + str(minabs[1]) )
 
    print2("sumBR     " + str(sumgrid_brown) + " " + str(absgrid_brown) + " " + " %5.1f" % (sumgrid_brown2/gridct,)) # + " " + str(np.sum(grid[:,0::2,0::2,:])) + " " + str(np.sum(grid[:,0::2,1::2,:]))  )
    #print2("sumXmomBR " + str(np.sum(grid_brown[:,:,:,Comp000] - grid[:,:,:,Comp180])) + " sumYmom " + str(np.sum(grid[:,:,:,Comp090] - grid[:,:,:,Comp270])) +  ' max ' + str(maxgrid) + ' min ' + str(np.min(grid)) +  ' 90pct ' + str(np.round(np.percentile(grid, 90))) + ' 10pct ' + str(np.round(np.percentile(grid, 10))) + ' mnab ' + str(minabs[0]) + ' ' + str(minabs[1]) )
    if t % 1000 == 0:
        import pdb; pdb.set_trace()        
 
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
        rmsarr = []
        #ferminoic = [] # gives percentage of nodes where abs(amplitude) < 2
        tarr = []
        gridsnapshot = []
        tarr2 = []
        
    if t % 100 == 0 and t > 0:
        filename1 = "blah.csv"
        rmsarr.append(np.sqrt(sumamp2)) # (np.sqrt(np.abs(np.sum(grid[:,:,:,:]*grid[:,:,:,:]))))
        #gridct = float(np.prod(list(grid.shape)))/2.0
        #fermipct = ((np.sum(np.abs(grid) <= 1) - gridct)) / float(gridct)
        #ferminoic.append(fermipct)
        tarr.append(t)

        #emaarrnp = np.vstack([np.array(tarr), np.array(emamaxarr), rmsarr])

        

        #np.savetxt(filename1, emaarrnp.T, fmt='%5.1f', delimiter=',')
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

python regress_on_etf_spy_lag_correllimited_float_wOFFSET_lessnoetfnospy.py -f /Volumes/PNT3/polygon_data/%/qt_@_%.bz2 --tk F --nonspy XLY --startdt 2022-01-01 --enddt 2023-09-01 --regfile trash_per4_smallnospyfixed_regfileY_@_@@_220101to230901_%sec_lagord6_err_%%_4wndw.txt --freq 5 --lag 5 --threshold -0.15 > /dev/null
./BR_vwapREV.exe --threshold 0.00015 --submitthreshold 0.00015 --TEMPregressionfile /Users/hrvojehrgovcic/quant/polygon/trash_per4_smallnospyfixed_regfileY_F_XLY_220101to230901_5sec_lagord6_err_n0p15_4wndw.txt   --output /Users/hrvojehrgovcic/quant/polygon/trash/temp.txt -y 2023 --startdt 2023-01-04 --enddt 2023-05-26 --exclude /Users/hrvojehrgovcic/quant/data/HolidayExcludeDates.csv --symbol F  --etf XLY --file /Volumes/PNT3/polygon_data/qtfmt_%_@.csv --ticksize 1 --dnin 50 --upin 50 --upout 30450 --probability 2.0 --latency 0.005 --starttimehistory 01:30 --starttime 09:35 --endtime 15:58 --endtimedealinitiation 15:40 --endtimehistory 16:00  --meanreversionfilter +1 --allowgap --dlysum --maxlimitorder 1 --seed 500 --fee 0.14 --rebate 0.4 --append --specificlookback 0 --tracefreq 24hr --fillonlyontrade --mindeals 5 --insidedistance 0 --fastema 100 --slowema 50 --nocancelbecausequantities --vwapinterval 2 > blah_F_2sec_vwap_0p00015_per4fixed2.txt
python simplesr.py blah_F_5sec_vwap_0p00015_per4fixed2.txt
pmset sleepnow

./BR_vwapREV.exe --threshold 0.00015 --submitthreshold 0.00015 --TEMPregressionfile /Users/hrvojehrgovcic/quant/polygon/TRcoeff_nogapremovebadcodes_per4_smallnospyfixed_regfileY_F_XLY_220101to230901_2sec_lagord6_err_n0p15_4wndw.txt --output /Users/hrvojehrgovcic/quant/polygon/trash/temp.txt -y 2023 --startdt 2023-01-04 --enddt 2023-05-26 --exclude /Users/hrvojehrgovcic/quant/data/HolidayExcludeDates.csv --symbol F  --etf XLY --file /Volumes/PNT3/polygon_data/qtfmt_%_@.csv --ticksize 1 --dnin 50 --upin 50 --upout 30450 --probability 2.0 --latency 0.005 --starttimehistory 01:30 --starttime 09:35 --endtime 15:58 --endtimedealinitiation 15:40 --endtimehistory 16:00  --meanreversionfilter +1 --allowgap --dlysum --maxlimitorder 1 --seed 500 --fee 0.14 --rebate 0.4 --append --specificlookback 0 --tracefreq 24hr --fillonlyontrade --mindeals 5 --insidedistance 0 --fastema 100 --slowema 50 --nocancelbecausequantities --vwapinterval 2 > blah_F_2sec_vwap_0p00015_per4filteredcodes.txt
python simplesr.py  blah_F_2sec_vwap_0p00015_per4filteredcodes.txt
pmset sleepnow


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

#plt.plot(a, b)
#plt.show()

from scipy.stats import linregress
linregress(a, b)





# python  ~/quant/latticegas_youngsdoubleslitexp_partialpooling_inbetween_bernouilli_modpool_FERMIattmpt2.py   --discrete --seed 24 --picklefile trash.pkl --alpha 1.0


%run   latticegas_youngsdoubleslitexp_partialpooling_inbetween_bernouilli_modpool_FERMIattmpt2.py    --discrete --seed 24 --picklefile trash.pkl --alpha 1.0
%run  latticegas_youngsdoubleslitexp_partialpooling_inbetween_bernouilli_modpool_FERMIattmpt4.py   --discrete --seed 45 --picklefile trash3.pkl --alpha 1.0

"""