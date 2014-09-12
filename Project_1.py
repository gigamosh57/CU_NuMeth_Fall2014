#### PSEUDOCODE FOR FINDING IDEAL BOUNDARIES FOR ROOT FINDING

########### Startup code

import numpy as np
# For pi
import math
# for plotting
import matplotlib.pyplot as plt
# for fsolve
from scipy.optimize import fsolve

########### SOLVER HERE
########### MODIFY THIS TO USE NEWTON'S METHOD INSTEAD OF BISECTION

# Bisection solver begins here

def proot(bounds,n = 1000,tol = 0.01): 
  it = 0
  i = 0
  while (i < n):
      i = i + 1
      it = i
      xmid = np.mean(bounds)
      t = [f(xmid),f(bounds[0])] 
      if (all( item > 0 for item in t) or all( item <= 0 for item in t)):
            bounds[0] = xmid
      else: 
           t = [f(xmid),f(bounds[1])] 
           if (all( item > 0 for item in t) or all( item <= 0 for item in t)):
                bounds[1] = xmid
      if (np.absolute(f(xmid)) < tol):
            i = n+1
  return it, xmid, f(xmid)
  
  
########### Define Function and initial params
f = lambda x : x * np.tan(x) - r
bd = [0,1.57]
r = 3

########### Solve problem
for r in np.arange(0,10.5,0.5):
  print proot(bd)
