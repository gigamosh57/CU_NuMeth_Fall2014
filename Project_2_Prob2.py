#### Page Weil
#### 10/8/14
#### CVEN 5537
#### PROJECT 2 , PROBLEM 1

import numpy as np
# For pi
import math
# for plotting
import matplotlib.pyplot as plt

####################################################################
######## Runge-Kutta ODE solver with Cash-Karp 4th-5th order params

def pagerkck4(feval,x0,tstart,tfinal,dt,ord = 4,tol=10**-6):
  # INITIALIZE ARRAY OF CASH-KARP FACTORS
  ck_a = [[0],[1/5.],[3/10.],[3/5.],[1.],[7/8.]]
	ck_b = np.array([[0,0,0,0,0],
				  [1/5.,0,0,0,0],
					  [3/40.,9/40.,0,0,0],
					  [3/10.,-9/10.,6/5.,0,0],
					  [-11/54.,5/2.,-70/27.,35/27.,0],
					  [1631/55296.,175/512.,575/13824.,44275/110592.,253/4096.],
					  ])
  ck_c =  [[37/378.],[0.],[250/621.],[125/594.],[0.],[512/1771.]]
  ck_cs = [[2825/27648.],[0.],[18575/48384.],[13525/55296.],[277/14336.],[1/4.]]
  
  ck_a = ck_a[range(0,ord)]
  ck_b = ck_b[range(0,ord),][:,range(0,ord)]
  ck_c = ck_c[range(0,ord)]
  ck_cs = ck_cs[range(0,ord)]
  
  #initialize time vector
  tvec=[tstart]
  #initialize x, initialize solution matrix xsol
  x=x0
  xsol=[] + [x,]
  ts = []
  
  # SINCE THIS FUNCTION USES ADAPTIVE TIMESTEPPING
  # THE TIME VECTOR IS NOT DEFINED INITIALLY
  t = tstart
  for t < tfinal:
    # RESET INITIAL DT
    h = dt
    
	xs = []
    ## WHILE ERROR IS ABOVE THRESHOLD
    while error > tol:
      
      # EVALUATE STAGES
      k = []
      
      # LOOPS THROUGH ALL K VALUES DEPENDING ON ORDER OF FUNCTION
      for klev in range(0,ord):
          
        # CALCULATES F FOR ALL FUNCTIONS IN FEVAL
        f = []
        
        # EVALUATES F DEPENDING ON LEVEL
        for fn in feval:
			f.append(fn(x+ck_b[range(0,(klev-1))]*k,t+ck_a[klev]))
        k = list.append(h*f)
		
	  # CALCULATE x values
      xnew = x + ck_a*k
	  xs = x + ck_cs*k
	  
	  # ESTIMATE ERROR
	  error = max(np.abs((xnew-xs)/xs))
	  hnew = h*(error/tol)**(0.2)
	  # ADJUST DT
	  h = hnew
    
    # STORE SOLUTION (we assume XS is better)
	xsol = xsol + [xs]
	
	# ADJUST TIME BY H
	t = t + h 
	tvec = tvec + [t]
	ts = ts + [h]

    
  return tvec,xsol,ts
######## Solver ends here
################################################

