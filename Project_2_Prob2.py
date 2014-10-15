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
  
  #initialize time vector
  tvec=[tstart]
  #initialize x, initialize solution matrix xsol
  x=x0
  xsol=[] + [x,]
  
  # SINCE THIS FUNCTION USES ADAPTIVE TIMESTEPPING
  # THE TIME VECTOR IS NOT DEFINED INITIALLY
  t = tstart
  for t < tfinal:
    # RESET INITIAL DT
    h = dt
    
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
      k2 = h*f*(x + 
    
     # ESTIMATE ERROR
    
     # ADJUST DT
    
    # ADJUST TIME BY H
    t = t + h 
    
  
    #stage 1
    k1=f(t,x)*dt;
    #stage 2
    k2=f(t+dt/2,x+k1/2)*dt;
    #stage 3
    k3=f(t+dt/2,x+k2/2)*dt;
    #stage 4
    k4=f(t+dt,x+k3)*dt;
    #calculate x at next time-step (using same variable name x, I don't need xold, xnew etc.)
    x=x+(k1+2*k2+2*k3+k4)/6;
    #update time vector and xsol matrix
    tvec=tvec + [t+dt]
    xsol=xsol + [x]
  return tvec,xsol
######## RK4 ODE solver ends here
################################################
