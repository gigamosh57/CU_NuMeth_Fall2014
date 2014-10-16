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

def pagerkck4(feval,x0,tstart,tfinal,dt,order = 4,tol=10**-6):
   # INITIALIZE ARRAY OF CASH-KARP FACTORS
   ck_a = [[0],[1/5.],[3/10.],[3/5.],[1.],[7/8.]]
   ck_b = np.array([[0,0,0,0,0],
               [1/5.,0,0,0,0],
               [3/40.,9/40.,0,0,0],
               [3/10.,-9/10.,6/5.,0,0],
               [-11/54.,5/2.,-70/27.,35/27.,0],
               [1631/55296.,175/512.,575/13824.,44275/110592.,253/4096.]])
   ck_c =  [[37/378.],[0.],[250/621.],[125/594.],[0.],[512/1771.]]
   ck_cs = [[2825/27648.],[0.],[18575/48384.],[13525/55296.],[277/14336.],[1/4.]]
  
   ck_a = ck_a[range(0,order)]
   ck_b = ck_b[range(0,order),][:,range(0,order)]
   ck_c = ck_c[range(0,order)]
   ck_cs = ck_cs[range(0,order)]
  
   #initialize time vector
   tvec=[tstart]
   #initialize x, initialize solution matrix xsol
   x=x0
   xsol=[] + [x,]
   ts = []
  
   # SINCE THIS FUNCTION USES ADAPTIVE TIMESTEPPING
   # THE TIME VECTOR IS NOT DEFINED INITIALLY
   t = tstart
   while t < tfinal:
      # RESET INITIAL DT
      h = dt
      xs = []
      ## WHILE ERROR IS ABOVE THRESHOLD
      while error > tol:
         # EVALUATE STAGES
         k = []
         # LOOPS THROUGH ALL K VALUES DEPENDING ON ORDER OF FUNCTION
         for klev in range(0,order):
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


################################################
######## Define inputs for the pred-prey model

import numpy as np
# For pi
import math
# for plotting
import matplotlib.pyplot as plt

a = 1
b = 0.5
c = 0.05
d = 0.02

feval = [lambda x,t: a*x[0,0]-c*x[0,0]*x[1,0], 
      lambda x,t: d*x[0,0]*x[1,0]-b*x[1,0]]
x0 =  [0,1]
tstart = 0
tfinal = 20
dt = (tfinal-tstart)/1000.
order = 4

#tvec,xsol,ts = pagerkck4(feval,x0,tstart,tfinal,dt,ord = 4,tol=10**-6)

################################################
######## test solver without function here
   # INITIALIZE ARRAY OF CASH-KARP FACTORS
ck_a = np.array([[0],[1/5.],[3/10.],[3/5.],[1.],[7/8.]])
ck_b = np.array([[0,0,0,0,0],
            [1/5.,0,0,0,0],
            [3/40.,9/40.,0,0,0],
            [3/10.,-9/10.,6/5.,0,0],
            [-11/54.,5/2.,-70/27.,35/27.,0],
            [1631/55296.,175/512.,575/13824.,44275/110592.,253/4096.]])
ck_c =  np.array([[37/378.],[0.],[250/621.],[125/594.],[0.],[512/1771.]])
ck_cs = np.array([[2825/27648.],[0.],[18575/48384.],[13525/55296.],[277/14336.],[1/4.]])

ck_a = ck_a[range(0,order),:]
ck_b = ck_b[range(0,order),][:,range(0,order)]
ck_c = ck_c[range(0,order),:]
ck_cs = ck_cs[range(0,order),:]

#initialize time vector
tvec=[tstart]
#initialize x, initialize solution matrix xsol
x0 = np.array([x0]).T
x=x0
xsol = np.empty((2,0),float)
xsol = np.append(xsol,x,axis=1)
ts = [0]

# SINCE THIS FUNCTION USES ADAPTIVE TIMESTEPPING
# THE TIME VECTOR IS NOT DEFINED INITIALLY
ord = 4
tol=10**-6
t = tstart

it1 = 0
while t < tfinal and it1 < 100:
   it1 = it1 + 1	
   #print(1)
   # RESET INITIAL DT
   h = dt
   xs = []
   ## WHILE ERROR IS ABOVE THRESHOLD
   error = tol + 1
   it2 = 0
   while error > tol and it2 < 1000:
      it2 = it2 + 1
      #print(2)
      ##IOND
      # EVALUATE STAGES
      f = np.empty((0,1),float)
      for fn in feval:
         #print(4)
         f = np.append(f,np.array([[fn(x,t)]]),axis=0)
      
      
      k = h*f
      #print(str(k))
      # LOOPS THROUGH ALL K VALUES DEPENDING ON ORDER OF FUNCTION
      for klev in range(1,order):
         #print(3)
         #print(str(klev))
         # CALCULATES F FOR ALL FUNCTIONS IN FEVAL
         f = np.empty((0,1),float)
         # EVALUATES F DEPENDING ON LEVEL
         for fn in feval:
            #print(5)
            cb = np.array([ck_b[klev,range(0,(klev))]]).T
            ca = ck_a[klev]
            ###
            ### THE BROKEN THING IS HERE
            ###
            #print(str(cb))
            #print(str(ca))
            #print(k)
            f = np.append(f,np.array([[fn(x+cb.T*k,t+ca)]]),axis = 0)
      
         k = np.append(k,h*f,axis = 1)
         
         ##IOND
      # CALCULATE x values
      xnew = np.sum(x + ck_c.T*k,axis = 1)
      xs = np.sum(x + ck_cs.T*k,axis=1)
      # ESTIMATE ERROR
      error = np.amax(np.absolute(xnew-xs))
      hnew = h*(tol/error)**(0.2)
      # ADJUST DT
      h = hnew
      
   
   x = np.reshape(xnew,(-1,1))
   # STORE SOLUTION (we assume XS is better)
   xsol = np.append(xsol,np.reshape(xs,(-1,1)),axis = 1)
   
   # ADJUST TIME BY H
   t = t + h 
   tvec = tvec + [t]
   ts = ts + [h]
######## Solver ends here
################################################

plt.plot(tvec,xsol[0,:])
plt.plot(tvec,xsol[1,:])
plt.show()



