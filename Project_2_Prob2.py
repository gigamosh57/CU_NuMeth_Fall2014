####################################################################
 ######## Runge-Kutta ODE solver with Cash-Karp 4th-5th order params
 
#def pagerkck4(feval,x0,tstart,tfinal,dt,ord = 4,tol=10**-6):

#### Solver code goes here

#   return tvec,xsol,ts
######## Solver ends here
################################################


################################################
######## Define inputs for the pred-prey model

import numpy as np
# For pi
import math
# for plotting
import matplotlib.pyplot as plt

a = 1    # 1    # 5
b = 0.5  # 0.5  # 1
c = 0.05 # 0.05 # 1
d = 0.02 # 0.02 # 1

feval = [lambda x,t: a*x[0,0]-c*x[0,0]*x[1,0], 
      lambda x,t: d*x[0,0]*x[1,0]-b*x[1,0]]
x0 =  [100,10]
tstart = 0
tfinal = 100
dt = (tfinal-tstart)/10000.
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
xsol = np.append(xsol,x0,axis=1)
ts = [0]

# SINCE THIS FUNCTION USES ADAPTIVE TIMESTEPPING
# THE TIME VECTOR IS NOT DEFINED INITIALLY
tol=10**-3
errtol = 10**-20
t = tstart
h = dt

it1 = 0
errvec = []
while t < tfinal:# and it1 < 20:
#   if round(t,2) == round(t,1): print('t: '+str(t))
   it1 = it1 + 1	
   #print(1)
   
   ## RESET VALUES
   error = tol*1.01
   it2 = 0
   xs = []
   xn = []
   hvec = []
   #h = dt
   
   xnvec = []
   xsvec = []
   errloop = 0
   #while error > tol and it2 < 5:
   while errloop == 0:# and it2 < 10:
      it2 = it2 + 1
      # EVALUATE STAGES
      #id 6      
      f = np.empty((0,1),float)
      for fn in feval:
         #print(4)
         f = np.append(f,fn(x,t))
      
      k = np.array([h*f]).T
      #print(str(k))
      
      # LOOPS THROUGH ALL K VALUES DEPENDING ON ORDER OF FUNCTION
      for klev in range(1,order):
         #indent 3
         # CALCULATES F FOR ALL FUNCTIONS IN FEVAL
         f = np.empty((0,1),float)
         # EVALUATES F DEPENDING ON LEVEL OF K
         for fn in feval:
            cb = np.array([ck_b[klev,range(0,klev)]])
            ca = ck_a[klev]
            kev = k[:,range(0,klev)]
            f = np.append(f,fn(x+np.array([np.sum(cb*kev,axis=1)]).T,t+ca*h))
         
         f = np.array([f]).T
         k = np.append(k,h*f,axis = 1)
         # End for klev
         #indent 3
      # id 6
      
      # CALCULATE x values
      xnew = np.sum(x + np.array([np.sum(ck_c.T*k,axis = 1)]).T,axis = 1)
      xs = np.sum(x + np.array([np.sum(ck_cs.T*k,axis=1)]).T , axis = 1)
      
      xnvec = np.append(xnvec,xnew)
      xsvec = np.append(xsvec,xs)
      
      # ESTIMATE ERROR
      abserror = np.absolute(xnew-xs)
      relerror = np.zeros((2,1))
      for i in range(0,len(xnew)): 
         if xnew[i] > errtol: relerror[i] = np.absolute((xnew[i] - xs[i])/xnew[i])
      
      #error = np.amin([np.absolute(np.amax(tol/abserror)),np.absolute(np.amax(tol/relerror))])
      error = np.amax(np.amax(abserror),np.amax(relerror))
      #error  = np.amax(np.absolute(xnew-xs))
      
      # ADJUST H BASED ON ERROR
      hnew = h*(tol/error)**(0.2)
      hvec = hvec + [hnew]
      if hnew > h*1.2: hnew=h*1.2
      h = hnew
      if error < tol: errloop = 1
      
      #print('error: ' + str(error))
      
      
      #print('klev: ' + str(klev))
      #print('k: ' + str(k))
      #error = tol - 0.1
      ### End while error > tol
   
   h = hnew
   #print("iter err: " + str(it2))
   #errvec = errvec + [error]      
   #print('h: ' + str(h))
   #print('error: ' + str(error))
   #print('f: ' + str(f))
   #print('xnew: ' + str(xnew))
   #print('xs: ' + str(xs))
   
   # STORE VALUES FROM ADAPTIVE TIMESTEPPING
   x = np.reshape(xnew,(-1,1))
   
   # STORE SOLUTION (we assume XNEW is better)
   xsol = np.append(xsol,x,axis = 1)
   #print('x = '+str(x))
   #print('t = '+str(round(t,4)))
   # ADJUST TIME BY H
   t = t + h 
   tvec = tvec + [t]
   ts = ts + [h]
   
   # Set starting h for next timestep (max increase)
   #h = 1.2*h
   
   ### End while time < tfinal
   
   
######## Solver ends here
################################################



pt = range(1,len(tvec)-1)
tvec = [tvec[i] for i in pt]
ts = [ts[i] for i in pt]
xsol1 = [xsol[0,i] for i in pt]
xsol2 = [xsol[1,i] for i in pt]
#errvec = [errvec[i] for i in pt]
plt.plot(tvec,xsol1)
plt.plot(tvec,xsol2)
plt.show() 


#plt.plot(tvec,xsol[0,:])
#plt.plot(tvec,xsol[1,:])
#plt.show()


