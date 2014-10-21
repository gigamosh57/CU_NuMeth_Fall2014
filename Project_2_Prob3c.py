####################################################################
######## INITIALIZE PYTHON
import numpy as np
# For pi
import math
# for plotting
import matplotlib.pyplot as plt
######## 
####################################################################

####################################################################
######## Bisection solver
def proot(f,bounds,n = 1000,tol = 0.01): 
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
  return xmid, f(xmid) # Location and value
######## 
####################################################################


####################################################################
######## Runge-Kutta ODE solver with Cash-Karp 4th-5th order params

def pagerkck4(feval,x0,tstart,tfinal,dt,order = 4,tol=10**-2,errtol = 10**-20):
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
   xsol = np.empty((len(x0),0),float)
   xsol = np.append(xsol,x0,axis=1)
   ts = [0]
   
   # SINCE THIS FUNCTION USES ADAPTIVE TIMESTEPPING
   # THE TIME VECTOR IS NOT DEFINED INITIALLY
   t = tstart
   h = dt
   errvec = []
   qvec = []
   Q = 0
   p = 1
   while t < tfinal:
      # Print progress bar
      if round(t,1) == round(t,2): print( "t: " + str(t))
      ## RESET VALUES
      error = tol*1.01
      errloop = 0
      
      while errloop == 0:
         
         # Initial Guess for V
         v = np.sqrt(2*grav*x[0])
         # Calculate Q and f
         global currH = x[0]
         fv = lambda x: np.sqrt(2*grav*currH/(1+x*L/D))
         fric,v = proot(fv,v)
         
         
         # EVALUATE STAGES
         f = np.empty((0,1),float)
         for fn in feval:
            f = np.append(f,fn(x,t))
         
         k = np.array([h*f]).T
         
         # LOOPS THROUGH ALL K VALUES DEPENDING ON ORDER OF FUNCTION
         for klev in range(1,order):
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
         
         # CALCULATE x values
         xnew = np.sum(x + np.array([np.sum(ck_c.T*k,axis = 1)]).T,axis = 1)
         xs = np.sum(x + np.array([np.sum(ck_cs.T*k,axis=1)]).T , axis = 1)
            
         # ESTIMATE ERROR
         abserror = np.absolute(xnew-xs)
         relerror = np.zeros((len(x0),1))
         for i in range(0,len(xnew)): 
            if xnew[i] > errtol: relerror[i] = np.absolute((xnew[i] - xs[i])/xnew[i])
         
         error = np.amax([np.amax(abserror)])#,np.amax(relerror)])
         
         # ADJUST H BASED ON ERROR
         hnew = h*(tol/error)**(0.2)
         if hnew > h*1.2: hnew=h*1.2
         h = hnew
         if error < tol: errloop = 1
         
         #print 'error: ' + str(error)
         ### End while errloop = 0
      
      # UPDATE ALL VALUES FOR NEXT TIMESTEP
      h = hnew
      x = np.reshape(xnew,(-1,1))
      xsol = np.append(xsol,x,axis = 1)
      t = t + h 
      tvec = tvec + [t]
      ts = ts + [h]
      ### End while time < tfinal
      
   return tvec,xsol,ts,prec

######## Solver ends here
################################################



####################################################################
######## INITIALIZE PYTHON
import numpy as np
# For pi
import math
# for plotting
import matplotlib.pyplot as plt
######## 
####################################################################


################################################
######## Define inputs for the hydrologic model problem

# Problem params

global Atank = 0.3  # m
global D =     0.01 # m
global nu = 10**-6  # m2/s
global L =     2    # m
global grav =  9.81 # m/s2

# Setup as:
#[0]   H
#[1]   dH/dt

x0 =  [3,
       0] 

feval = [lambda x,t:  x[0,0],
         lambda x,t: -Q/Atank]

##### FROM HERE
### check if statements
### Setup p as a global variable
### Add routine to calculate precip at each timestep
### plots



# Model params
tstart = 0
tfinal = 1800
dt = (tfinal-tstart)/1000.
order = 3
tol=10**-2
errtol = 10**-20

######## End function input here
################################################


################################################
######## test solver with function here:


tvecf,xsolf,tsf,prec = pagerkck4(feval,x0,tstart,tfinal,dt,tol = tol,order = 3)
fig = plt.figure()
axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])
axes.set_xlabel('T')
axes.set_ylabel('mm')

uplot = xsolf[0,:]
lplot = xsolf[2,:]
rplot = xsolf[4,:]
qplot = ku*uplot + kl*lplot + rplot
axes.plot(tvecf,uplot,label = "u")
axes.plot(tvecf,lplot,label = "l")
axes.plot(tvecf,rplot,label = "r")
axes.plot(tvecf,qplot,label = "q")
axes.plot(tvecf,prec,label = "precip")

axes.legend()
plt.show() 

######## End test solver with function
################################################
