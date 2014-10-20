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
   while t < tfinal:
      
      ## RESET VALUES
      error = tol*1.01
      errloop = 0
      while errloop == 0:
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
      
   return tvec,xsol,ts

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

umax = 10 # mm
lmax = 20 # mm
ku = 0.3  # day^-1 
kl = 0.05 # day^-1 
kp = 0.1  # day^-1 
tc = 0.5  # day

# du/dl
# dl/dl

# Setup as:
#[0]   u
#[1]   du/dt
#[2]   l
#[3]   du/dt
#[4]   r

x0 =  [1,
       0,
       12,
       0,
       0] 

feval = [lambda x,t: x[0,0],
         lambda x,t: -ku*x[0,0] - kp*(x[2,0]-(x[2,0]/lmax)**3)*x[0,0]-x[4,0]-p,
         lambda x,t: x[2,0],
         lambda x,t: kp*(x[2,0] - (x[2,0]/lmax)**3)*x[0,0]-kl*x[2,0]
         lambda x,t: 0 if x[0,0] < umax else (x[0,0]-umax)/tc
         ]

##### FROM HERE
### check if statements
### Setup p as a global variable
### Add routine to calculate precip at each timestep
### plots



# Model params
tstart = 0
tfinal = 100
dt = (tfinal-tstart)/1000.
order = 3
tol=10**-3
errtol = 10**-20

######## End function input here
################################################


################################################
######## test solver with function here:


tvecf,xsolf,tsf = pagerkck4(feval,x0,tstart,tfinal,dt,tol = tol,order = 3)
plt.plot(tvecf,xsolf[0,:])
plt.plot(tvecf,xsolf[2,:])
plt.show() 

######## End test solver with function
################################################
