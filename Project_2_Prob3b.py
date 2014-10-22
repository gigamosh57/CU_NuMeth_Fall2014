#### Page Weil
#### 10/21/14
#### CVEN 5537
#### PROJECT 2 , PROBLEM 3B
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
      xnew = x
      xsol = np.empty((len(x0),0),float)
      xsol = np.append(xsol,x0,axis=1)
      ts = [0]
      
      # SINCE THIS FUNCTION USES ADAPTIVE TIMESTEPPING
      # THE TIME VECTOR IS NOT DEFINED INITIALLY
      t = tstart
      h = dt
      errvec = []
      global r
      p = 0
      rvec = [0]
      prec = [p] 
      while t < tfinal:
            if round(t,1) == round(t,2): print( "t: " + str(round(t,1)) + ", h: " + str(round(h,2)))
            
            ## RESET VALUES
            error = tol*1.01
            errloop = 0
            while errloop == 0:
               p = precip(t)
               if xnew[0] <= umax:
                  r = 0
               else:
                  r = (xnew[0]-umax)/tc
               
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
               
               error = np.amin([np.amax(abserror),np.amax(relerror)])
               
               # ADJUST H BASED ON ERROR
               hnew = h*(tol/error)**(0.2)
               if hnew > h*1.2: hnew=h*1.2
               h = hnew
               if hnew >= h: errloop = 1
               
               #print 'error: ' + str(error)
               ### End while errloop = 0
            prec = prec + [p]
            rvec = rvec + [r]
            # UPDATE ALL VALUES FOR NEXT TIMESTEP
            h = hnew
            x = np.reshape(xnew,(-1,1))
            xsol = np.append(xsol,x,axis = 1)
            t = t + h 
            tvec = tvec + [t]
            ts = ts + [h]
            ### End while time < tfinal
      
      return tvec,xsol,ts,prec,rvec

######## Solver ends here
################################################

####################################################################
######## Precip function
def precip(t):
      ## CALCULATE PRECIP
   global p
   if t >= 0 and t <= 1:
      p = 1
   elif t > 1 and t <= 2:
      p = 5
   elif t > 2 and t <= 3:
      p = 5
   elif t > 8 and t <= 9:
      p = 3
   elif t > 9 and t <= 10:
      p = 10
   elif t > 10 and t <= 11:
      p = 12
   elif t > 11 and t <= 12:
      p = 6
   elif t > 12 and t <= 13:
      p = 1
   else:
      p = 0
   return p
######## end precip   
####################################################################

################################################
######## Define inputs for the hydrologic model problem

# Problem params
global umax
umax = 10 # mm
lmax = 20 # mm
ku = 0.3  # day^-1 
kl = 0.05 # day^-1 
kp = 0.1 # day^-1 
global tc
tc = 0.5  # a

# Setup as:
#[0]   du/dt
#[1]   dl/dt

x0 =  [1,
       12]

feval = [lambda x,t: -ku*x[0,0] - kp*(1-(x[1,0]/lmax)**3)*x[0,0] - r + p,
         lambda x,t: kp*(1 - (x[1,0]/lmax)**3)*x[0,0]-kl*x[1,0]]

# Model params
tstart = 0
tfinal = 20
dt = (tfinal-tstart)/1000.
order = 4
tol=10**-
errtol = 10**-20

######## End function input here
################################################

################################################
######## test solver with function here:

tvecf,xsolf,tsf,prec,rplot = pagerkck4(feval,x0,tstart,tfinal,dt,tol = tol,order = order)
fig = plt.figure()
axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])
axes.set_xlabel('time')
axes.set_ylabel('mm')

uplot = xsolf[0,:]
lplot = xsolf[1,:]
splot = uplot + lplot
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

################################################
######## Discussion
# The equation is not set up to limit the maximum value of u to less than umax
# As soon as u passes umax, the runoff kicks in and u begins to drop.
# This hydrograph shows an immediate response to precipitation by the upper zone
# u (blue line is u volume) and a very muted response by the lower zone (green line is l volume).
# 
# both upper and lower zones fill quickly and release water slowly.
# Runoff occurs when the upper zone is full, even if the lower zone is not.  This is because
# the lower zone cannot transmit water through interflow or percolation quickly enough
# to outpace the precip
