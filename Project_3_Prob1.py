### Page Weil
### CVEN 5537
### 11/3/14
### Project 3, Problem 1

## EF Matlab code from class ported to Python:
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

start1 = datetime.now()

# Euler forward Matlab Code
#set parameters b, deltax and J = number of nodes 

# setting b=1 "dummy value" just for illustration
b=1 
deltax=0.01 

# L=1, hence 1/deltax
J=1/deltax+1 

# mu is set here......
mu=0.4 

#finding deltat based on mu
deltat=mu*deltax**2/b

#the total time I want is 0.5 units, so number of timesteps nt = 0.5/deltat
nt=round(0.5/deltat,0)

#define initial condition - example 1
u=np.ones((J,1))

#define x vector for plotting and plot initial condition
x=np.arange(0,1+deltax,deltax)
#plt.plot(x,u)
#plt.show()

unext = u

# hold on%holds the plot frame so we can add other plots (see below)
#time loop starts here
for it in np.arange(0,nt)+1:
    #in EF, we first do interior nodes and then the B.C.
    #EF for interior nodes
    for j in np.arange(1,J-1): # this is actually J-1
        unext[j,0]=u[j,0]+mu*(u[j-1,0]-2*u[j,0]+u[j+1,0])
    
    #B.C.
    #boundary conditions - homogeneous first type at both ends for this
    #example, so quite easy
    unext[0,0]=0
    unext[J-1,0]=0
    #update u for next step as unext from this step
    u = unext
    #plot at selected time-steps

end1 = datetime.now()

###################################
###### Euler forward Matlab Code (With matrix operations instead of loops)
start2 = datetime.now()
#set parameters b, deltax and J = number of nodes 

# setting b=1 "dummy value" just for illustration
b=1 
deltax=0.01 

# L=1, hence 1/deltax
J=1/deltax+1 

# mu is set here......
mu=0.4 

#finding deltat based on mu
deltat=mu*deltax**2/b

#the total time I want is 0.5 units, so number of timesteps nt = 0.5/deltat
nt=round(0.5/deltat,0)

#define initial condition - example 1
u=np.ones((J,1))

#define x vector for plotting and plot initial condition
x=np.arange(0,1+deltax,deltax)
#plt.plot(x,u)
#plt.show()

unext = u

# hold on%holds the plot frame so we can add other plots (see below)
#time loop starts here
for it in np.arange(0,nt)+1:
    #in EF, we first do interior nodes and then the B.C.
    #EF for interior nodes
    unext[1:(J-1),:] = np.array([u[1:(J-1),0]+mu*(u[0:(J-2),0]-2*u[1:(J-1),0]+u[2:(J),0])]).T
    
    #B.C.
    #boundary conditions - homogeneous first type at both ends for this
    unext[0,0]=0
    unext[J-1,0]=0
    #update u for next step as unext from this step
    u = unext
    #plot at selected time-steps

end2 = datetime.now()

print("No matrix ops: "+str(end1-start1)+", with matrix ops: "+str(end2-start2))

# Results: No matrix ops: 0:00:48.520919, with matrix ops: 0:00:02.227821
