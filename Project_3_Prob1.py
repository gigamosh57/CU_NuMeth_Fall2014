### Page Weil
### CVEN 5537
### 11/3/14
### Project 3, Problem 1

## EF Matlab code from class ported to Python:
import numpy as np
import matplotlib.pyplot as plt

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
plt.plot(x,u)
plt.show()

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

# Turn plot off  
# if (it==fix(nt/10)) || (it==fix(nt/5)) || (it==fix(nt/4)) || (it==fix(nt/2)) || (it == nt)
#     plot(x,u)
# end
# end
# hold off

