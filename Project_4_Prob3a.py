### Page Weil
### CVEN 5537
### 11/14/14
### Project 4, Problem 3 - Sheet Pile - GS/SOR

####### FIX HOW SUBSETS ARE SELECTED AND BOUNDARIES CALCULATED
####### SOMETHING IS BROKEN IN ZONE 5


import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

###################################
###### Begin define problem variables

NX = 100
NY = 101
# TOLERANCE FOR OLD DELTA ERROR METHOD
# tol = 10**-6
# TOLERANCE FOR NEW MASS BALANCE ERROR METHOD
tol = 10**-3.
maxiter = 50000
u0 = 0.
u1 = 1.
jbot = 91 #11, 31, 51, 71, 91
plen = NY-jbot
sploc = int(NX*2/4)-1

# Grid size

xmax = 1.
ymax = 1.
x = np.arange(0,xmax+xmax/(NX-1),xmax/(NX-1))
y = np.arange(0,ymax+ymax/(NY-1),ymax/(NY-1))

k=np.pi
w = 0.5*(1-k/NX)

# initialize matrices
u = np.zeros((NY,NX))
unew = np.zeros((NY,NX))
delta = np.zeros((NY,NX))

# generate set of coordinates for type 1 BCs
# left set
t1z1 = np.hstack((np.zeros((sploc+1,1)),np.zeros((sploc+1,1))))
t1z1[:,1] = np.arange(0,sploc+1)
t1z1[:,0] = NY-1
# right set
t1z2 = np.hstack((np.zeros((NX-sploc-1,1)),np.zeros((NX-sploc-1,1))))
t1z2[:,1] = np.arange(sploc+1,NX)
t1z2[:,0] = NY-1

# generate set of coordinates for type 2 BCs
# zone 1 (left wall)
t2z1 = np.hstack((np.zeros((NY-2,1)),np.zeros((NY-2,1))))
t2z1[:,0] = np.arange(1,NY-1)
t2z1[:,1]= 0

# zone 3 (bottom wall)
t2z3 = np.hstack((np.zeros((NX-2,1)),np.zeros((NX-2,1))))
t2z3[:,1] = np.arange(1,NX-1)
t2z3[:,0]= 0

# zone 5 (right wall)
t2z5 = np.hstack((np.zeros((NY-1,1)),np.zeros((NY-1,1))))
t2z5[:,0] = np.arange(1,NY)
t2z5[:,1]= NX-1

# zone 6 (left side of sheet pile)
t2z6 = np.hstack((np.zeros((plen-1,1)),np.zeros((plen-1,1))))
t2z6[:,0] = np.arange(jbot,NY-1)
t2z6[:,1] = sploc

# zone 7 (right side of sheet pile)
t2z7 = np.hstack((np.zeros((plen-1,1)),np.zeros((plen-1,1))))
t2z7[:,0] = np.arange(jbot,NY-1)
t2z7[:,1] = sploc+1

# homogenous corner point 2
t2z2 = np.array([[0,0]])

# homogenous corner point 4
t2z4 = np.array([[0,NX-1]])

# Initial conditions for type 1 BCs
for a in t1z1:
   u[a[0],a[1]] = u0

for a in t1z2:
   u[a[0],a[1]] = u1

# type 1 corners
u[0,NX-1]=u0
u[NY-1,NX-1]=u1
###### End define problem variables
###################################

###################################
###### Begin GS/SOR solver

# while error > tol
err1 = 10.
it=0
while err1 > tol:
  it = it + 1
  # Loop through interior nodes
  for i in range(0,NY-1):
    for j in range(0,NX):
       if [i,j] in t2z1.tolist():
       # type 2, zone 1          
          delta[i,j] = 4*w/3*(u[i+1,j]+u[i-1,j]+u[i,j+1]-3.*u[i,j])
       elif [i,j] in t2z2.tolist(): 
       # type 2, zone 2          
          delta[i,j] = 4*w/2*(u[i,j+1]+u[i+1,j]-2.*u[i,j])
       elif [i,j] in t2z3.tolist(): 
       # type 2, zone 3          
          delta[i,j] = 4*w/3*(u[i+1,j]+u[i,j+1]+u[i,j-1]-3.*u[i,j])
       elif [i,j] in t2z4.tolist(): 
       # type 2, zone 4          
          delta[i,j] = 4*w/2*(u[i+1,j]+u[i,j-1]-2.*u[i,j])
       elif [i,j] in t2z5.tolist():
          # type 2, zone 5
          delta[i,j] = 4*w/3*(u[i+1,j]+u[i-1,j]+u[i,j-1]-3.*u[i,j])
       elif [i,j] in t2z6.tolist():
          # type 2, zone 6
          delta[i,j] = 4*w/3*(u[i+1,j]+u[i-1,j]+u[i,j-1]-3.*u[i,j])
       elif [i,j] in t2z7.tolist():
          # type 2, zone 7
          delta[i,j] = 4*w/3*(u[i+1,j]+u[i-1,j]+u[i,j+1]-3.*u[i,j])
       else:
          # interior
          delta[i,j] = w*(u[i-1,j]+u[i,j-1]+u[i,j+1]+u[i+1,j]-4.*u[i,j])
       u[i,j] = u[i,j]+delta[i,j]
  # calculate error
  # OLD DELTA ERROR
  # err = np.max(delta)
  # NEW ERROR BASED ON INFLOW-OUTFLOW CONVERGENCE
  outflow = np.sum(u[NY-2,range(0,sploc+1)]-u[NY-1,range(0,sploc+1)])
  inflow = np.sum(u[NY-1,range(sploc+1,NX)]-u[NY-2,range(sploc+1,NX)])
  err1 = abs(outflow-inflow)
  print("outflow = "+str(outflow)+", inflow = "+str(inflow)+", err1= "+str(err1)+", logerr1= "+str(np.log(err1)))
  if round(it,-2) == it:
     print("it: ",str(it),", log(err): ",str(np.log(err1)))
  if it > maxiter: 
    err1 = tol
    print("Hit ",str(maxiter)," iter")


# REMOVE THIS TEMPORARILY

if 0 == 1:
   # Set values for exterior nodes
   # type 2 zone 1
   for a in t2z1:
      u[a[0],a[1]-1] = u[a[0],a[1]]
   
   # type 2 zone 2
   for a in t2z2:
      u[a[0],a[1]-1] = u[a[0],a[1]]
      u[a[0]-1,a[1]] = u[a[0],a[1]]
      u[a[0]-1,a[1]-1] = u[a[0],a[1]]
   
   # type 2 zone 3
   for a in t2z3:
      u[a[0]-1,a[1]] = u[a[0],a[1]]
   
   # type 2 zone 4
   for a in t2z4:
      u[a[0],a[1]+1] = u[a[0],a[1]]
      u[a[0]-1,a[1]] = u[a[0],a[1]]
      u[a[0]-1,a[1]+1] = u[a[0],a[1]]
   
   # type 2 zone 5
   for a in t2z5:
      u[a[0],a[1]+1] = u[a[0],a[1]]
   
   # type 2 zone 6
   for a in t2z6:
      u[a[0],a[1]+1] = u[a[0],a[1]]
   
   # type 2 zone 7
   for a in t2z7:
      u[a[0],a[1]-1] = u[a[0],a[1]]
   
   print("GS took ",str(it)," iterations")
# END REMOVE THIS TEMPORARILY

###### End GS/SOR solver
###################################

###################################
###### Begin plotting code here

# Make plottable arrays (X,Y)
# Coordinate grids 
X, Y = np.meshgrid(x, y)

CS = plt.contourf(X, Y, u, 25, cmap=plt.cm.bone)

plt.title("Problem 4, jbot= "+str(jbot))
plt.xlabel('x')
plt.ylabel('y')
plt.show()

#outflow - inflow
print("jbot= "+str(jbot)+", outflow = "+str(outflow)+", inflow = "+str(inflow)+", net out= "+str(outflow - inflow))
###### End plotting code
###################################

# Plot net mass balance vs sheet pile location
if 1==0:
   netout = [0.00698,0.01,0.0138,0.021,0.642]
   #outflow = [0.184,0.259,0.343,0.481,0.818]
   jbot = (NX-np.array([11,31,51,71,91.]))/NX
   
   # Plot outflow vs sheet pile location (how does sheetpile impact total outflow)
   plt.plot(jbot,outflow)
   plt.xlabel('Sheetpile height as pct of total')
   plt.ylabel('Outflow')
   plt.show()
   
   # Plot net mass balance vs sheet pile location
   plt.plot(jbot,netout)
   plt.xlabel('Sheetpile height as pct of total')
   plt.ylabel('Net outflow')
   plt.show()
   


