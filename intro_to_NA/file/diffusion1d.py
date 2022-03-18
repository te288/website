## Import Modules
import numpy as np
from scipy.stats import hmean
import matplotlib.pyplot as plt

def PlotPressure(x, P, t, L):
  # Function to Plot and Save Graph
  # fig = plt.figure()
  plt.plot(x, P_new, label='t={0:05.2f}'.format(t)) 
  plt.xlabel('x[m]')
  plt.ylabel('Pressure [Pa]')
  plt.xlim(0, L)
  plt.ylim(-1,1)
  plt.grid()
  plt.title('Pressure Diffusion 1D@{0:05.2f}[s]'.format(t))
  return

## Input Parmeters
L = 2*np.pi # Length of Reservoir[m]
N = 100 # Number of Control Volume[-]
k = 0.2*np.ones(N) # Permiability[m^2]
phi = np.ones(N) # Porosity[m^2]
c = np.ones(N) # Compressibility[Pa^-1]
mu= 1 # Viscosity of Fluid[Pa^-1]
dx = L/N # Size of Control Volume[m]
x = np.zeros(N) # x coordinate[m]
x[0] = dx/2
for i in range(1, N):
  x[i] = x[i-1] + dx

## Parameters for Simulation
tmax = 20 # Time to stop simlation [s]
dt   = 0.005   # dt [s]
nout = 500  # output result ever nout step[s]
# Variable to decide Boundary condition
# 1 -> Neumann, 0 -> Dirichlet
B_right = 0 # Boundary at right(x = 0) 
B_left  = 0 # Boundary at Left (x = L)
Pb_right = 0 # Pressure Value at x = 0
Pb_left  = 0 # Pressure Value at x = L

## Initial Conditions
# P_init = np.ones(N)      # Initial Pressure
P_init = np.sin(x)

## Simulation
P_old  = np.copy(P_init) # Pressure at n-th step
P_new  = np.copy(P_init) # Pressure at n+1-th step
t = 0
n = 0

# Plot Initial Condition
fig = plt.figure()
PlotPressure(x, P_new, t, L)
while True:
  for i in range(1, N-1): # for P[1] ~ P[N-2]
    alpha = #-# Write Your Code Here #-#
    A = #-# Write Your Code Here #-#
    B = #-# Write Your Code Here #-#
    C = #-# Write Your Code Here #-#
    P_new[i] = A*P_old[i+1] + B*P_old[i] + C*P_old[i-1]
  
  # P[0]
  if B_left == 1: # Neumann Condition
    #-# Write Your Code Here #-#

  else: # Dirichlet Condition 
    #-# Write Your Code Here #-#

  # P[N-1]
  if B_right  == 1: # Neumann Condition
    #-# Write Your Code Here #-#
    
  else: # Dirichlet Condition 
    #-# Write Your Code Here #-#
    
  # Update Values, time step and Add plot
  P_old = np.copy(P_new)
  t = t + dt
  n = n + 1
  if t >= tmax:
    break
  if n%nout == 0:
    print('{0:04d}th Time step {1:05.2f}'.format(n, t))
    PlotPressure(x, P_new, t, L)

    
PlotPressure(x, P_new, t, L)
plt.show()

