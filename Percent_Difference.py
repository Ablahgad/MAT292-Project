import scipy.io as sio
import numpy as np
import math
from Finite_difference import vx, vy, X, Y

'''
Run with:
    domain_x = 30
    domain_y = 10
    divisions_x = 25
    divisions_y = 21
'''

# Loading experimental data
experimental = sio.loadmat("Experimental_Data/velocity_cruise.mat")

# Experimental vx --> u
# Experimental vy --> v

u = experimental["u"]
v = experimental["v"]


# Our solver blows up too fast to consider any time too far into the future so just comparing one time step towards the start, as the fish just starts to move
v = v[1]
v = v.T
u = u[1]
u = u.T

# Scale experimental data to match dimensions of our calculated data
u = u.reshape(100//2, 2, 84//2, 2).mean(axis=(1, 3))
v = v.reshape(100//2, 2, 84//2, 2).mean(axis=(1, 3))
u = u*1000
v = v*1000

p_diffx = (vx-u)/u
p_diffy = (vy-v)/v
p_diff = np.sqrt(np.square(p_diffx) + np.square(p_diffx**2))*100

np.savetxt('Percent_Difference/Percent_difference_data.txt', p_diff, delimiter='   ')
np.savetxt('Percent_Difference/X.txt', X, delimiter='   ')
np.savetxt('Percent_Difference/Y.txt', Y, delimiter='   ')