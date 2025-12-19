import scipy.io as sio
import numpy as np
from Finite_difference import vx, vy

# Loading experimental data
experimental = sio.loadmat("Experimental_Data/velocity_cruise.mat")

# Experimental vx --> u
# Experimental vy --> v

u = experimental["u"]
v = experimental["v"]