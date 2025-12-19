'''
This is the model code specified for the finite difference numerical analysis model
'''


#import libraries
from scipy.integrate import quad
import numpy as np
import math

'''
SETTINGS
'''
duration = 0.002 # How long the CFD runs until, in seconds

dt = 0.001 # Time step for movement of fish tail in seconds

v_i = 0 # Initial speed of the flow in the swimming tube

# real domain dimensions (mm)
domain_x = 30 # full x dimension of tube is 150mm but full scale is too large to see the fish properly
domain_y = 10 # full y dimension of tube is 45.7
# Note that the actual domain will be from [-domain_x, domain_x] and the range will be [-domain_y, domain_y]

divisions_x = 15
divisions_y = 15
 # Keep in mind the number of divisions in the each direction will actually be double these values due to mirroring across the x and y axes


'''
CONSTANTS
'''
rho = 0.001 #density in kg/mm^3
kv = 1.004 # kinematic viscosity in mm/s

# nu is time penalty in Brinkman Penalization, determines how fast the velocity changes to the body velocity, the larger it is the more fluid it lets into the fish, smaller it is the faster the velocity change (recommended between 0.1 and 1 -- it is set lower here because we wanted to emphacise the impact of the tail on the fluid)
nu = 0.005


# Creates list of all x values and y values from the domain and number of divisions

d_x = (domain_x/(divisions_x - 1))
d_y = (domain_y/(divisions_y - 1))
x = [d_x*i for i in range(divisions_x)]
x_alt = [i*-1 for i in x]
x.reverse()
for i in x_alt:
    x.append(i)
x.reverse()
# x contains all of the x coordinates of the nodes

y = [d_y*i for i in range(divisions_y)]
y_alt = [i*-1 for i in y]
y.reverse()
for i in y_alt:
    y.append(i)
y.reverse()
# y contains all of the y coordinates of the nodes

'''
INITIALIZING ARRAYS
'''

X, Y = np.meshgrid(x, y, indexing='ij') #array that holds x and y postion values in mm for each node

# Speed of the water at each node
vx = np.full((divisions_x*2, divisions_y*2), v_i)

vy = np.zeros((divisions_x*2, divisions_y*2))

# Speed of the body of the fish
vbodyy = np.zeros((divisions_x*2, divisions_y*2))

vbodyx = np.zeros((divisions_x*2, divisions_y*2))

# Pressure
P = np.zeros((divisions_x*2, divisions_y*2))

# Equal to one inside the fish and zero outside of the fish
chi = np.zeros((divisions_x*2, divisions_y*2))

# Calculated new velocities before pressure correction
v_starx = np.zeros((divisions_x*2, divisions_y*2))

v_stary = np.zeros((divisions_x*2, divisions_y*2))

# Forces on the fluid at each node
Fx = np.zeros((divisions_x*2, divisions_y*2))

Fy = np.zeros((divisions_x*2, divisions_y*2))

# Constant to assist with pressure correction
b = np.zeros((divisions_x*2, divisions_y*2))

'''
Initializing parameters for tail equations
'''

fish_length = 34
head = -1*fish_length/2
tail = fish_length
    # Head is the starting point of the head at the left, head+tail is the end point of the tail
    # Ex. head = -4.5, tail = 9 means that the fish starts at -4.5 and ends at 5

thickness = 3.5
    # Thickness is the width of the thickest point of the fish, ex. thickness = 1

'''
USEFUL FUNCTIONS FOR MODELLING MOTION OF FISH
'''
# Defining the function of x with respect to the parameter s
def f_x(s, head, tail):
    return head+tail*s

# Defining the function of s with respect to x, same head and tail as in f_x
def f_s(x, head, tail):
    return (x-head)/tail

def f_y(s, thickness): # Defines where the thickness of the body of the fish along its length
    return thickness/0.2*(0.2969*math.sqrt(s) - 0.126*s - 0.3516*s**2 + 0.2843*s**3-0.1036*s**4)
    # Based off of NACA thickness equation

def y_tail(x, t): # Defines the y value of each x point along the tail at any given time
    return (-1.2958*x + 0.5105*x**2)*math.sin(2*math.pi*(x/29.766+0.25*t))/10 + thickness/4

def v_tail(x, t): # Defines the speed of any point along the tail at any given time, the derivative of y_tail with respect to time
    return (-1.2958*x + 0.5105*x**2)*2*math.pi*0.25*math.cos(2*math.pi*(x/29.766+0.25*t))/10

def dy_dx_tail(x, t): # Used to calculate the length of the tail using the arc length formula, the derivative of y_tail with respect to x
    return (-1.2958*x + 0.5105*x**2)*2*math.pi/29.766*math.cos(2*math.pi*(x/29.766+0.25*t))/10 + (-1.2958 + 0.5105*x*2)*math.sin(2*math.pi*(x/29.766+0.25*t))/10

def tail_length_helper(x, t): # the function within the integral for the arc length equation -- used to calculate length of tail
    return math.sqrt(1+dy_dx_tail(x, t)**2) #using the v_tail formula here as the derivative of

'''
MAIN LOOP
'''

# list of file names that will be read by MATLAB
L_file_names = [] # List of files with info on nodes outside of the fish
L_fish_file_names = [] # List of files with info on nodes inside of the fish


for t in np.arange(0, duration, dt):

    file_text = ""  # Holds the information for nodes outside of the fish
    fish_text = "" # Holds the information for nodes inside of the fish
    index = 0
    for i in range(1, len(x)-1):
        for j in range(1, len(y)-1):

            '''
            Updating where the fish is
            '''
            # resetting forces and body velocity for each node
            Fx[i, j] = 0
            Fy[i, j] = 0
            vbodyx[i, j] = 0
            vbodyy[i, j] = 0

            s_i = f_s(x[i], head, tail)
            if s_i < 1 and s_i > 0:
                if x[i] > 0:
                    y_i = y_tail(x[i], t)

                    #using arc length formula to find length of tail at x coordinate x[i] starting from the base of the tail (x=0)
                    tail_length = quad(tail_length_helper, 0, x[i], args=(t))[0]

                    if (tail_length<=fish_length/2): # ensures length of tail remains equal to 1/2 fish_length even as it moves
                        y_i_bot =  y_tail(x[i], t)-f_y(s_i, thickness)-thickness/4 # If all conditions for the tail to exist are satisfied the thickness of the tail is subtracted from the y tail position and centered to the body of the fish by subtracting thickness/4 to create the lower bounds of the tail
                    else:
                        y_i_bot = y_i # If tail is too long it will have no thickness at this x value (won't appear or affect results)

                else:
                    y_i = f_y(s_i, thickness)
                    y_i_bot = -1*f_y(s_i, thickness)
                    # If the x coordinate isn't where the tail is meant to be (x<0) then the lower bound is just the upper bound reflected in the x-axis to form the head of the fish


                if y[j]>y_i_bot and y[j]<y_i: #If coordinate is within the body of the fish
                    chi[i, j] = 1


                    if x[i] > 0: # If coordinate is also specfically part of the tail
                        vbodyy[i, j] = v_tail(x[i], t) # Setting the body speed of the tail


                    '''
                    BRINKMAN PENALIZATION
                    '''
                    Fx[i, j] = -1*chi[i, j]/nu*(vx[i, j]-vbodyx[i, j])
                    Fy[i, j] = -1*chi[i, j]/nu*(vy[i, j]-vbodyy[i, j])

                    # Add node to data to be read by MATLAB
                    fish_text += f"{X[i, j]} {Y[i, j]} {chi[i, j]} {vx[i, j]} {vy[i, j]} {vbodyx[i, j]} {vbodyy[i, j]} \n"

                else: # Otherwise the node is not inside the fish so chi is equal to zero
                    chi[i, j] = 0
            else:
                chi[i, j] = 0

            '''
            FINITE DIFFERENCE ESTIMATIONS OF SLOPES
            '''
            dvx_dx = (vx[i+1, j] - vx[i-1, j])/2/d_x
            dvx_dy = (vx[i, j+1] - vx[i, j-1])/2/d_y


            dvy_dx = (vy[i+1, j] - vy[i-1, j])/2/d_x
            dvy_dy = (vy[i, j+1] - vy[i, j-1])/2/d_y


            # Laplacian operator
            D2vx = (vx[i+1, j] -2*vx[i, j] + vx[i-1, j])/d_x**2 + (vx[i, j+1] -2*vx[i, j] + vx[i, j-1])/d_y**2
            D2vy = (vy[i+1, j] -2*vy[i, j] + vy[i-1, j])/d_x**2 + (vy[i, j+1] -2*vy[i, j] + vy[i, j-1])/d_y**2

            dP_dx = (P[i+1, j] - P[i-1, j])/2/d_x
            dP_dy = (P[i, j+1] - P[i, j-1])/2/d_y

            '''
            ACTUAL NAVIER-STOKES
            '''
            dv_dtx = -1*(vx[i, j]*dvx_dx + vy[i, j]*dvx_dy) - 1/rho*dP_dx + kv*(D2vx) + 1/rho*Fx[i, j]
            dv_dty = -1*(vx[i, j]*dvy_dx + vy[i, j]*dvy_dy) - 1/rho*dP_dy + kv*(D2vy) + 1/rho*Fy[i, j]


            v_starx[i, j] = vx[i, j] + dt*dv_dtx
            v_stary[i, j] = vy[i, j] + dt*dv_dty


            '''
            JACOBIAN ITERATION
            '''

            for n in range(1, len(x)-1):
                for m in range(1, len(y)-1):
                    b[n, m] = rho/dt*((v_starx[n+1, m]-v_starx[n-1, m])/2/d_x + (v_stary[n, m+1]-v_stary[n, m-1])/2/d_y)

            for k in range(10):
                for n in range(1, len(x)-1):
                    for m in range(1, len(y)-1):
                        P[n, m] = (P[n+1, m] + P[n-1, m] + P[n, m+1] + P[n, m-1] - b[n, m]*d_x**2)/4
                        vx[0, m] = 26

            '''
            PRESSURE CORRECTED V
            '''
            vx[i, j] = v_starx[i, j] - dt/rho*(P[i+1, j] - P[i-1, j])/2/d_x
            vy[i, j] = v_stary[i, j] - dt/rho*(P[i+1, j] - P[i-1, j])/2/d_y

            if chi[i, j] == 0: # Adding nodes that are outside of the fish to their file_text variable to be read by MATLAB
                file_text += f"{X[i, j]} {Y[i, j]} {chi[i, j]} {vx[i, j]} {vy[i, j]} {vbodyx[i, j]} {vbodyy[i, j]} \n"

    # Saves node information to text files
    with open("dataout/nodes" +str(t)+ ".txt", "w") as file:
        file.write(file_text)
    with open("dataout/fish" +str(t)+ ".txt", "w") as file:
        file.write(fish_text)

    L_file_names.append("dataout/nodes" +str(t)+ ".txt")
    L_fish_file_names.append("dataout/fish" +str(t)+ ".txt")

# Creates variable used by MATLAB to get all the files names
L_files = []
L_files.append(L_file_names)
L_files.append(L_fish_file_names)
