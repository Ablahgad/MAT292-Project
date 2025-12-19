'''
This is the model code specified for the finite difference numerical analysis model
'''


#import libraries
from scipy.integrate import quad
import numpy as np
import math

'''
CONSTANTS
'''
rho = 1 #density
kv = 1 # kinematic viscosity

# nu is time penalty in Brinkman Penalization, determines how fast the velocity changes to the body velocity, the larger it is the more fluid it lets into the fish, smaller it is the faster the velocity change (recommended between 0.1 and 1)
nu = 0.005

#input section - real domain dimensions (mm), number of nodes along x, number of nodes along y
domain_x = 51.6 # full x dimension of tube is 150mm but full scale is too large to see the fish properly
domain_y = 26.7 # full y dimension of tube is 45.7

divisions_x = 14
divisions_y = 7

D_nodes = {} #dictionary that stores all the nodes, keys are str(x_coordinate) + " " + str(y_coordinate)

#create matrix with x and y intersections (nodal locations)
nodes = np.zeros(((divisions_x * divisions_y), 2))

#initialize loop variables

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

vx = np.full((divisions_x*2, divisions_y*2), 26)

vy = np.zeros((divisions_x*2, divisions_y*2))

vbodyy = np.zeros((divisions_x*2, divisions_y*2))

vbodyx = np.zeros((divisions_x*2, divisions_y*2))

P = np.zeros((divisions_x*2, divisions_y*2))

chi = np.zeros((divisions_x*2, divisions_y*2))

v_starx = np.zeros((divisions_x*2, divisions_y*2))

v_stary = np.zeros((divisions_x*2, divisions_y*2))

Fx = np.zeros((divisions_x*2, divisions_y*2))

Fy = np.zeros((divisions_x*2, divisions_y*2))

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


# Defining the function of x with respect to the parameter s
def f_x(s, head, tail):
    return head+tail*s


# Defining the function of s with respect to x, same head and tail as in f_x
def f_s(x, head, tail):
    return (x-head)/tail

def f_y(s, thickness):
    return thickness/0.2*(0.2969*math.sqrt(s) - 0.126*s - 0.3516*s**2 + 0.2843*s**3-0.1036*s**4)
    # Based off of NACA thickness equation

def y_tail(x, t):
    return (-1.2958*x + 0.5105*x**2)*math.sin(2*math.pi*(x/29.766+0.25*t))/10 + thickness/4
def v_tail(x, t):
    return (-1.2958*x + 0.5105*x**2)*2*math.pi*0.25*math.cos(2*math.pi*(x/29.766+0.25*t))/10 + 20

def dy_dx_tail(x, t):
    return (-1.2958*x + 0.5105*x**2)*2*math.pi/29.766*math.cos(2*math.pi*(x/29.766+0.25*t))/10 + (-1.2958 + 0.5105*x*2)*math.sin(2*math.pi*(x/29.766+0.25*t))/10
def tail_length_helper(x, t):
    return 1+dy_dx_tail(x, t)**2 #using the v_tail formula here as the derivative of

'''
MAIN LOOP
'''

dt = 0.001 #time step for movement of fish tail

L_file_names = []
L_fish_file_names = []
# list of file names that will be read by MATLAB

# changing chi for whether or not the datapoint is within the fish
for t in np.arange(0, 0.006, dt):

    file_text = ""
    fish_text = ""
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

                    #using arc length formula to find length of tail at x coordinate b starting from the base of the tail (x=0)
                    tail_length = quad(tail_length_helper, 0, x[i], args=(t))[0]

                    if (tail_length<=fish_length/2): #ensuring that the displayed tail doesn't randomly grow longer and shorter
                        y_i_bot =  y_tail(x[i], t)-f_y(s_i, thickness)-thickness/4
                    else:
                        y_i_bot = y_i

                else:
                    y_i = f_y(s_i, thickness)
                    y_i_bot = -1*f_y(s_i, thickness)


                if y[j]>y_i_bot and y[j]<y_i: #If coordinate is within the body of the fish
                    chi[i, j] = 1


                    if x[i] > 0:
                        vbodyy[i, j] = v_tail(x[i], t)

                    '''
                    BRINKMAN PENALIZATION
                    '''
                    # chi equal to one inside the body of the fish, zero outside the fish
                    Fx[i, j] = -1*chi[i, j]/nu*(vx[i, j]-vbodyx[i, j])
                    Fy[i, j] = -1*chi[i, j]/nu*(vy[i, j]-vbodyy[i, j])

                    fish_text += f"{X[i, j]} {Y[i, j]} {chi[i, j]} {vx[i, j]} {vy[i, j]} {vbodyx[i, j]} {vbodyy[i, j]} \n"

                else:
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

            if chi[i, j] == 0:
                file_text += f"{X[i, j]} {Y[i, j]} {chi[i, j]} {vx[i, j]} {vy[i, j]} {vbodyx[i, j]} {vbodyy[i, j]} \n"


    with open("dataout/nodes" +str(t)+ ".txt", "w") as file:
        file.write(file_text)
    with open("dataout/fish" +str(t)+ ".txt", "w") as file:
        file.write(fish_text)

    L_file_names.append("dataout/nodes" +str(t)+ ".txt")
    L_fish_file_names.append("dataout/fish" +str(t)+ ".txt")

L_files = []
L_files.append(L_file_names)
L_files.append(L_fish_file_names)
