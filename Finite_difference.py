'''
This is the standard code for the model and motion of the fish, as well as the basic information sent to MATLAB for visualization. Each of the specific numerical methods will expand on this code to include the actual fluid dynamics expected of the interaction between the fish and water.
'''


#import libraries
from scipy.integrate import quad
import numpy as np
import math

# Setting up Node object
class Node:
    def __init__(self, x, y):
        self.x = x
        self.y = y

        self.vx = 26
        self.vy = 0

        self.vbody = 0

        self.P = 0 #Pressure to zero, every resulting pressure is relative to the initial pressure within the swimming tube

        self.chi = 0

        self.v_starx = 0
        self.v_stary = 0

        self.Fx = 0
        self.Fy = 0

    def __str__(self):
        return f"{self.x} {self.y} {self.chi} {self.vx} {self.vy} {self.vbody}\n" # Add any new variables that are important to the visualization here, MATLAB only has what is printed out as a text file from this code
'''
CONSTANTS
'''
rho = 1 #density
kv = 1 # kinematic viscosity

# nu is time penalty in Brinkman Penalization, determines how fast the velocity changes to the body velocity, the larger it is the more fluid it lets into the fish, smaller it is the faster the velocity change (recommended between 0.1 and 1)
nu = 1

#input section - real domain dimensions (mm), number of nodes along x, number of nodes along y
domain_x = 30 # full x dimension of tube is 150mm
domain_y = 10 # full y dimension of tube is 45.7

divisions_x = 20
divisions_y = 20

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
INITIALIZATION LOOP
'''
#loop through all columns
for x_i in x:
    #loop through all rows
    for y_i in y:
        i_j = str(x_i) + " "+ str(y_i)
        D_nodes[i_j]=Node(x_i, y_i)


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
    return (-1.2958*x + 0.5105*x**2)*2*math.pi*0.25*math.cos(2*math.pi*(x/29.766+0.25*t))/10

def dy_dx_tail(x, t):
    return (-1.2958*x + 0.5105*x**2)*2*math.pi/29.766*math.cos(2*math.pi*(x/29.766+0.25*t))/10 + (-1.2958 + 0.5105*x*2)*math.sin(2*math.pi*(x/29.766+0.25*t))/10
def tail_length_helper(x, t):
    return 1+dy_dx_tail(x, t)**2 #using the v_tail formula here as the derivative of

'''
MAIN LOOP
'''

dt = 0.002 #time step for movement of fish tail

L_file_names = []
L_fish_file_names = []
# list of file names that will be read by MATLAB

# changing chi for whether or not the datapoint is within the fish
for t in np.arange(0, 0.5, dt):

    file_text = ""
    fish_text = ""
    index = 0
    for i in range(1, len(x)-1):
        for j in range(1, len(y)-1):
            i_j = str(x[i]) + ' ' + str(y[j])
            i1_j = str(x[i+1]) + ' ' + str(y[j])
            i_j1 = str(x[i]) + ' ' + str(y[j+1])
            im1_j = str(x[i-1]) + ' ' + str(y[j])
            i_jm1 = str(x[i]) + ' ' + str(y[j-1])
            im1_jm1 = str(x[i-1]) + ' ' + str(y[j-1])
            i1_j1 = str(x[i+1]) + ' ' + str(y[j+1])

            '''
            Updating where the fish is
            '''

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


                if y[j]>y_i_bot and y[j]<y_i: #If
                    D_nodes[i_j].chi = 1


                    if x[i] > 0:
                        D_nodes[i_j].vbody = v_tail(x[i], t)

                    '''
                    BRINKMAN PENALIZATION
                    '''
                    # chi equal to one inside the body of the fish, zero outside the fish
                    D_nodes[i_j].Fx = -1*D_nodes[i_j].chi/nu*(D_nodes[i_j].vx-D_nodes[i_j].vbody)
                    D_nodes[i_j].Fy = -1*D_nodes[i_j].chi/nu*(D_nodes[i_j].vy-D_nodes[i_j].vbody)

                    fish_text += str(D_nodes[i_j])

                else:
                    D_nodes[i_j].chi = 0
            else:
                D_nodes[i_j].chi = 0

            '''
            FINITE DIFFERENCE ESTIMATIONS OF SLOPES
            '''
            dvx_dx = (D_nodes[i1_j].vx - D_nodes[im1_j].vx)/2/d_x
            dvx_dy = (D_nodes[i_j1].vx - D_nodes[i_jm1].vx)/2/d_y


            dvy_dx = (D_nodes[i1_j].vy - D_nodes[im1_j].vy)/2/d_x
            dvy_dy = (D_nodes[i_j1].vy - D_nodes[i_jm1].vy)/2/d_y


            # Laplacian operator
            D2vx = (D_nodes[i1_j].vx -2*D_nodes[i_j].vx + D_nodes[im1_j].vx)/d_x**2 + (D_nodes[i_j1].vx -2*D_nodes[i_j].vx + D_nodes[i_jm1].vx)/d_y**2
            D2vy = (D_nodes[i1_j].vy -2*D_nodes[i_j].vy + D_nodes[im1_j].vy)/d_x**2 + (D_nodes[i_j1].vy -2*D_nodes[i_j].vy + D_nodes[i_jm1].vy)/d_y**2

            dP_dx = (D_nodes[i1_j].P - D_nodes[im1_j].P)/2/d_x

            '''
            ACTUAL NAVIER-STOKES
            '''
            dv_dtx = -1*(D_nodes[i_j].vx*dvx_dx + D_nodes[i_j].vy*dvx_dy) - 1/rho*dP_dx + kv*(D2vx) + 1/rho*D_nodes[i_j].Fx
            dv_dty = -1*(D_nodes[i_j].vx*dvy_dx + D_nodes[i_j].vy*dvy_dy) - 1/rho*dP_dx + kv*(D2vy) + 1/rho*D_nodes[i_j].Fy


            D_nodes[i_j].v_starx = D_nodes[i_j].vx + dt*dv_dtx
            D_nodes[i_j].v_stary = D_nodes[i_j].vy + dt*dv_dty


            '''
            JACOBIAN ITERATION
            '''

            b = rho/dt*((D_nodes[i1_j].v_starx-D_nodes[im1_j].v_starx)/2/d_x + (D_nodes[i_j1].v_stary-D_nodes[i_jm1].v_stary)/2/d_y)

            for m in range(50):
                D_nodes[i_j].P = (D_nodes[i1_j].P + D_nodes[im1_j].P + D_nodes[i_j1].P + D_nodes[i_jm1].P - b*d_x**2)/4



            '''
            PRESSURE CORRECTED V
            '''
            D_nodes[i_j].vx = D_nodes[i_j].v_starx - dt/rho*(D_nodes[i1_j].P - D_nodes[im1_j].P)/2/d_x
            D_nodes[i_j].vy = D_nodes[i_j].v_stary - dt/rho*(D_nodes[i1_j].P - D_nodes[im1_j].P)/2/d_y


            file_text += str(D_nodes[i_j])


    with open("dataout/nodes" +str(t)+ ".txt", "w") as file:
        file.write(file_text)
    with open("dataout/fish" +str(t)+ ".txt", "w") as file:
        file.write(fish_text)

    L_file_names.append("dataout/nodes" +str(t)+ ".txt")
    L_fish_file_names.append("dataout/fish" +str(t)+ ".txt")

L_files = []
L_files.append(L_file_names)
L_files.append(L_fish_file_names)
