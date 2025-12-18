'''
This is the standard code for the model and motion of the fish, as well as the basic information sent to MATLAB for visualization. Each of the specific numerical methods will expand on this code to include the actual fluid dynamics expected of the interaction between the fish and water.
'''


#import libraries
import numpy as np
import math

# Setting up Node object
class Node:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __str__(self):
        return f"{self.x} {self.y} {self.chi} {self.vx} {self.vy} \n" # Add any new variables that are important to the visualization here, MATLAB only has what is printed out as a text file from this code

#input section - real domain dimensions (mm), number of nodes along x, number of nodes along y
domain_x = 50 # full x dimension of tube is 150mm
domain_y = 10 # full y dimension of tube is 45.7

divisions_x = 30
divisions_y = 30

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
        i_dict = str(x_i) + " "+ str(y_i)
        D_nodes[i_dict]=Node(x_i, y_i)

        D_nodes[i_dict].chi = 1

        D_nodes[i_dict].vx = 1
        D_nodes[i_dict].vy = 1

# Initializing parameters for tail equations
head = -35/2
tail = 35
    # Head is the starting point of the head at the left, head+tail is the end point of the tail
    # Ex. head = -4.5, tail = 9 means that the fish starts at -4.5 and ends at 5
thickness = 3
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


scale_modifier = thickness

def y_tail(x, t):
    return (2.7826*x - 0.1485*x**2)*math.sin(2*math.pi*(x/29.766+0.25*t))+thickness/4

def v_tail(x, t):
    return (2.7826*x - 0.1485*x**2)*2*math.pi*0.25*math.cos(2*math.pi*(x/29.766+0.25*t))+thickness/4

dt = 4 #time step for movement of fish tail

L_file_names = []
# list of file names that will be read by MATLAB

# changing chi for whether or not the datapoint is within the fish
for t in np.arange(0, 5, dt):

    file_text = ""
    index = 0
    for xi in x:
        for yi in y:
            i_dict = str(xi) + " "+ str(yi)

            D_nodes[i_dict].vx = 1
            D_nodes[i_dict].vy = 1

            s_i = f_s(xi, head, tail)
            if s_i < 1 and s_i > 0:
                if xi > 0:
                    y_i = y_tail(xi, t)
                    y_i_bot =  y_tail(xi, t)-f_y(s_i, thickness)-thickness/4
                else:
                    y_i = f_y(s_i, thickness)
                    y_i_bot = -1*f_y(s_i, thickness)
                if yi>y_i_bot and yi<y_i:
                    D_nodes[i_dict].chi = 1
                    D_nodes[i_dict].vx = 0
                    D_nodes[i_dict].vy = 0
                    if xi > 0:
                        D_nodes[i_dict].vbody = v_tail(xi, t)
                else:
                    D_nodes[i_dict].chi = 0
            else:
                D_nodes[i_dict].chi = 0
            index += 1

            file_text += str(D_nodes[i_dict])

    with open("dataout/nodes" +str(t)+ ".txt", "w") as file:
        file.write(file_text)

    L_file_names.append("dataout/nodes" +str(t)+ ".txt")


    with open("dataout/filenames.txt", "w") as file:
        file.write(str(L_file_names))