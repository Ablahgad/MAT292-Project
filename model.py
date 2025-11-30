#import libraries
import numpy as np
import math

# Setting up Node object
class Node:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __str__(self):
        return f"{self.x} {self.y} {self.chi} {self.vx} {self.vy} \n"

#input section - real domain dimensions (mm), number of nodes along x, number of nodes along y
domain_x = 6 # supposed to be 150
domain_y = 2 # supposed to be 45.7

divisions_x = 10
divisions_y = 10

# Create list of Node objects
L_nodes=[]
for i in range(divisions_y*2):
    L=[0]*divisions_x*2
    L_nodes.append(L)

L_nodes_print=[]
for i in range(divisions_y*2):
    L=[0]*divisions_x*2
    L_nodes_print.append(L)

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

y = [d_y*i for i in range(divisions_y)]
y_alt = [i*-1 for i in y]
y.reverse()
for i in y_alt:
    y.append(i)
y.reverse()

index = 0
#loop through all columns
for x_i in x:
    #loop through all rows
    for y_i in y:
        L_nodes[index%(divisions_x*2)][index//(divisions_y*2)]=Node(x_i, y_i)
        L_nodes[index%(divisions_x*2)][index//(divisions_y*2)].chi = 1

        L_nodes[index%(divisions_x*2)][index//(divisions_y*2)].vx = 1
        L_nodes[index%(divisions_x*2)][index//(divisions_y*2)].vy = 1

        #Initialize pressures and temperatures ****


        L_nodes_print[index%(divisions_x*2)][index//(divisions_y*2)]=str(L_nodes[index%(divisions_x*2)][index//(divisions_y*2)])
        index += 1

# Initializing parameters for tail equations
head = -4.5
tail = 9
thickness = 1


# Defining the function of x with respect to the parameter s
def f_x(s, head, tail):
    return head+tail*s
    # Head is the starting point of the head at the left, head+tail is the end point of the tail
    # Ex. head = -4.5, tail = 9 means that the fish starts at -4.5 and ends at 5

# Defining the function of s with respect to x, same head and tail as in f_x
def f_s(x, head, tail):
    return (x-head)/tail

def f_y(s, thickness):
    return thickness/0.2*(0.2969*math.sqrt(s) - 0.126*s - 0.3516*s**2 + 0.2843*s**3-0.1036*s**4)
    # Thickness is the width of the thickest point of the fish, ex. thickness = 1
    # Based off of NACA thickness equation


scale_modifier = thickness/10

def y_tail(x, t):
    return scale_modifier*(2.7826*x - 0.1485*x**2)*math.sin(2*math.pi*(x/29.766+0.25*t))+thickness/4


dt = 0.5 #time step for movement of fish tail

L_file_names = []

# changing chi for whether or not the datapoint is within the fish
for t in np.arange(0, 5, dt):

    file_text = ""
    index = 0
    for xi in x:
        for yi in y:
            i = index%(divisions_x*2)
            j = index//(divisions_y*2)

            s_i = f_s(xi, head, tail)
            if s_i < 1 and s_i > 0:
                if xi > 0:
                    y_i = y_tail(xi, t)
                    y_i_bot =  y_tail(xi, t)-f_y(s_i, thickness)-thickness/4
                else:
                    y_i = f_y(s_i, thickness)
                    y_i_bot = -1*f_y(s_i, thickness)
                if yi>y_i_bot and yi<y_i:
                    L_nodes[i][j].chi = 1
                    L_nodes[i][j].vx = 0
                    L_nodes[i][j].vy = 0
                else:
                    L_nodes[i][j].chi = 0
            else:
                L_nodes[i][j].chi = 0
            index += 1

            file_text += str(L_nodes[i][j])

    with open("dataout/nodes" +str(t)+ ".txt", "w") as file:
        file.write(file_text)

    L_file_names.append("dataout/nodes" +str(t)+ ".txt")


    with open("dataout/filenames.txt", "w") as file:
        file.write(str(L_file_names))