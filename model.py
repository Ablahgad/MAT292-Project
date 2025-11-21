#import libraries
import numpy as np
import math

# Setting up Node object
class Node:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __str__(self):
        return f"{self.x} {self.y} {self.chi} \n"

# Create list of Node objects
L_nodes = []

#input section - real domain dimensions (mm), number of nodes along x, number of nodes along y
domain_x = 6 # supposed to be 150
domain_y = 2 # supposed to be 45.7

divisions_x = 10
divisions_y = 10

#create matrix with x and y intersections (nodal locations)
nodes = np.zeros(((divisions_x * divisions_y), 2))

#initialize loop variables
index = 0
d_x = (domain_x/(divisions_x - 1))
d_y = (domain_y/(divisions_y - 1))

#loop through all columns
for i in range(divisions_x):
    #loop through all rows
    for j in range(divisions_y):
        L_nodes.append(Node(i * d_x, j * d_y))
        L_nodes.append(Node(i * -1*d_x, j * -1*d_y))
        L_nodes.append(Node(i * d_x, j * -1*d_y))
        L_nodes.append(Node(i * -1*d_x, j * d_y))


        #Initialize pressures and temperatures ****


        x = i * d_x
        y = j * d_y
        nodes[index] = (x,y)
        index += 1



#print values in matrix
print(nodes)
print(nodes.shape)

print("file saved")
np.savetxt(r"nodes.txt", nodes, fmt="%.6f", comments='')


# Initializing parameter for
s = [0]
ds = 1/divisions_x
sum_s=0
for i in range(divisions_x):
    sum_s += ds
    s.append(sum_s)

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

# L_x = [L_nodes[i].x for i in range(len(L_nodes))]
# L_y = [L_nodes[i].y for i in range(len(L_nodes))]

scale_modifier = thickness/10

def y_tail(x, t):
    return scale_modifier*(2.7826*x - 0.1485*x**2)*math.sin(2*math.pi*(x/29.766+0.25*t))+thickness/4


dt = 0.1 #time step for movement of fish tail
# t = 0.01

L_file_names = []

for t in np.arange(0, 5, dt):

    file_text = ""

    for node in L_nodes:
        s_i = f_s(node.x, head, tail)
        if s_i < 1 and s_i > 0:
            if node.x > 0:
                y_i = y_tail(node.x, t)
                y_i_bot =  y_tail(node.x, t)-f_y(s_i, thickness)-thickness/4
            else:
                y_i = f_y(s_i, thickness)
                y_i_bot = -1*f_y(s_i, thickness)
            if node.y>y_i_bot and node.y<y_i:
                node.chi = 1
            else:
                node.chi = 0
        else:
            node.chi = 0

        file_text += str(node)

    with open("dataout/nodes" +str(t)+ ".txt", "w") as file:
        file.write(file_text)

    L_file_names.append("dataout/nodes" +str(t)+ ".txt")


    # with open("dataout/filenames.txt", "w") as file:
        # file.write(str(L_file_names))