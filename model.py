#import libraries
import numpy as np
import math

# Setting up Node object
class Node:
    def __init__(self, x, y):
        self.x = x
        self.y = y

# Create list of Node objects
L_nodes = []

#input section - real domain dimensions (mm), number of nodes along x, number of nodes along y
domain_x = 150
domain_y = 45.7

divisions_x = 100
divisions_y = 100

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

# L_x = [L_nodes[i].x for i in range(len(L_nodes))]
# L_y = [L_nodes[i].y for i in range(len(L_nodes))]

head = -4.5
tail = 9
thickness = 1

for node in L_nodes:
    s_i = f_s(node.x, head, tail)
    y_i = f_y(s_i, thickness)
    if node.y>(-1)*y_i and node.y<y_i:
        node.chi = 1
    else:
        node.chi = 0

