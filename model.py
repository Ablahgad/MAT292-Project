#import libraries
import numpy as np

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
        x = i * d_x
        y = j * d_y
        nodes[index] = (x,y)
        index += 1

#print values in matrix
print(nodes)
print(nodes.shape)

print("file saved")
np.savetxt(r"C:\Users\22cia\Documents\nodes.txt", nodes, fmt="%.6f", header="x y", comments='')
