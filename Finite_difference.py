#Assume each Node starts with an x, y, chi, vx, vy, Px, Py

from model import D_nodes, d_x, d_y, dt, x, y

'''
CONSTANTS
'''
rho = 1 #density
kv = 1 # kinematic viscosity

# nu is time penalty in Brinkman Penalization, determines how fast the velocity changes to the body velocity, the larger it is the more fluid it lets into the fish, smaller it is the faster the velocity change (recommended between 0.1 and 1)
nu = 0.5


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
        BRINKMAN PENALIZATION
        '''
        # chi equal to one inside the body of the fish, zero outside the fish
        L_nodes[i][j].Fx = -1*L_nodes[i][j].chi/nu*(L_nodes[i][j].vx-v_body)
        L_nodes[i][j].Fy = -1*L_nodes[i][j].chi/nu*(L_nodes[i][j].vy-v_body)

        # need info or equation for v_body --> use period of the tail***

        '''
        FINITE DIFFERENCE ESTIMATIONS OF SLOPES
        '''
        dv_dx = (L_nodes[i+1][j].vx - L_nodes[i-1][j].vx)/2/d_x
        dv_dy = (L_nodes[i][j+1].vy - L_nodes[i][j-1].vy)/2/d_y

        # Laplacian operator
        D2v = (L_nodes[i+1][j].vx -2*L_nodes[i][j].vx + L_nodes[i-1][j].vx)/d_x**2 + (L_nodes[i][j+1].vy -2*L_nodes[i][j].vy + L_nodes[i][j-1].vy)/d_y**2

        dP_dx = (L_nodes[i+1][j].Px - L_nodes[i-1][j].Px)/2/d_x
        dP_dy = (L_nodes[i][j+1].Py - L_nodes[i][j-1].Py)/2/d_y

        '''
        ACTUAL NAVIER-STOKES
        '''
        dv_dtx = -1*(L_nodes[i][j].vx*dv_dx + L_nodes[i][j].vy*dv_dy) - 1/rho*(dP_dx+dP_dy) + kv*(D2v) + 1/rho*L_nodes[i][j].Fx
        dv_dtx = -1*(L_nodes[i][j].vx*dv_dx + L_nodes[i][j].vy*dv_dy) - 1/rho*(dP_dx+dP_dy) + kv*(D2v) + 1/rho*L_nodes[i][j].Fy


        L_nodes[i][j].v_starx = L_nodes[i][j].vx + dt*dv_dtx
        L_nodes[i][j].v_stary = L_nodes[i][j].vy + dt*dv_dty


        '''
        JACOBIAN ITERATION
        '''

        b = rho/dt*((L_nodes[i+1][j].v_starx-L_nodes[i-1][j].v_starx)/2/dx + L_nodes[i][j+1].v_stary-L_nodes[i][j-1].v_stary)/2/d_y)


        L_nodes[i][j].Px = (L_nodes[i+1][j].Px + L_nodes[i-1][j].Px + L_nodes[i][j+1].Px + L_nodes[i][j-1].Px - b*d_x**2)/4
        L_nodes[i][j].Py = (L_nodes[i+1][j].Py + L_nodes[i-1][j].Py + L_nodes[i][j+1].Py + L_nodes[i][j-1].Py - b*d_x**2)/4


        '''
        PRESSURE CORRECTED V
        '''

        L_nodes[i][j].vx = L_nodes[i][j].v_starx - dt/rho*(L_nodes[i+1][j].Px - L_nodes[i-1][j].Px)/2/d_x
        L_nodes[i][j].vy = L_nodes[i][j].v_stary - dt/rho*(L_nodes[i+1][j].Py - L_nodes[i-1][j].Py)/2/d_y




