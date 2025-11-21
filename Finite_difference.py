#Assume each Node starts with an x, y, chi, vx, vy, Px, Py

from model import L_nodes

'''
CONSTANTS
'''
dt = 1 #time step
rho = 1 #density
kv = 1 # kinematic viscosity

dx = 1
dy = 1


for i in range(len(L_nodes)):


    '''
    BRINKMAN PENALIZATION
    '''
    # chi equal to one inside the body of the fish, zero outside the fish
    # nu is time penalty, determines how fast the velocity changes to the body velocity, the larger it is the more fluid it lets into the fish, smaller it is the faster the velocity change (recommended between 0.1 and 1)

        node.Fx = -1*chi/nu*(node.vx-v_body)
        node.Fy = -1*chi/nu*(node.vy-v_body)

    # need info or equation for v_body ***

    '''
    FINITE DIFFERENCE ESTIMATIONS OF SLOPES
    '''
    dv_dx = (v[0][i+1][j] - v[0][i-1][j])/2/dx
    dv_dy = (v[1][i][j+1] - v[1][i][j-1])/2/dy

    # Laplacian operator
    D2v = (v[0][i+1][j] -2*v[0][i][j] + v[0][i-1][j])/dx**2 + (v[1][i][j+1] -2*v[1][i][j] + v[1][i][j-1])/dy**2

    dP_dx = (P[0][i+1][j] - P[0][i-1][j])/2/dx
    dP_dy = (P[1][i][j+1] - P[1][i][j-1])/2/dy

    '''
    ACTUAL NAVIER-STOKES
    '''
    dv_dt = -1*(v[0]*dv_dx + v[1]*dv_dy) - 1/rho*(dP_dx+dP_dy) + kv*(D2v) + 1/rho*F

    v_star = v + dt*dv_dt #returns both x and y components


    '''
    JACOBIAN ITERATION
    '''

    b[i][j] = rho/dt*((v_star[0][i+1][j]-v_star[0][i-1][j])/2/dx + v_star[1][i][j+1]-v_star[1][i][j-1])/2/dy)

    P[i, j] = (P[i+1][j] + P[i-1][j] + P[i][j+1] + P[i][j-1] - b[i][j]*dx**2)/4


    '''
    PRESSURE CORRECTED V
    '''

    v_next = v_star - dt/rho*(P[i+1, j] - P[i-1, j])/2/dx


