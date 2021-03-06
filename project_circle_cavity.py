#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 15:48:13 2020

@author: sixtine
"""

# problems : wall is not intirely symetric 
#           periodical BC ?
#           BC de la dérivée 

import numpy
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D



def build_up_b(rho, dt, dx, dy, u, v):
    b = numpy.empty_like(u)
    b[1:-1, 1:-1] = (rho * (1 / dt * ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx) + 
                                      (v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)) - 
                            ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx))**2 - 
                            2 * ((u[2:, 1:-1] - u[0:-2, 1:-1]) / (2 * dy) * 
                                 (v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dx))- 
                            ((v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy))**2)) 
    
    # # Periodic BC Pressure @ x = xmax
    # b[1:-1, -1] = (rho * (1 / dt * ((u[1:-1, 0] - u[1:-1,-2]) / (2 * dx) + \
    #                                 (v[2:, -1] - v[0:-2, -1]) / (2 * dy)) - \
    #                       ((u[1:-1, 0] - u[1:-1, -2]) / (2 * dx))**2 - \
    #                       2 * ((u[2:, -1] - u[0:-2, -1]) / (2 * dy) * \
    #                             (v[1:-1, 0] - v[1:-1, -2]) / (2 * dx)) - \
    #                       ((v[2:, -1] - v[0:-2, -1]) / (2 * dy))**2)) \

    # # Periodic BC Pressure @ x = xmin
    # b[1:-1, 0] = (rho * (1 / dt * ((u[1:-1, 1] - u[1:-1, -1]) / (2 * dx) + \
    #                                 (v[2:, 0] - v[0:-2, 0]) / (2 * dy)) - \
    #                       ((u[1:-1, 1] - u[1:-1, -1]) / (2 * dx))**2 - \
    #                       2 * ((u[2:, 0] - u[0:-2, 0]) / (2 * dy) * \
    #                           (v[1:-1, 1] - v[1:-1, -1]) / (2 * dx))- \
    #                       ((v[2:, 0] - v[0:-2, 0]) / (2 * dy))**2)) \
    
    return( b)

def pressure_poisson_periodic(p, dx, dy):
    pn = numpy.empty_like(p)
    
    for q in range(nit):
        pn = p.copy()
        p[1:-1, 1:-1] = (((pn[1:-1, 2:] + pn[1:-1, 0:-2]) * dy**2 +
                          (pn[2:, 1:-1] + pn[0:-2, 1:-1]) * dx**2) /
                         (2 * (dx**2 + dy**2)) -
                         dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * b[1:-1, 1:-1])
        p =  BC(i_list,j_list_top,j_list_bottom,n_y_bottom, n_y_top,n_x,p) 
        

        # # Periodic BC Pressure @ x = xmax
        # p[1:-1, -1] = (((pn[1:-1, 0] + pn[1:-1, -2])* dy**2 +
        #                 (pn[2:, -1] + pn[0:-2, -1]) * dx**2) /
        #                (2 * (dx**2 + dy**2)) -
        #                dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * b[1:-1, -1])

        # # Periodic BC Pressure @ x = xmin
        # p[1:-1, 0] = (((pn[1:-1, 1] + pn[1:-1, -1])* dy**2 +
        #                (pn[2:, 0] + pn[0:-2, 0]) * dx**2) /
        #               (2 * (dx**2 + dy**2)) -
        #               dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * b[1:-1, 0])
        
        # # Wall boundary conditions, pressure
        # p[-1, :] =p[-2, :]  # dp/dy = 0 at y = 2
        # p[0, :] = p[1, :]  # dp/dy = 0 at y = 0
        
        # Circular Boundary condition
        #p = BC(p,R,ymax,ymin,xmax,xmin,dx,dy,n) 
     
    return p

def geometry (A ,lamb,ymax,ymin):
    
    x_list = x.copy()
    y_list_top = A * numpy.sin(numpy.pi * x_list / lamb) + (ymax - A)
    y_list_bottom = A * numpy.sin(numpy.pi * x_list / lamb + numpy.pi) + (ymin + A)
        
    for n in range(nx - 1):
        i_list[n] = int(x_list[n]/dx)
        j_list_top[n] = int(y_list_top[n]/dy)
        j_list_bottom[n]= int(y_list_bottom[n]/dy)
        
    #segment lenght       
    l_top = numpy.sqrt((x_list[1:]-x_list[:-1])**2 + (y_list_top[1:]-y_list_top[:-1])**2)
    l_bottom = numpy.sqrt((x_list[1:]-x_list[:-1])**2 + (y_list_bottom[1:]-y_list_bottom[:-1])**2)

        
    #normal vectors 
    n_x = (x_list[:-1] - x_list[1:])/l_top
    n_y_top = (y_list_top[:-1] - y_list_top[1:])/l_top
    n_y_bottom = (y_list_bottom[:-1] - y_list_bottom[1:])/l_bottom
    
    return x_list,y_list_top,y_list_bottom,i_list,j_list_top,j_list_bottom,n_y_bottom, n_y_top,n_x
        
def BC(i_list,j_list_top,j_list_bottom,n_y_bottom, n_y_top,n_x,u): # The fonction input =0 at the wall, and outside 
    un = numpy.empty_like(u)
    
    for q in range(nx-1):
        un = u.copy()
        #print(q)
        u[int(j_list_top[q]):,int(i_list[q])] = 0 
        u[:int(j_list_bottom[q]),int(i_list[q])] = 0 
        #u[int(j_list_top[q])   ,int(i_list[q])] = 1/((1/n_x) + (1/n_y_top))    * ((1/n_x*un[int(j_list_top[q]),   int(i_list[q])-1 ]) + (1/n_y_top   *un[int(j_list_top[q]-1),int(i_list[q])]))
        #u[int(j_list_bottom[q]),int(i_list[q])] = 1/((1/n_x) + (1/n_y_bottom)) * ((1/n_x*un[int(j_list_bottom[q]),int(i_list[q])-1 ]) + (1/n_y_bottom*un[int(j_list_top[q] - 1),int(i_list[q])]))                                                                   
    
<<<<<<< Updated upstream
    for q in range(n): #top circular boundaries
        ya = ymax - R
        xa = (2*q + 1) *  R + xmin 
        while i < (xmin + 2 * (q +1) * R) / dx :
            u[int((numpy.sqrt(R**2 - (i*dx - xa )**2 ) + ya)/dy) : ,int(i) ] = 0
            i = i + 1 
            
    for q in range(n): #bottom circular boundaries               
        ya = ymin + R 
        xa = (2*q + 1) * R + xmin 
        while j < (xmin + 2 * (q + 1) * R) / dx:
            u[: int(( - numpy.sqrt(R**2 - (j*dx - xa )**2 ) + ya)/dy) + 1, int(j) ] = 0
            j = j + 1 
            
    u[ int((ymax-R)/dy) : , -1] = 0 
    u[: int((ymin+R)/dy) , -1] = 0 
=======
>>>>>>> Stashed changes
    return(u)


##variable declarations
nx = 50
ny = 50
nt = 10
nit = 50 
c = 1
n= 3 # number of half circle
R = 3 # radius of the curve of the wall 
xmin = 0
xmax = n*2*R
ymin =0 
ymax = 7
dx = (xmax - xmin ) / (nx - 1)
dy = (ymax - ymin) / (ny - 1)
x = numpy.linspace(xmin, xmax, nx)
y = numpy.linspace(ymin, ymax, ny)
X, Y = numpy.meshgrid(x, y)
A = 1 # sinusoidal amplitude of the wall
lamb = 2 # wavelenght of sinusoidal wall 

##physical variables
rho = 1
nu = .1
F = 1
dt = .01

#initial conditions 
u = numpy.zeros((ny, nx))
un = numpy.zeros((ny, nx))

v = numpy.zeros((ny, nx))
vn = numpy.zeros((ny, nx))

p = numpy.ones((ny, nx))
pn = numpy.ones((ny, nx))

b = numpy.zeros((ny, nx))

i_list = numpy.zeros(nx)
j_list_top = numpy.zeros(ny)
j_list_bottom = numpy.zeros(ny)

udiff = 1
stepcount = 0

while udiff > .05:
    un = u.copy()
    vn = v.copy()
    b = build_up_b(rho, dt, dx, dy, u, v)
    p = pressure_poisson_periodic(p, dx, dy)

    u[1:-1, 1:-1] = (un[1:-1, 1:-1] -
                     un[1:-1, 1:-1] * dt / dx * 
                    (un[1:-1, 1:-1] - un[1:-1, 0:-2]) -
                     vn[1:-1, 1:-1] * dt / dy * 
                    (un[1:-1, 1:-1] - un[0:-2, 1:-1]) -
                     dt / (2 * rho * dx) * 
                    (p[1:-1, 2:] - p[1:-1, 0:-2]) +
                     nu * (dt / dx**2 * 
                    (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) +
                     dt / dy**2 * 
                    (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1])) + 
                     F * dt)

    v[1:-1, 1:-1] = (vn[1:-1, 1:-1] -
                     un[1:-1, 1:-1] * dt / dx * 
                    (vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) -
                     vn[1:-1, 1:-1] * dt / dy * 
                    (vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) -
                     dt / (2 * rho * dy) * 
                    (p[2:, 1:-1] - p[0:-2, 1:-1]) +
                     nu * (dt / dx**2 *
                    (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2]) +
                     dt / dy**2 * 
                    (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1])))

    # # Periodic BC u @ x = xmax     
    # u[1:-1, -1] = (un[1:-1, -1] - un[1:-1, -1] * dt / dx * 
    #               (un[1:-1, -1] - un[1:-1, -2]) -
    #                vn[1:-1, -1] * dt / dy * 
    #               (un[1:-1, -1] - un[0:-2, -1]) -
    #                dt / (2 * rho * dx) *
    #               (p[1:-1, 0] - p[1:-1, -2]) + 
    #                nu * (dt / dx**2 * 
    #               (un[1:-1, 0] - 2 * un[1:-1,-1] + un[1:-1, -2]) +
    #                dt / dy**2 * 
    #               (un[2:, -1] - 2 * un[1:-1, -1] + un[0:-2, -1])) + F * dt)

    # # Periodic BC u @ x = xmin
    # u[1:-1, 0] = (un[1:-1, 0] - un[1:-1, 0] * dt / dx *
    #              (un[1:-1, 0] - un[1:-1, -1]) -
    #               vn[1:-1, 0] * dt / dy * 
    #              (un[1:-1, 0] - un[0:-2, 0]) - 
    #               dt / (2 * rho * dx) * 
    #              (p[1:-1, 1] - p[1:-1, -1]) + 
    #               nu * (dt / dx**2 * 
    #              (un[1:-1, 1] - 2 * un[1:-1, 0] + un[1:-1, -1]) +
    #               dt / dy**2 *
    #              (un[2:, 0] - 2 * un[1:-1, 0] + un[0:-2, 0])) + F * dt)

    # # Periodic BC v @ x = xmax
    # v[1:-1, -1] = (vn[1:-1, -1] - un[1:-1, -1] * dt / dx *
    #               (vn[1:-1, -1] - vn[1:-1, -2]) - 
    #                vn[1:-1, -1] * dt / dy *
    #               (vn[1:-1, -1] - vn[0:-2, -1]) -
    #                dt / (2 * rho * dy) * 
    #               (p[2:, -1] - p[0:-2, -1]) +
    #                nu * (dt / dx**2 *
    #               (vn[1:-1, 0] - 2 * vn[1:-1, -1] + vn[1:-1, -2]) +
    #                dt / dy**2 *
    #               (vn[2:, -1] - 2 * vn[1:-1, -1] + vn[0:-2, -1])))

    # # Periodic BC v @ x = xmin
    # v[1:-1, 0] = (vn[1:-1, 0] - un[1:-1, 0] * dt / dx *
    #              (vn[1:-1, 0] - vn[1:-1, -1]) -
    #               vn[1:-1, 0] * dt / dy *
    #              (vn[1:-1, 0] - vn[0:-2, 0]) -
    #               dt / (2 * rho * dy) * 
    #              (p[2:, 0] - p[0:-2, 0]) +
    #               nu * (dt / dx**2 * 
    #              (vn[1:-1, 1] - 2 * vn[1:-1, 0] + vn[1:-1, -1]) +
                 #  dt / dy**2 * 
                 # (vn[2:, 0] - 2 * vn[1:-1, 0] + vn[0:-2, 0])))


    # # Wall BC: u,v = 0 @ y = ymin,ymax
    # u[0, :] = 0
    # u[-1, :] = 0
    # v[0, :] = 0
    # v[-1, :]=0
    
    # Circular boundary condition
    x_list,y_list_top,y_list_bottom,i_list,j_list_top,j_list_bottom,n_y_bottom, n_y_top,n_x = geometry (A ,lamb,ymax,ymin)
    u =  BC(i_list,j_list_top,j_list_bottom,n_y_bottom, n_y_top,n_x,u) 
    v =  BC(i_list,j_list_top,j_list_bottom,n_y_bottom, n_y_top,n_x,v)
    
    udiff = (numpy.sum(u) - numpy.sum(un)) / numpy.sum(u)
    stepcount += 1

fig = pyplot.figure(figsize = (11,7), dpi=100)
pyplot.quiver(X[::3, ::3], Y[::3, ::3], u[::3, ::3], v[::3, ::3]);

fig = pyplot.figure(figsize = (11,7), dpi=100)
pyplot.quiver(X, Y, u, v);

fig = pyplot.figure(figsize=(12,6), dpi=100)
fig.add_subplot(111,aspect='equal')
pyplot.fill(x,y,fill=False,lw=3)
pyplot.contourf(X,Y,p,15,alpha=0.5)
pyplot.colorbar()
pyplot.quiver(X,Y,u,v)
pyplot.title('Velocity & pressure through constriction')
pyplot.show()

#%% Test BC

u = numpy.ones((ny,nx))
u = BC(u,R,ymax,ymin,xmax,xmin,dx,dy,n)

v = numpy.ones((ny,nx))
v = BC(v,R,ymax,ymin,xmax,xmin,dx,dy,n)

fig = pyplot.figure(figsize = (11,7), dpi=100)
pyplot.quiver(X, Y, u, v);

#%% Test 
<<<<<<< Updated upstream
m = [0,1,2,4]
=======

u = numpy.linspace(1, 3, 3)
u[int(2.0)]
>>>>>>> Stashed changes

