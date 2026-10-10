# Imports
import math
from numba import njit, prange, set_num_threads
import numpy as np
import matplotlib.pyplot as plt
import cma
import pickle
import os



# Mac
set_num_threads(6)

# COMSOL related
cellsx=20
cellsy=20
hn=0.04
E=-54.7*hn**4+52.06*hn**3-15.05*hn**2+4.932*hn-0.003074

# Solver settings
divisions=7
maxiterations=1e7
diagfreq=1e4
converged=1e-6*diagfreq # currently n K difference over i iterations
damping=0.425
edge=0.5
mutationiters=5



filename="inputthetas1306_optimised4.txt"
outputfilename="inputthetas1306_optimised4.txt"

pointsx=cellsx*(divisions-1)+1
pointsy=cellsy*(divisions-1)+1

xa_lookup=np.empty(pointsx,np.int32)
xb_lookup=np.empty(pointsx,np.int32)
xc_lookup=np.empty(pointsx,np.int32)
xd_lookup=np.empty(pointsx,np.int32)
for x in range(pointsx):
    xa_lookup[x]=math.ceil(x/(divisions-1))-1
    xb_lookup[x]=math.floor(x/(divisions-1))
    xc_lookup[x]=xa_lookup[x] # ------------ could optimise a bit
    xd_lookup[x]=xb_lookup[x]

ya_lookup=np.empty(pointsy,np.int32)
yb_lookup=np.empty(pointsy,np.int32)
yc_lookup=np.empty(pointsy,np.int32)
yd_lookup=np.empty(pointsy,np.int32)
for y in range(pointsy):
    ya_lookup[y]=math.floor(y/(divisions-1))
    yb_lookup[y]=ya_lookup[y]
    yc_lookup[y]=math.ceil(y/(divisions-1))-1
    yd_lookup[y]=yc_lookup[y]



@njit(parallel=True, fastmath=True)
def solve(
    points, newpoints, storage,
    cells_k_xx, cells_k_xy, cells_k_yx, cells_k_yy,
    xa_lookup, xb_lookup, xc_lookup, xd_lookup,
    ya_lookup, yb_lookup, yc_lookup, yd_lookup):
    for i in range(maxiterations):
        regressionlong=0.0
        for y in prange(pointsy):
            for x in range(pointsx):
                if x==0 or x==pointsx-1:
                    continue
                # Bottom stencil
                if y==0:
                    xa=xa_lookup[x]; xb=xb_lookup[x]
                    ya=ya_lookup[y]; yb=yb_lookup[y]
                    pc=points[y,x]
                    pu=points[y+1,x]; pl=points[y,x-1]; pr=points[y,x+1]
                    pur=points[y+1,x+1]; pul=points[y+1,x-1]
                    newpoints[y][x] = pc + damping*(0.5*(pu-pc)*(cells_k_yy[ya,xa]+cells_k_yy[yb,xb]) + 0.125*(pur+pr-pul-pl)*(cells_k_yx[ya,xa]+cells_k_yx[yb,xb]) + edge*((pl-pc)*(cells_k_xx[ya,xa]) + 0.5*(pl+pc-pul-pu)*(cells_k_xy[ya,xa]) + (pr-pc)*(cells_k_xx[yb,xb]) + 0.5*(pu+pur-pc-pr)*(cells_k_xy[yb,xb])))
                # Top stencil
                elif y==pointsy-1:
                    xc=xc_lookup[x]; xd=xd_lookup[x]
                    yc=yc_lookup[y]; yd=yd_lookup[y]
                    pc=points[y,x]
                    pd=points[y-1,x]; pl=points[y,x-1]; pr=points[y,x+1]
                    pdr=points[y-1,x+1]; pdl=points[y-1,x-1]
                    newpoints[y][x] = pc + damping*(((pl-pc)*(cells_k_xx[yc,xc]) + 0.5*(pdl+pd-pl-pc)*(cells_k_xy[yc,xc]) + (pr-pc)*(cells_k_xx[yd,xd]) + 0.5*(pc+pr-pd-pdr)*(cells_k_xy[yd,xd]))*edge + 0.5*(pd-pc)*(cells_k_yy[yc,xc]+cells_k_yy[yd,xd]) + 0.125*(pl+pdl-pr-pdr)*(cells_k_yx[yc,xc]+cells_k_yx[yd,xd]))
                # General stencil
                else:
                    xa=xa_lookup[x]; xb=xb_lookup[x]; xc=xc_lookup[x]; xd=xd_lookup[x]
                    ya=ya_lookup[y]; yb=yb_lookup[y]; yc=yc_lookup[y]; yd=yd_lookup[y]
                    pc=points[y,x]
                    pu=points[y+1,x]; pd=points[y-1,x]; pl=points[y,x-1]; pr=points[y,x+1]
                    pur=points[y+1,x+1]; pul=points[y+1,x-1]; pdr=points[y-1,x+1]; pdl=points[y-1,x-1]
                    newpoints[y][x] = pc + damping*(0.5*(pu-pc)*(cells_k_yy[ya,xa]+cells_k_yy[yb,xb]) + 0.125*(pur+pr-pul-pl)*(cells_k_yx[ya,xa]+cells_k_yx[yb,xb]) + 0.5*(pl-pc)*(cells_k_xx[ya,xa]+cells_k_xx[yc,xc]) + 0.125*(pdl+pd-pul-pu)*(cells_k_xy[ya,xa]+cells_k_xy[yc,xc]) + 0.5*(pr-pc)*(cells_k_xx[yb,xb]+cells_k_xx[yd,xd]) + 0.125*(pu+pur-pd-pdr)*(cells_k_xy[yb,xb]+cells_k_xy[yd,xd]) + 0.5*(pd-pc)*(cells_k_yy[yc,xc]+cells_k_yy[yd,xd]) + 0.125*(pl+pdl-pr-pdr)*(cells_k_yx[yc,xc]+cells_k_yx[yd,xd]))
                # Diagnostics
                if i==0 or (i+1)%diagfreq==0:
                    regressionlong+=abs(newpoints[y,x]-storage[y,x])
                    storage[y,x]=newpoints[y,x]
        if (i+1)%diagfreq==0:
            regressionlongavg = regressionlong/(pointsx*pointsy)
            #print(i+1, regressionlongavg)
            if regressionlongavg<converged:
                print("\nConverged at iteration", i+1, regressionlongavg)
                break
        points,newpoints=newpoints,points
    return points



def run_simulation(theta):
    theta = np.deg2rad(theta)

    cells_k_xx=np.cos(theta)**2+E*np.sin(theta)**2
    cells_k_xy=(1-E)*np.cos(theta)*np.sin(theta)
    cells_k_yx=cells_k_xy.copy()
    cells_k_yy=np.sin(theta)**2+E*np.cos(theta)**2

    points=np.zeros((pointsy,pointsx))
    for x in range(pointsx):
        points[:,x]=400-x*100/(pointsx-1)

    newpoints=points.copy()
    storage=points.copy()

    points = solve(points, newpoints, storage,
                   cells_k_xx, cells_k_xy, cells_k_yx, cells_k_yy,
                   xa_lookup, xb_lookup, xc_lookup, xd_lookup,
                   ya_lookup, yb_lookup, yc_lookup, yd_lookup)
    
    return points



def objective(theta_flat):

    theta = theta_flat.reshape((cellsy, cellsx))

    points = run_simulation(theta)

    cy = points.shape[0] // 2
    cx = points.shape[1] // 2

    return points[cy, cx]



def optimise():

    if os.path.exists("cma_state.pkl"):

        print("Continuing previous optimisation...")

        with open("cma_state.pkl", "rb") as f:
            es = pickle.load(f)

    else:

        print("Starting new optimisation...")

        theta0 = np.loadtxt(filename).flatten()

        es = cma.CMAEvolutionStrategy(
            theta0,
            5.0,
            {
                "popsize": 10,
                "bounds": [-360, 360]
            }
        )


    for generation in range(mutationiters):

        solutions = es.ask()

        scores = [objective(theta) for theta in solutions]

        es.tell(solutions, scores)

        print(
            f"Generation {generation+1}: "
            f"Best = {es.result.fbest:.6f}"
        )


    best_theta = es.result.xbest.reshape((cellsy,cellsx))

    np.savetxt(outputfilename, best_theta, fmt="%.6f")


    with open("cma_state.pkl", "wb") as f:
        pickle.dump(es, f)
        


if __name__ == "__main__":
    optimise()

#this is a diffuser optimiser for the centre point