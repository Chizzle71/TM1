# Imports
import math
from numba import njit, prange, set_num_threads
import numpy as np
import matplotlib.pyplot as plt
import cma
import pickle
import os



# Temperature requirements

# Features
ftw=1 # Weighting for objective function

ft1on=True # True = feature is made, False = no feature
ft1positive=True # True = increase temperature, False = decrease temperature ####edited
ft1point=False # True = point, False = circular area
ft1location=[0.15,0.15] # Centre (x,y) of low T point from bottom left corner of domain in metres
ft1radius=0.045 # radius of circle to decrease temperature of

ft2on=False
ft2positive=False
ft2point=False
ft2location=[0.15,0.04]
ft2radius=0.02

ft3on=False
ft3positive=True
ft3point=True
ft3location=[0,0]
ft3radius=0

ftgon=True
ftgexcludep=[0.15,0.15]
ftgexcluder=0.045

# Cloaks
clw=5 # Weighting for objective functin

cloak1location=[0.15,0.15] # Centre (x,y) of cloak from bottom left corner of domain in metres
cloak1radius=0.045 # Radius of circular cloaked region in metres

cloak2location=[0,0]
cloak2radius=0

cloak3location=[0,0]
cloak3radius=0



# Cores
set_num_threads(18)

# COMSOL related
cellsx=30
cellsy=30
side=0.01 # metres
thickness=0.0001 # metres
hn=thickness/side
E=-54.7*hn**4+52.06*hn**3-15.05*hn**2+4.932*hn-0.003074

# Solver settings
divisions=7
maxiterations=1e7
diagfreq=1e3
converged=1e-5*diagfreq # currently n K difference over i iterations
damping=0.425
edge=0.5
mutationiters=100



filename="inputthetas1_mut1.txt"
outputfilename="inputthetas1_mut1.txt"

pointsx=cellsx*(divisions-1)+1
pointsy=cellsy*(divisions-1)+1

xa_lookup=np.empty(pointsx,np.int32)
xb_lookup=np.empty(pointsx,np.int32)
xc_lookup=np.empty(pointsx,np.int32)
xd_lookup=np.empty(pointsx,np.int32)
for x in range(pointsx):
    xa_lookup[x]=math.ceil(x/(divisions-1))-1
    xb_lookup[x]=math.floor(x/(divisions-1))
    xc_lookup[x]=xa_lookup[x]
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

    totaltemp1=0 # select all points in cloaked circle and average their temps
    npoints1=0
    var1=0 # iterate through all points in cloaked circle and just do abs(temp-average) for each
    if cloak1radius != 0:
        for y in range(pointsy):
            for x in range(pointsx):
                if (cloak1location[0]-x*side/(divisions-1))**2+(cloak1location[1]-y*side/(divisions-1))**2 <= cloak1radius**2:
                    totaltemp1 += points[y,x]
                    npoints1 += 1
        for y in range(pointsy):
            for x in range(pointsx):
                if (cloak1location[0]-x*side/(divisions-1))**2+(cloak1location[1]-y*side/(divisions-1))**2 <= cloak1radius**2:
                    if npoints1 != 0:
                        var1 += abs(points[y,x]-totaltemp1/npoints1)/npoints1
        if npoints1 == 0:
            print("cloak 1 too small to reach a grid point --> Increase number of divisions or increase cloak size")

    totaltemp2=0
    npoints2=0
    var2=0
    if cloak2radius != 0:
        for y in range(pointsy):
            for x in range(pointsx):
                if (cloak2location[0]-x*side/(divisions-1))**2+(cloak2location[1]-y*side/(divisions-1))**2 <= cloak2radius**2:
                    totaltemp2 += points[y,x]
                    npoints2 += 1
        for y in range(pointsy):
            for x in range(pointsx):
                if (cloak2location[0]-x*side/(divisions-1))**2+(cloak2location[1]-y*side/(divisions-1))**2 <= cloak2radius**2:
                    if npoints2 != 0:
                        var2 += abs(points[y,x]-totaltemp2/npoints2)/npoints2
        if npoints2 == 0:
            print("cloak 2 too small to reach a grid point --> Increase number of divisions or increase cloak size")

    totaltemp3=0
    npoints3=0
    var3=0
    if cloak3radius != 0:
        for y in range(pointsy):
            for x in range(pointsx):
                if (cloak3location[0]-x*side/(divisions-1))**2+(cloak3location[1]-y*side/(divisions-1))**2 <= cloak3radius**2:
                    totaltemp3 += points[y,x]
                    npoints3 += 1
        for y in range(pointsy):
            for x in range(pointsx):
                if (cloak3location[0]-x*side/(divisions-1))**2+(cloak3location[1]-y*side/(divisions-1))**2 <= cloak3radius**2:
                    if npoints3 != 0:
                        var3 += abs(points[y,x]-totaltemp3/npoints3)/npoints3
        if npoints3 == 0:
            print("cloak 3 too small to reach a grid point --> Increase number of divisions or increase cloak size")

    var4=0
    ft1count=0
    if ft1on and ft1point:
        if ft1positive:
            var4 = abs(350-points[round(ft1location[1]*(divisions-1)/side),round(ft1location[0]*(divisions-1)/side)]) #points[y,x]
        else:
            var4 = abs(points[round(ft1location[1]*(divisions-1)/side),round(ft1location[0]*(divisions-1)/side)]-300)
    elif ft1on and not ft1point:
        if ft1positive:
            for y in range(pointsy):
                for x in range(pointsx):
                    if (ft1location[0]-x*side/(divisions-1))**2+(ft1location[1]-y*side/(divisions-1))**2 <= ft1radius**2:
                        ft1count += 1
                        var4 += abs(350-points[y,x])#################
        else:
            for y in range(pointsy):
                for x in range(pointsx):
                    if (ft1location[0]-x*side/(divisions-1))**2+(ft1location[1]-y*side/(divisions-1))**2 <= ft1radius**2:
                        ft1count += 1
                        var4 += abs(points[y,x]-300)
        if ft1count != 0:
            var4 = var4/ft1count
        else:
            print("Feature 1 too small to reach a grid point --> Increase number of divisions or increase feature radius")

    var5=0
    ft2count=0
    if ft2on and ft2point:
        if ft2positive:
            var5 = abs(400-points[round(ft2location[1]*(divisions-1)/side),round(ft2location[0]*(divisions-1)/side)])
        else:
            var5 = abs(points[round(ft2location[1]*(divisions-1)/side),round(ft2location[0]*(divisions-1)/side)]-300)
    elif ft2on and not ft2point:
        if ft2positive:
            for y in range(pointsy):
                for x in range(pointsx):
                    if (ft2location[0]-x*side/(divisions-1))**2+(ft2location[1]-y*side/(divisions-1))**2 <= ft2radius**2:
                        ft2count += 1
                        var5 += abs(400-points[y,x])
        else:
            for y in range(pointsy):
                for x in range(pointsx):
                    if (ft2location[0]-x*side/(divisions-1))**2+(ft2location[1]-y*side/(divisions-1))**2 <= ft2radius**2:
                        ft2count += 1
                        var5 += abs(points[y,x]-300)
        if ft2count != 0:
            var5 = var5/ft2count
        else:
            print("Feature 2 too small to reach a grid point --> Increase number of divisions or increase feature radius")

    var6=0
    ft3count=0
    if ft3on and ft3point:
        if ft3positive:
            var6 = abs(400-points[round(ft3location[1]*(divisions-1)/side),round(ft3location[0]*(divisions-1)/side)])
        else:
            var6 = abs(points[round(ft3location[1]*(divisions-1)/side),round(ft3location[0]*(divisions-1)/side)]-300)
    elif ft3on and not ft3point:
        if ft3positive:
            for y in range(pointsy):
                for x in range(pointsx):
                    if (ft3location[0]-x*side/(divisions-1))**2+(ft3location[1]-y*side/(divisions-1))**2 <= ft3radius**2:
                        ft3count += 1
                        var6 += abs(400-points[y,x])
        else:
            for y in range(pointsy):
                for x in range(pointsx):
                    if (ft3location[0]-x*side/(divisions-1))**2+(ft3location[1]-y*side/(divisions-1))**2 <= ft3radius**2:
                        ft3count += 1
                        var6 += abs(points[y,x]-300)
        if ft3count != 0:
            var6 = var6/ft3count
        else:
            print("Feature 3 too small to reach a grid point --> Increase number of divisions or increase feature radius")

    var7=0
    ftgcount=0
    if ftgon:
        for y in range(pointsy):
            for x in range(pointsx):
                if (ftgexcludep[0]-x*side/(divisions-1))**2+(ftgexcludep[1]-y*side/(divisions-1))**2 > ftgexcluder**2:
                    ftgcount += 1
                    var7 += abs(points[y,x]-(400-x*100/(pointsx-1)))
        if ftgcount != 0:
            var7 = var7/ftgcount


    vartotal=clw*(var1+var2+var3)+ftw*(2*var4+var5+var6+var7)#########

    return vartotal



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
            5,
            {
                "popsize": 12
            }
        )

    for generation in range(mutationiters):

        solutions = es.ask()
        scores = [objective(theta) for theta in solutions]
        es.tell(solutions, scores)

        print(
            f"Generation {es.countiter} | "
            f"Best {es.result.fbest:.6f} | "
            f"sigma {es.sigma:.4f}"
        )

    best_theta = es.result.xbest.reshape((cellsy,cellsx))
    np.savetxt(outputfilename, best_theta, fmt="%.6f")

    with open("cma_state.pkl", "wb") as f:
        pickle.dump(es, f)



if __name__ == "__main__":
    optimise()