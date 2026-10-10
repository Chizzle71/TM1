import math
from numba import njit, prange, set_num_threads
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps
from matplotlib.patches import Circle, Rectangle, RegularPolygon
list(colormaps)

plt.rcParams.update({
    'font.family': 'Times New Roman',
    'font.size': 12*1.8
})
plt.rcParams['mathtext.fontset'] = 'stix'

# Cores
set_num_threads(18)

# COMSOL related
cellsx=30
cellsy=30
hn=0.01
E=-54.7*hn**4+52.06*hn**3-15.05*hn**2+4.932*hn-0.003074

# Solver settings
divisions=12
maxiterations=1e7
diagfreq=1e4
converged=1e-7*diagfreq # Currently 0.001K difference over 10000 iterations
damping=0.425
edge=0.5



filename="inputthetas1_mut1.txt"

pointsx=cellsx*(divisions-1)+1
pointsy=cellsy*(divisions-1)+1

cells_theta=np.loadtxt(filename)
theta=np.deg2rad(cells_theta)

cells_k_xx=np.cos(theta)**2+E*np.sin(theta)**2
cells_k_xy=(1-E)*np.cos(theta)*np.sin(theta)
cells_k_yx=cells_k_xy.copy()
cells_k_yy=np.sin(theta)**2+E*np.cos(theta)**2

points=np.zeros((pointsy,pointsx))
for x in range(pointsx):
    points[:,x]=400-x*100/(pointsx-1)

newpoints=points.copy()
storage=points.copy()

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

@njit(parallel=True, fastmath=True) # Floating point reassociation and commutativity on
def solve(points,newpoints,storage):
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
            print(i+1, regressionlongavg)
            if regressionlongavg<converged:
                print("Converged at iteration", i+1)
                break
        points,newpoints=newpoints,points
    return newpoints

points=solve(points,newpoints,storage)

# Cloaks
cloak_locations = [
    (0.15, 0.15),
]
cloak_radii = [
    0.045,
]

# Cold features
cold_points = [
    #(x,y)
]

cold_regions = [
    #(x,y,r)
]

# Hot features
hot_points = [
    #(x,y)
]

hot_regions = [
    #(x,y,r)
]


fig, ax = plt.subplots(figsize=(10, 8))

im = ax.imshow(
    points,
    origin='lower',
    cmap='seismic', # gray will give the default one
    vmin=300,
    vmax=400,
    extent=[0, 0.3, 0, 0.3]
)
#plt.colorbar(im, ax=ax, label="Temperature")

# Cloaks
for (cx, cy), r in zip(cloak_locations, cloak_radii):
    ax.add_patch(
        Circle((cx, cy), r,
               edgecolor="black",#####
               linewidth=2,#####
               fill=False)
    )

# Cold point features
for x, y in cold_points:
    ax.plot(x, y, marker='x', color='blue',
            markersize=10, markeredgewidth=2)

# Cold circular features
for x, y, r in cold_regions:
    ax.add_patch(
        Circle((x, y), r,
               edgecolor="blue",
               linestyle="--",
               linewidth=2,
               fill=False)
    )

# Hot point features
for x, y in hot_points:
    ax.plot(x, y, marker='x', color='red',
            markersize=10, markeredgewidth=2)

# Hot circular features
for x, y, r in hot_regions:
    ax.add_patch(
        Circle((x, y), r,
               edgecolor="red",
               linestyle="--",
               linewidth=2,
               fill=False)
    )

ax.set_xlabel("$x$ (m)")
ax.set_ylabel("$y$ (m)")
ax.set_aspect("equal")
ax.set_title(f"2D case 1")

ax.set_xticks(np.arange(0, 0.3 + 0.001, 0.06))
ax.set_yticks(np.arange(0, 0.3 + 0.001, 0.06))

plt.savefig("2D1.png", dpi=600, bbox_inches="tight")
plt.show()