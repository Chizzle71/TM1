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
ft1positive=True # True = increase temperature, False = decrease temperature
ft1point=False # True = point, False = sphere
ft1location=[0.8,0.8,0.8] # Centre (x,y,z) of low T point from bottom left near corner of domain in metres
ft1radius=0.24 # radius of sphere to decrease temperature of

ft2on=False
ft2positive=False
ft2point=False
ft2location=[0.48,0.48,0.48]
ft2radius=0.16

ft3on=False
ft3positive=True
ft3point=True
ft3location=[0,0,0]
ft3radius=0

ftgon=True
ftgexcludep=[0.8,0.8,0.8]
ftgexcluder=0.24

# Cloaks
clw=5 # Weighting for objective functin

cloak1location=[0.8,0.8,0.8] # Centre (x,y,z) of cloak from bottom left near corner of domain in metres
cloak1radius=0.24 # Radius of spherical cloaked region in metres

cloak2location=[0,0,0]
cloak2radius=0

cloak3location=[0,0,0]
cloak3radius=0

 

# Cores
set_num_threads(18)

# COMSOL related
cellsx=20
cellsy=20
cellsz=20
side=0.08
rns=0.015 # 0.015 for machine learning comparison but thats too anisotropic tbh - using second example from paper now
E = -43520*rns*rns*rns*rns+11610*rns*rns*rns-1166*rns*rns+57.02*rns-0.6073 # mostly accurate - each loop is less accurate
if (rns > 0.0875):
  E = -35.5*rns*rns*rns*rns+19.8*rns*rns*rns+4.196*rns*rns+1.528*rns+0.5088
  if (rns > 0.1375):
    E = -602.6*rns*rns*rns*rns+532.4*rns*rns*rns-178*rns*rns+27.62*rns-0.7644
    if (rns > 0.2875):
      E = 0.02667*rns+0.9913
      if (rns > 0.325):
        E = 1

# Solver settings
divisions=5 # For output files in standard format then d must be 5, 9, 13, 17, 21...
maxiterations=1e4
diagfreq=1e3
converged=1e-5*diagfreq # Currently n K difference over i iterations
damping=0.425
edge=0.5
corner=0.25
mutationiters=100



# Initialise
cells_theta1 = np.empty((cellsz, cellsy, cellsx))
cells_theta2 = np.empty((cellsz, cellsy, cellsx))
theta1 = np.empty((cellsz, cellsy, cellsx))
theta2 = np.empty((cellsz, cellsy, cellsx))

pointsx=cellsx*(divisions-1)+1
pointsy=cellsy*(divisions-1)+1
pointsz=cellsz*(divisions-1)+1

points=np.zeros((pointsz,pointsy,pointsx))
for z in range(pointsz):
    points[z,:,:]=400-z*100/(pointsz-1)
pointsprev=points.copy()

xu_lookup=np.empty(pointsx,np.int32)
xl_lookup=np.empty(pointsx,np.int32)
for x in range(pointsx):
    xu_lookup[x]=math.floor(x/(divisions-1))
    xl_lookup[x]=math.ceil(x/(divisions-1))-1

yu_lookup=np.empty(pointsy,np.int32)
yl_lookup=np.empty(pointsy,np.int32)
for y in range(pointsy):
    yu_lookup[y]=math.floor(y/(divisions-1))
    yl_lookup[y]=math.ceil(y/(divisions-1))-1

zu_lookup=np.empty(pointsz,np.int32)
zl_lookup=np.empty(pointsz,np.int32)
for z in range(pointsz):
    zu_lookup[z]=math.floor(z/(divisions-1))
    zl_lookup[z]=math.ceil(z/(divisions-1))-1



@njit(parallel=True, fastmath=True)
def solve(points, newpoints, storage,
    cells_k_xx, cells_k_xy, cells_k_yx, cells_k_yy, cells_k_xz, cells_k_zx, cells_k_zz, cells_k_zy, cells_k_yz,
    xu_lookup, xl_lookup, yu_lookup, yl_lookup, zu_lookup, zl_lookup):
    for i in range(maxiterations):
        regressionlong=0.0
        for z in prange(pointsz):
            for y in range(pointsy):
                for x in range(pointsx):
                    if z==0 or z==pointsz-1:
                        continue
                    # General stencil
                    if x!=0 and x!=pointsx-1 and y!=0 and y!=pointsy-1:
                        xu=xu_lookup[x]; xl=xl_lookup[x]
                        yu=yu_lookup[y]; yl=yl_lookup[y]
                        zu=zu_lookup[z]; zl=zl_lookup[z]
                        pc=points[z,y,x] #1,6,12,8
                        pu=points[z+1,y,x]; pd=points[z-1,y,x]; pl=points[z,y,x-1]; pr=points[z,y,x+1]; pf=points[z,y+1,x]; pb=points[z,y-1,x]
                        pfl=points[z,y+1,x-1]; pfr=points[z,y+1,x+1]; pbl=points[z,y-1,x-1]; pbr=points[z,y-1,x+1]
                        puf=points[z+1,y+1,x]; pur=points[z+1,y,x+1]; pub=points[z+1,y-1,x]; pul=points[z+1,y,x-1]
                        pdf=points[z-1,y+1,x]; pdr=points[z-1,y,x+1]; pdb=points[z-1,y-1,x]; pdl=points[z-1,y,x-1]
                        newpoints[z][y][x] = pc + damping*0.0625*(4*(pu-pc)*(cells_k_zz[zu,yu,xl]+cells_k_zz[zu,yu,xu]+cells_k_zz[zu,yl,xl]+cells_k_zz[zu,yl,xu]) + (pf+puf-pb-pub)*(cells_k_zy[zu,yu,xl]+cells_k_zy[zu,yu,xu]+cells_k_zy[zu,yl,xl]+cells_k_zy[zu,yl,xu]) + (pr+pur-pl-pul)*(cells_k_zx[zu,yu,xl]+cells_k_zx[zu,yu,xu]+cells_k_zx[zu,yl,xl]+cells_k_zx[zu,yl,xu])
                                                           + 4*(pd-pc)*(cells_k_zz[zl,yu,xl]+cells_k_zz[zl,yu,xu]+cells_k_zz[zl,yl,xl]+cells_k_zz[zl,yl,xu]) + (pb+pdb-pf-pdf)*(cells_k_zy[zl,yu,xl]+cells_k_zy[zl,yu,xu]+cells_k_zy[zl,yl,xl]+cells_k_zy[zl,yl,xu]) + (pl+pdl-pr-pdr)*(cells_k_zx[zl,yu,xl]+cells_k_zx[zl,yu,xu]+cells_k_zx[zl,yl,xl]+cells_k_zx[zl,yl,xu])
                                                           + 4*(pf-pc)*(cells_k_yy[zu,yu,xl]+cells_k_yy[zu,yu,xu]+cells_k_yy[zl,yu,xl]+cells_k_yy[zl,yu,xu]) + (pr+pfr-pl-pfl)*(cells_k_yx[zu,yu,xl]+cells_k_yx[zu,yu,xu]+cells_k_yx[zl,yu,xl]+cells_k_yx[zl,yu,xu]) + (pu+puf-pd-pdf)*(cells_k_yz[zu,yu,xl]+cells_k_yz[zu,yu,xu]+cells_k_yz[zl,yu,xl]+cells_k_yz[zl,yu,xu])
                                                           + 4*(pb-pc)*(cells_k_yy[zu,yl,xl]+cells_k_yy[zu,yl,xu]+cells_k_yy[zl,yl,xl]+cells_k_yy[zl,yl,xu]) + (pl+pbl-pr-pbr)*(cells_k_yx[zu,yl,xl]+cells_k_yx[zu,yl,xu]+cells_k_yx[zl,yl,xl]+cells_k_yx[zl,yl,xu]) + (pd+pdb-pu-pub)*(cells_k_yz[zu,yl,xl]+cells_k_yz[zu,yl,xu]+cells_k_yz[zl,yl,xl]+cells_k_yz[zl,yl,xu])
                                                           + 4*(pr-pc)*(cells_k_xx[zu,yu,xu]+cells_k_xx[zu,yl,xu]+cells_k_xx[zl,yu,xu]+cells_k_xx[zl,yl,xu]) + (pu+pur-pd-pdr)*(cells_k_xz[zu,yu,xu]+cells_k_xz[zu,yl,xu]+cells_k_xz[zl,yu,xu]+cells_k_xz[zl,yl,xu]) + (pf+pfr-pb-pbr)*(cells_k_xy[zu,yu,xu]+cells_k_xy[zu,yl,xu]+cells_k_xy[zl,yu,xu]+cells_k_xy[zl,yl,xu])
                                                           + 4*(pl-pc)*(cells_k_xx[zu,yu,xl]+cells_k_xx[zu,yl,xl]+cells_k_xx[zl,yu,xl]+cells_k_xx[zl,yl,xl]) + (pd+pdl-pu-pul)*(cells_k_xz[zu,yu,xl]+cells_k_xz[zu,yl,xl]+cells_k_xz[zl,yu,xl]+cells_k_xz[zl,yl,xl]) + (pb+pbl-pf-pfl)*(cells_k_xy[zu,yu,xl]+cells_k_xy[zu,yl,xl]+cells_k_xy[zl,yu,xl]+cells_k_xy[zl,yl,xl]))
                    # 4 side stencils
                    elif x!=0 and x!=pointsx-1 and y==0 and y!=pointsy-1:
                        xu=xu_lookup[x]; xl=xl_lookup[x]
                        yu=yu_lookup[y]
                        zu=zu_lookup[z]; zl=zl_lookup[z]
                        pc=points[z,y,x]
                        pu=points[z+1,y,x]; pd=points[z-1,y,x]; pl=points[z,y,x-1]; pr=points[z,y,x+1]; pf=points[z,y+1,x]
                        pfl=points[z,y+1,x-1]; pfr=points[z,y+1,x+1]
                        puf=points[z+1,y+1,x]; pur=points[z+1,y,x+1]; pul=points[z+1,y,x-1]
                        pdf=points[z-1,y+1,x]; pdr=points[z-1,y,x+1]; pdl=points[z-1,y,x-1]
                        newpoints[z][y][x] = pc + damping*(0.5*(pu-pc)*(cells_k_zz[zu,yu,xl]+cells_k_zz[zu,yu,xu])*edge + 0.25*(pf+puf-pc-pu)*(cells_k_zy[zu,yu,xl]+cells_k_zy[zu,yu,xu])*edge + 0.5*0.25*(pr+pur-pl-pul)*(cells_k_zx[zu,yu,xl]+cells_k_zx[zu,yu,xu])*edge
                                                           + 0.5*(pd-pc)*(cells_k_zz[zl,yu,xl]+cells_k_zz[zl,yu,xu])*edge + 0.25*(pc+pd-pf-pdf)*(cells_k_zy[zl,yu,xl]+cells_k_zy[zl,yu,xu])*edge + 0.5*0.25*(pl+pdl-pr-pdr)*(cells_k_zx[zl,yu,xl]+cells_k_zx[zl,yu,xu])*edge
                                                           + 0.25*(pf-pc)*(cells_k_yy[zu,yu,xl]+cells_k_yy[zu,yu,xu]+cells_k_yy[zl,yu,xl]+cells_k_yy[zl,yu,xu]) + 0.25*0.25*(pr+pfr-pl-pfl)*(cells_k_yx[zu,yu,xl]+cells_k_yx[zu,yu,xu]+cells_k_yx[zl,yu,xl]+cells_k_yx[zl,yu,xu]) + 0.25*0.25*(pu+puf-pd-pdf)*(cells_k_yz[zu,yu,xl]+cells_k_yz[zu,yu,xu]+cells_k_yz[zl,yu,xl]+cells_k_yz[zl,yu,xu])
                                                           + 0.5*(pr-pc)*(cells_k_xx[zu,yu,xu]+cells_k_xx[zl,yu,xu])*edge + 0.5*0.25*(pu+pur-pd-pdr)*(cells_k_xz[zu,yu,xu]+cells_k_xz[zl,yu,xu])*edge + 0.25*(pf+pfr-pc-pr)*(cells_k_xy[zu,yu,xu]+cells_k_xy[zl,yu,xu])*edge
                                                           + 0.5*(pl-pc)*(cells_k_xx[zu,yu,xl]+cells_k_xx[zl,yu,xl])*edge + 0.5*0.25*(pd+pdl-pu-pul)*(cells_k_xz[zu,yu,xl]+cells_k_xz[zl,yu,xl])*edge + 0.25*(pc+pl-pf-pfl)*(cells_k_xy[zu,yu,xl]+cells_k_xy[zl,yu,xl])*edge)
                    elif x!=0 and x!=pointsx-1 and y!=0 and y==pointsy-1:
                        xu=xu_lookup[x]; xl=xl_lookup[x]
                        yl=yl_lookup[y]
                        zu=zu_lookup[z]; zl=zl_lookup[z]
                        pc=points[z,y,x]
                        pu=points[z+1,y,x]; pd=points[z-1,y,x]; pl=points[z,y,x-1]; pr=points[z,y,x+1]; pb=points[z,y-1,x]
                        pbl=points[z,y-1,x-1]; pbr=points[z,y-1,x+1]
                        pur=points[z+1,y,x+1]; pub=points[z+1,y-1,x]; pul=points[z+1,y,x-1]
                        pdr=points[z-1,y,x+1]; pdb=points[z-1,y-1,x]; pdl=points[z-1,y,x-1]
                        newpoints[z][y][x] = pc + damping*(0.5*(pu-pc)*(cells_k_zz[zu,yl,xl]+cells_k_zz[zu,yl,xu])*edge + 0.25*(pc+pu-pb-pub)*(cells_k_zy[zu,yl,xl]+cells_k_zy[zu,yl,xu])*edge + 0.5*0.25*(pr+pur-pl-pul)*(cells_k_zx[zu,yl,xl]+cells_k_zx[zu,yl,xu])*edge
                                                           + 0.5*(pd-pc)*(cells_k_zz[zl,yl,xl]+cells_k_zz[zl,yl,xu])*edge + 0.25*(pb+pdb-pc-pd)*(cells_k_zy[zl,yl,xl]+cells_k_zy[zl,yl,xu])*edge + 0.5*0.25*(pl+pdl-pr-pdr)*(cells_k_zx[zl,yl,xl]+cells_k_zx[zl,yl,xu])*edge
                                                           + 0.25*(pb-pc)*(cells_k_yy[zu,yl,xl]+cells_k_yy[zu,yl,xu]+cells_k_yy[zl,yl,xl]+cells_k_yy[zl,yl,xu]) + 0.25*0.25*(pl+pbl-pr-pbr)*(cells_k_yx[zu,yl,xl]+cells_k_yx[zu,yl,xu]+cells_k_yx[zl,yl,xl]+cells_k_yx[zl,yl,xu]) + 0.25*0.25*(pd+pdb-pu-pub)*(cells_k_yz[zu,yl,xl]+cells_k_yz[zu,yl,xu]+cells_k_yz[zl,yl,xl]+cells_k_yz[zl,yl,xu])
                                                           + 0.5*(pr-pc)*(cells_k_xx[zu,yl,xu]+cells_k_xx[zl,yl,xu])*edge + 0.5*0.25*(pu+pur-pd-pdr)*(cells_k_xz[zu,yl,xu]+cells_k_xz[zl,yl,xu])*edge + 0.25*(pc+pr-pb-pbr)*(cells_k_xy[zu,yl,xu]+cells_k_xy[zl,yl,xu])*edge
                                                           + 0.5*(pl-pc)*(cells_k_xx[zu,yl,xl]+cells_k_xx[zl,yl,xl])*edge + 0.5*0.25*(pd+pdl-pu-pul)*(cells_k_xz[zu,yl,xl]+cells_k_xz[zl,yl,xl])*edge + 0.25*(pb+pbl-pc-pl)*(cells_k_xy[zu,yl,xl]+cells_k_xy[zl,yl,xl])*edge)
                    elif x==0 and x!=pointsx-1 and y!=0 and y!=pointsy-1:
                        xu=xu_lookup[x]
                        yu=yu_lookup[y]; yl=yl_lookup[y]
                        zu=zu_lookup[z]; zl=zl_lookup[z]
                        pc=points[z,y,x]
                        pu=points[z+1,y,x]; pd=points[z-1,y,x]; pr=points[z,y,x+1]; pf=points[z,y+1,x]; pb=points[z,y-1,x]
                        pfr=points[z,y+1,x+1]; pbr=points[z,y-1,x+1]
                        puf=points[z+1,y+1,x]; pur=points[z+1,y,x+1]; pub=points[z+1,y-1,x]
                        pdf=points[z-1,y+1,x]; pdr=points[z-1,y,x+1]; pdb=points[z-1,y-1,x]
                        newpoints[z][y][x] = pc + damping*(0.5*(pu-pc)*(cells_k_zz[zu,yu,xu]+cells_k_zz[zu,yl,xu])*edge + 0.5*0.25*(pf+puf-pb-pub)*(cells_k_zy[zu,yu,xu]+cells_k_zy[zu,yl,xu])*edge + 0.25*(pr+pur-pc-pu)*(cells_k_zx[zu,yu,xu]+cells_k_zx[zu,yl,xu])*edge
                                                          + 0.5*(pd-pc)*(cells_k_zz[zl,yu,xu]+cells_k_zz[zl,yl,xu])*edge + 0.5*0.25*(pb+pdb-pf-pdf)*(cells_k_zy[zl,yu,xu]+cells_k_zy[zl,yl,xu])*edge + 0.25*(pc+pd-pr-pdr)*(cells_k_zx[zl,yu,xu]+cells_k_zx[zl,yl,xu])*edge
                                                          + 0.5*(pf-pc)*(cells_k_yy[zu,yu,xu]+cells_k_yy[zl,yu,xu])*edge + 0.25*(pr+pfr-pc-pf)*(cells_k_yx[zu,yu,xu]+cells_k_yx[zl,yu,xu])*edge + 0.5*0.25*(pu+puf-pd-pdf)*(cells_k_yz[zu,yu,xu]+cells_k_yz[zl,yu,xu])*edge
                                                          + 0.5*(pb-pc)*(cells_k_yy[zu,yl,xu]+cells_k_yy[zl,yl,xu])*edge + 0.25*(pc+pb-pr-pbr)*(cells_k_yx[zu,yl,xu]+cells_k_yx[zl,yl,xu])*edge + 0.5*0.25*(pd+pdb-pu-pub)*(cells_k_yz[zu,yl,xu]+cells_k_yz[zl,yl,xu])*edge
                                                          + 0.25*(pr-pc)*(cells_k_xx[zu,yu,xu]+cells_k_xx[zu,yl,xu]+cells_k_xx[zl,yu,xu]+cells_k_xx[zl,yl,xu]) + 0.25*0.25*(pu+pur-pd-pdr)*(cells_k_xz[zu,yu,xu]+cells_k_xz[zu,yl,xu]+cells_k_xz[zl,yu,xu]+cells_k_xz[zl,yl,xu]) + 0.25*0.25*(pf+pfr-pb-pbr)*(cells_k_xy[zu,yu,xu]+cells_k_xy[zu,yl,xu]+cells_k_xy[zl,yu,xu]+cells_k_xy[zl,yl,xu]))
                    elif x!=0 and x==pointsx-1 and y!=0 and y!=pointsy-1:
                        xl=xl_lookup[x]
                        yu=yu_lookup[y]; yl=yl_lookup[y]
                        zu=zu_lookup[z]; zl=zl_lookup[z]
                        pc=points[z,y,x]
                        pu=points[z+1,y,x]; pd=points[z-1,y,x]; pl=points[z,y,x-1]; pf=points[z,y+1,x]; pb=points[z,y-1,x]
                        pfl=points[z,y+1,x-1]; pbl=points[z,y-1,x-1]
                        puf=points[z+1,y+1,x]; pub=points[z+1,y-1,x]; pul=points[z+1,y,x-1]
                        pdf=points[z-1,y+1,x]; pdb=points[z-1,y-1,x]; pdl=points[z-1,y,x-1]
                        newpoints[z][y][x] = pc + damping*(0.5*(pu-pc)*(cells_k_zz[zu,yu,xl]+cells_k_zz[zu,yl,xl])*edge + 0.5*0.25*(pf+puf-pb-pub)*(cells_k_zy[zu,yu,xl]+cells_k_zy[zu,yl,xl])*edge + 0.25*(pc+pu-pl-pul)*(cells_k_zx[zu,yu,xl]+cells_k_zx[zu,yl,xl])*edge
                                                           + 0.5*(pd-pc)*(cells_k_zz[zl,yu,xl]+cells_k_zz[zl,yl,xl])*edge + 0.5*0.25*(pb+pdb-pf-pdf)*(cells_k_zy[zl,yu,xl]+cells_k_zy[zl,yl,xl])*edge + 0.25*(pl+pdl-pc-pd)*(cells_k_zx[zl,yu,xl]+cells_k_zx[zl,yl,xl])*edge
                                                           + 0.5*(pf-pc)*(cells_k_yy[zu,yu,xl]+cells_k_yy[zl,yu,xl])*edge + 0.25*(pc+pf-pl-pfl)*(cells_k_yx[zu,yu,xl]+cells_k_yx[zl,yu,xl])*edge + 0.5*0.25*(pu+puf-pd-pdf)*(cells_k_yz[zu,yu,xl]+cells_k_yz[zl,yu,xl])*edge
                                                           + 0.5*(pb-pc)*(cells_k_yy[zu,yl,xl]+cells_k_yy[zl,yl,xl])*edge + 0.25*(pl+pbl-pc-pb)*(cells_k_yx[zu,yl,xl]+cells_k_yx[zl,yl,xl])*edge + 0.5*0.25*(pd+pdb-pu-pub)*(cells_k_yz[zu,yl,xl]+cells_k_yz[zl,yl,xl])*edge
                                                           + 0.25*(pl-pc)*(cells_k_xx[zu,yu,xl]+cells_k_xx[zu,yl,xl]+cells_k_xx[zl,yu,xl]+cells_k_xx[zl,yl,xl]) + 0.25*0.25*(pd+pdl-pu-pul)*(cells_k_xz[zu,yu,xl]+cells_k_xz[zu,yl,xl]+cells_k_xz[zl,yu,xl]+cells_k_xz[zl,yl,xl]) + 0.25*0.25*(pb+pbl-pf-pfl)*(cells_k_xy[zu,yu,xl]+cells_k_xy[zu,yl,xl]+cells_k_xy[zl,yu,xl]+cells_k_xy[zl,yl,xl]))
                    # 4 corner stencils
                    elif x==0 and x!=pointsx-1 and y==0 and y!=pointsy-1:
                        xu=xu_lookup[x]
                        yu=yu_lookup[y]
                        zu=zu_lookup[z]; zl=zl_lookup[z]
                        pc=points[z,y,x]
                        pu=points[z+1,y,x]; pd=points[z-1,y,x]; pr=points[z,y,x+1]; pf=points[z,y+1,x]
                        pfr=points[z,y+1,x+1]
                        puf=points[z+1,y+1,x]; pur=points[z+1,y,x+1]
                        pdf=points[z-1,y+1,x]; pdr=points[z-1,y,x+1]
                        newpoints[z][y][x] = pc + damping*((pu-pc)*(cells_k_zz[zu,yu,xu])*corner + 0.5*(pf+puf-pc-pu)*(cells_k_zy[zu,yu,xu])*corner + 0.5*(pr+pur-pc-pu)*(cells_k_zx[zu,yu,xu])*corner
                                                           + (pd-pc)*(cells_k_zz[zl,yu,xu])*corner + 0.5*(pc+pd-pf-pdf)*(cells_k_zy[zl,yu,xu])*corner + 0.5*(pc+pd-pr-pdr)*(cells_k_zx[zl,yu,xu])*corner
                                                           + 0.5*(pf-pc)*(cells_k_yy[zu,yu,xu]+cells_k_yy[zl,yu,xu])*edge + 0.25*(pr+pfr-pc-pf)*(cells_k_yx[zu,yu,xu]+cells_k_yx[zl,yu,xu])*edge + 0.5*0.25*(pu+puf-pd-pdf)*(cells_k_yz[zu,yu,xu]+cells_k_yz[zl,yu,xu])*edge
                                                           + 0.5*(pr-pc)*(cells_k_xx[zu,yu,xu]+cells_k_xx[zl,yu,xu])*edge + 0.5*0.25*(pu+pur-pd-pdr)*(cells_k_xz[zu,yu,xu]+cells_k_xz[zl,yu,xu])*edge + 0.25*(pf+pfr-pc-pr)*(cells_k_xy[zu,yu,xu]+cells_k_xy[zl,yu,xu])*edge)
                    elif x!=0 and x==pointsx-1 and y==0 and y!=pointsy-1:
                        xl=xl_lookup[x]
                        yu=yu_lookup[y]
                        zu=zu_lookup[z]; zl=zl_lookup[z]
                        pc=points[z,y,x]
                        pu=points[z+1,y,x]; pd=points[z-1,y,x]; pl=points[z,y,x-1]; pf=points[z,y+1,x]
                        pfl=points[z,y+1,x-1]
                        puf=points[z+1,y+1,x]; pul=points[z+1,y,x-1]
                        pdf=points[z-1,y+1,x]; pdl=points[z-1,y,x-1]
                        newpoints[z][y][x] = pc + damping*((pu-pc)*(cells_k_zz[zu,yu,xl])*corner + 0.5*(pf+puf-pc-pu)*(cells_k_zy[zu,yu,xl])*corner + 0.5*(pc+pu-pl-pul)*(cells_k_zx[zu,yu,xl])*corner
                                                           + (pd-pc)*(cells_k_zz[zl,yu,xl])*corner + 0.5*(pc+pd-pf-pdf)*(cells_k_zy[zl,yu,xl])*corner + 0.5*(pl+pdl-pc-pd)*(cells_k_zx[zl,yu,xl])*corner
                                                           + 0.5*(pf-pc)*(cells_k_yy[zu,yu,xl]+cells_k_yy[zl,yu,xl])*edge + 0.25*(pc+pf-pl-pfl)*(cells_k_yx[zu,yu,xl]+cells_k_yx[zl,yu,xl])*edge + 0.5*0.25*(pu+puf-pd-pdf)*(cells_k_yz[zu,yu,xl]+cells_k_yz[zl,yu,xl])*edge
                                                           + 0.5*(pl-pc)*(cells_k_xx[zu,yu,xl]+cells_k_xx[zl,yu,xl])*edge + 0.5*0.25*(pd+pdl-pu-pul)*(cells_k_xz[zu,yu,xl]+cells_k_xz[zl,yu,xl])*edge + 0.25*(pc+pl-pf-pfl)*(cells_k_xy[zu,yu,xl]+cells_k_xy[zl,yu,xl])*edge)
                    elif x==0 and x!=pointsx-1 and y!=0 and y==pointsy-1:
                        xu=xu_lookup[x]
                        yl=yl_lookup[y]
                        zu=zu_lookup[z]; zl=zl_lookup[z]
                        pc=points[z,y,x]
                        pu=points[z+1,y,x]; pd=points[z-1,y,x]; pr=points[z,y,x+1]; pb=points[z,y-1,x]
                        pbr=points[z,y-1,x+1]
                        pur=points[z+1,y,x+1]; pub=points[z+1,y-1,x]
                        pdr=points[z-1,y,x+1]; pdb=points[z-1,y-1,x]
                        newpoints[z][y][x] = pc + damping*((pu-pc)*(cells_k_zz[zu,yl,xu])*corner + 0.5*(pc+pu-pb-pub)*(cells_k_zy[zu,yl,xu])*corner + 0.5*(pr+pur-pc-pu)*(cells_k_zx[zu,yl,xu])*corner
                                                           + (pd-pc)*(cells_k_zz[zl,yl,xu])*corner + 0.5*(pb+pdb-pc-pd)*(cells_k_zy[zl,yl,xu])*corner + 0.5*(pc+pd-pr-pdr)*(cells_k_zx[zl,yl,xu])*corner
                                                           + 0.5*(pb-pc)*(cells_k_yy[zu,yl,xu]+cells_k_yy[zl,yl,xu])*edge + 0.25*(pc+pb-pr-pbr)*(cells_k_yx[zu,yl,xu]+cells_k_yx[zl,yl,xu])*edge + 0.5*0.25*(pd+pdb-pu-pub)*(cells_k_yz[zu,yl,xu]+cells_k_yz[zl,yl,xu])*edge
                                                           + 0.5*(pr-pc)*(cells_k_xx[zu,yl,xu]+cells_k_xx[zl,yl,xu])*edge + 0.5*0.25*(pu+pur-pd-pdr)*(cells_k_xz[zu,yl,xu]+cells_k_xz[zl,yl,xu])*edge + 0.25*(pc+pr-pb-pbr)*(cells_k_xy[zu,yl,xu]+cells_k_xy[zl,yl,xu])*edge)
                    elif x!=0 and x==pointsx-1 and y!=0 and y==pointsy-1: # can make this else to optimise a tad
                        xl=xl_lookup[x]
                        yl=yl_lookup[y]
                        zu=zu_lookup[z]; zl=zl_lookup[z]
                        pc=points[z,y,x]
                        pu=points[z+1,y,x]; pd=points[z-1,y,x]; pl=points[z,y,x-1]; pb=points[z,y-1,x]
                        pbl=points[z,y-1,x-1]; 
                        pub=points[z+1,y-1,x]; pul=points[z+1,y,x-1]
                        pdb=points[z-1,y-1,x]; pdl=points[z-1,y,x-1]
                        newpoints[z][y][x] = pc + damping*((pu-pc)*(cells_k_zz[zu,yl,xl])*corner + 0.5*(pc+pu-pb-pub)*(cells_k_zy[zu,yl,xl])*corner + 0.5*(pc+pu-pl-pul)*(cells_k_zx[zu,yl,xl])*corner
                                                           + (pd-pc)*(cells_k_zz[zl,yl,xl])*corner + 0.5*(pb+pdb-pc-pd)*(cells_k_zy[zl,yl,xl])*corner + 0.5*(pl+pdl-pc-pd)*(cells_k_zx[zl,yl,xl])*corner
                                                           + 0.5*(pb-pc)*(cells_k_yy[zu,yl,xl]+cells_k_yy[zl,yl,xl])*edge + 0.25*(pl+pbl-pc-pb)*(cells_k_yx[zu,yl,xl]+cells_k_yx[zl,yl,xl])*edge + 0.5*0.25*(pd+pdb-pu-pub)*(cells_k_yz[zu,yl,xl]+cells_k_yz[zl,yl,xl])*edge
                                                           + 0.5*(pl-pc)*(cells_k_xx[zu,yl,xl]+cells_k_xx[zl,yl,xl])*edge + 0.5*0.25*(pd+pdl-pu-pul)*(cells_k_xz[zu,yl,xl]+cells_k_xz[zl,yl,xl])*edge + 0.25*(pb+pbl-pc-pl)*(cells_k_xy[zu,yl,xl]+cells_k_xy[zl,yl,xl])*edge)
                    # Diagnostics
                    if i==0 or (i+1)%diagfreq==0:
                        regressionlong+=abs(newpoints[z,y,x]-storage[z,y,x])
                        storage[z,y,x]=newpoints[z,y,x]        
        if (i+1)%diagfreq==0:
            regressionlongavg = regressionlong/(pointsz*pointsx*pointsy)
            #print(i+1, regressionlongavg)
            if regressionlongavg<converged:
                print("\nConverged at iteration", i+1, regressionlongavg)
                break
        points = newpoints
    return points



def run_simulation(theta1, theta2):
    theta1 = np.deg2rad(theta1)
    theta2 = np.deg2rad(theta2)

    c1 = np.cos(theta1)
    c2 = np.cos(theta2)
    s1 = np.sin(theta1)
    s2 = np.sin(theta2)
    cells_k_xx = s1*s1*c2*c2+E*c1*c1*c2*c2+E*s2*s2
    cells_k_yy = s1*s1*s2*s2+E*c1*c1*s2*s2+E*c2*c2
    cells_k_zz = c1*c1+E*s1*s1
    cells_k_yx = (s1*s1+E*c1*c1-E)*s2*c2
    cells_k_xy = cells_k_yx.copy()
    cells_k_zx = (1-E)*s1*c1*c2
    cells_k_xz = cells_k_zx.copy()
    cells_k_zy = (1-E)*s1*c1*s2
    cells_k_yz = cells_k_zy.copy()

    global pointsprev
    points = pointsprev.copy()
    newpoints=points.copy()
    storage=points.copy()

    points = solve(points, newpoints, storage,
                   cells_k_xx, cells_k_xy, cells_k_yx, cells_k_yy, cells_k_xz, cells_k_zx, cells_k_zz, cells_k_zy, cells_k_yz,
                   xu_lookup, xl_lookup, yu_lookup, yl_lookup, zu_lookup, zl_lookup)
    
    pointsprev[:] = points
    return points



def objective(line):
    n = cellsx*cellsy*cellsz
    theta1 = line[:n].reshape((cellsz,cellsy,cellsx))
    theta2 = line[n:].reshape((cellsz,cellsy,cellsx))

    global points
    points = run_simulation(theta1, theta2)

    totaltemp1=0 # select all points in cloaked sphere and average their temps
    npoints1=0
    var1=0 # iterate through all points in cloaked sphere and just do abs(temp-average) for each
    if cloak1radius != 0:
        for z in range(pointsz):
            for y in range(pointsy):
                for x in range(pointsx):
                    if (cloak1location[0]-x*side/(divisions-1))**2+(cloak1location[1]-y*side/(divisions-1))**2+(cloak1location[2]-z*side/(divisions-1))**2 <= cloak1radius**2:
                        totaltemp1 += points[z,y,x]
                        npoints1 += 1
        for z in range(pointsz):
            for y in range(pointsy):
                for x in range(pointsx):
                    if (cloak1location[0]-x*side/(divisions-1))**2+(cloak1location[1]-y*side/(divisions-1))**2+(cloak1location[2]-z*side/(divisions-1))**2 <= cloak1radius**2:
                        if npoints1 != 0:
                            var1 += abs(points[z,y,x]-totaltemp1/npoints1)/npoints1
        if npoints1 == 0:
            print("cloak 1 too small to reach a grid point --> Increase number of divisions or increase cloak size")

    totaltemp2=0
    npoints2=0
    var2=0
    if cloak2radius != 0:
        for z in range(pointsz):
            for y in range(pointsy):
                for x in range(pointsx):
                    if (cloak2location[0]-x*side/(divisions-1))**2+(cloak2location[1]-y*side/(divisions-1))**2+(cloak2location[2]-z*side/(divisions-1))**2 <= cloak2radius**2:
                        totaltemp2 += points[z,y,x]
                        npoints2 += 1
        for z in range(pointsz):
            for y in range(pointsy):
                for x in range(pointsx):
                    if (cloak2location[0]-x*side/(divisions-1))**2+(cloak2location[1]-y*side/(divisions-1))**2+(cloak2location[2]-z*side/(divisions-1))**2 <= cloak2radius**2:
                        if npoints2 != 0:
                            var2 += abs(points[z,y,x]-totaltemp2/npoints2)/npoints2
        if npoints2 == 0:
            print("cloak 2 too small to reach a grid point --> Increase number of divisions or increase cloak size")

    totaltemp3=0
    npoints3=0
    var3=0
    if cloak3radius != 0:
        for z in range(pointsz):
            for y in range(pointsy):
                for x in range(pointsx):
                    if (cloak3location[0]-x*side/(divisions-1))**2+(cloak3location[1]-y*side/(divisions-1))**2+(cloak3location[2]-z*side/(divisions-1))**2 <= cloak3radius**2:
                        totaltemp3 += points[z,y,x]
                        npoints3 += 1
        for z in range(pointsz):
            for y in range(pointsy):
                for x in range(pointsx):
                    if (cloak3location[0]-x*side/(divisions-1))**2+(cloak3location[1]-y*side/(divisions-1))**2+(cloak3location[2]-z*side/(divisions-1))**2 <= cloak3radius**2:
                        if npoints3 != 0:
                            var3 += abs(points[z,y,x]-totaltemp3/npoints3)/npoints3
        if npoints3 == 0:
            print("cloak 3 too small to reach a grid point --> Increase number of divisions or increase cloak size")

    var4=0
    ft1count=0
    if ft1on and ft1point:
        if ft1positive:
            var4 = abs(350-points[round(ft1location[2]*(divisions-1)/side),round(ft1location[1]*(divisions-1)/side),round(ft1location[0]*(divisions-1)/side)]) #points[y,x]
        else:
            var4 = abs(points[round(ft1location[2]*(divisions-1)/side),round(ft1location[1]*(divisions-1)/side),round(ft1location[0]*(divisions-1)/side)]-300)
    elif ft1on and not ft1point:
        if ft1positive:
            for z in range(pointsz):
                for y in range(pointsy):
                    for x in range(pointsx):
                        if (ft1location[0]-x*side/(divisions-1))**2+(ft1location[1]-y*side/(divisions-1))**2+(ft1location[2]-z*side/(divisions-1))**2 <= ft1radius**2:
                            ft1count += 1
                            var4 += abs(350-points[z,y,x])#####
        else:
            for z in range(pointsz):
                for y in range(pointsy):
                    for x in range(pointsx):
                        if (ft1location[0]-x*side/(divisions-1))**2+(ft1location[1]-y*side/(divisions-1))**2+(ft1location[2]-z*side/(divisions-1))**2 <= ft1radius**2:
                            ft1count += 1
                            var4 += abs(points[z,y,x]-350)
        if ft1count != 0:
            var4 = var4/ft1count
        else:
            print("Feature 1 too small to reach a grid point --> Increase number of divisions or increase feature radius")

    var5=0
    ft2count=0
    if ft2on and ft2point:
        if ft2positive:
            var5 = abs(400-points[round(ft2location[2]*(divisions-1)/side),round(ft2location[1]*(divisions-1)/side),round(ft2location[0]*(divisions-1)/side)])
        else:
            var5 = abs(points[round(ft2location[2]*(divisions-1)/side),round(ft2location[1]*(divisions-1)/side),round(ft2location[0]*(divisions-1)/side)]-300)
    elif ft2on and not ft2point:
        if ft2positive:
            for z in range(pointsz):
                for y in range(pointsy):
                    for x in range(pointsx):
                        if (ft2location[0]-x*side/(divisions-1))**2+(ft2location[1]-y*side/(divisions-1))**2+(ft2location[2]-z*side/(divisions-1))**2 <= ft2radius**2:
                            ft2count += 1
                            var5 += abs(400-points[z,y,x])
        else:
            for z in range(pointsz):
                for y in range(pointsy):
                    for x in range(pointsx):
                        if (ft2location[0]-x*side/(divisions-1))**2+(ft2location[1]-y*side/(divisions-1))**2+(ft2location[2]-z*side/(divisions-1))**2 <= ft2radius**2:
                            ft2count += 1
                            var5 += abs(points[z,y,x]-300)
        if ft2count != 0:
            var5 = var5/ft2count
        else:
            print("Feature 2 too small to reach a grid point --> Increase number of divisions or increase feature radius")

    var6=0
    ft3count=0
    if ft3on and ft3point:
        if ft3positive:
            var6 = abs(400-points[round(ft3location[2]*(divisions-1)/side),round(ft3location[1]*(divisions-1)/side),round(ft3location[0]*(divisions-1)/side)])
        else:
            var6 = abs(points[round(ft3location[2]*(divisions-1)/side),round(ft3location[1]*(divisions-1)/side),round(ft3location[0]*(divisions-1)/side)]-300)
    elif ft3on and not ft3point:
        if ft3positive:
            for z in range(pointsz):
                for y in range(pointsy):
                    for x in range(pointsx):
                        if (ft3location[0]-x*side/(divisions-1))**2+(ft3location[1]-y*side/(divisions-1))**2+(ft3location[2]-z*side/(divisions-1))**2 <= ft3radius**2:
                            ft3count += 1
                            var6 += abs(400-points[z,y,x])
        else:
            for z in range(pointsz):
                for y in range(pointsy):
                    for x in range(pointsx):
                        if (ft3location[0]-x*side/(divisions-1))**2+(ft3location[1]-y*side/(divisions-1))**2+(ft3location[2]-z*side/(divisions-1))**2 <= ft3radius**2:
                            ft3count += 1
                            var6 += abs(points[z,y,x]-300)
        if ft3count != 0:
            var6 = var6/ft3count
        else:
            print("Feature 3 too small to reach a grid point --> Increase number of divisions or increase feature radius")

    var7=0
    ftgcount=0
    if ftgon:
      for z in range(pointsy):
        for y in range(pointsy):
          for z in range(pointsz):
            if (ftgexcludep[0]-x*side/(divisions-1))**2+(ftgexcludep[1]-y*side/(divisions-1))**2+(ftgexcludep[2]-z*side/(divisions-1))**2 > ftgexcluder**2:
              ftgcount+=1
              var7 += abs(points[z,y,x]-(400-z*100/(pointsz-1)))#####
      if ftgcount != 0:
        var7 = var7/ftgcount

    vartotal=clw*(var1+var2+var3)+ftw*(2*var4+var5+var6+var7)######

    return vartotal



def optimise():
    if os.path.exists("cma_state.pkl"):
        print("Continuing previous optimisation...")
        with open("cma_state.pkl", "rb") as f:
            es = pickle.load(f)

    else:
        print("Starting new optimisation...")
        for z in range(cellsz):
            theta1[z] = np.loadtxt(f"input1optimised/input1thetas{z+1}.txt")
            theta2[z] = np.loadtxt(f"input2optimised/input2thetas{z+1}.txt")
        theta10 = theta1.flatten()
        theta20 = theta2.flatten()
        theta1020 = np.concatenate((theta10,theta20))
        es = cma.CMAEvolutionStrategy(
            theta1020,
            5.0,
            {
                "popsize": 6
            }
        )

    for generation in range(mutationiters):
        solutions = es.ask()
        scores = [objective(line) for line in solutions]
        es.tell(solutions, scores)

        print(
            f"Generation {es.countiter} | "
            f"Best {es.result.fbest:.6f} | "
            f"sigma {es.sigma:.4f}"
        )

    best = es.result.xbest
    n = cellsx*cellsy*cellsz
    best_theta1 = best[:n].reshape((cellsz,cellsy,cellsx))
    best_theta2 = best[n:].reshape((cellsz,cellsy,cellsx))

    for z in range(cellsz):
        np.savetxt(f"input1optimised/input1thetas{z+1}.txt", best_theta1[z], fmt="%.6f")
        np.savetxt(f"input2optimised/input2thetas{z+1}.txt", best_theta2[z], fmt="%.6f")

    with open("cma_state.pkl", "wb") as f:
        pickle.dump(es, f)
        


if __name__ == "__main__":
    optimise()