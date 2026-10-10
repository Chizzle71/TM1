import math
import numpy as np
import matplotlib.pyplot as plt
from numba import njit, prange, set_num_threads
from matplotlib import colormaps, cm
from matplotlib.patches import Circle, Rectangle, RegularPolygon
from matplotlib.colors import Normalize
from matplotlib.ticker import FuncFormatter
from mpl_toolkits.mplot3d import Axes3D, proj3d
list(colormaps)

plt.rcParams.update({
    'font.family': 'Times New Roman',
    'font.size': 12*1.8
})

# Cores
set_num_threads(18)

# COMSOL related
cellsx=20
cellsy=20
cellsz=20
side=0.08
rns=0.015
E = -43520*rns*rns*rns*rns+11610*rns*rns*rns-1166*rns*rns+57.02*rns-0.6073 # mostly accurate - each indent is less accurate
if (rns > 0.0875):
  E = -35.5*rns*rns*rns*rns+19.8*rns*rns*rns+4.196*rns*rns+1.528*rns+0.5088
  if (rns > 0.1375):
    E = -602.6*rns*rns*rns*rns+532.4*rns*rns*rns-178*rns*rns+27.62*rns-0.7644
    if (rns > 0.2875):
      E = 0.02667*rns+0.9913
      if (rns > 0.325):
        E = 1

# Solver settings
divisions=9 # For output files in standard format then d must be 5, 9, 13, 17, 21...
maxiterations=1e7
diagfreq=1e4
converged=1e-6*diagfreq # Currently n K difference over i iterations
damping=0.425
edge=0.5
corner=0.25

# Simulation inputs
simulation="1"
folder=f"{simulation}"

cells_theta1 = np.empty((cellsz, cellsy, cellsx))
cells_theta2 = np.empty((cellsz, cellsy, cellsx))
theta1 = np.empty((cellsz, cellsy, cellsx))
theta2 = np.empty((cellsz, cellsy, cellsx))

for z in range(cellsz):
    theta1[z] = np.deg2rad(np.loadtxt(f"input1optimised/input1thetas{z+1}.txt"))
    theta2[z] = np.deg2rad(np.loadtxt(f"input2optimised/input2thetas{z+1}.txt"))

# Initialise
pointsx=cellsx*(divisions-1)+1
pointsy=cellsy*(divisions-1)+1
pointsz=cellsz*(divisions-1)+1

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

points=np.zeros((pointsz,pointsy,pointsx))
for z in range(pointsz):
    points[z,:,:]=400-z*100/(pointsz-1)

newpoints=points.copy()
storage=points.copy()

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

# Solve
@njit(parallel=True, fastmath=True) # Floating point reassociation and commutativity on
def solve(points,newpoints,storage):
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
            print(i+1, regressionlongavg)
            if regressionlongavg<converged:
                print("Converged at iteration", i+1)
                break
        points,newpoints=newpoints,points
    return newpoints

points=solve(points,newpoints,storage)









# Font
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman', 'Times', 'Liberation Serif', 'DejaVu Serif']
plt.rcParams['mathtext.fontset'] = 'stix'

#Physical extents of the domain
Lx = cellsx * side
Ly = cellsy * side
Lz = cellsz * side
x_coords = np.linspace(0, Lx, pointsx)
y_coords = np.linspace(0, Ly, pointsy)
z_coords = np.linspace(0, Lz, pointsz)
cmap = cm.seismic#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
vmin, vmax = 300, 400
norm = Normalize(vmin=vmin, vmax=vmax)
ISO_ELEV = 35.264 # true isometric elevation angle
ISO_AZIM = 45
Y_DRAW = Ly / 2 # every slice is drawn here
SLICE_LINE=4

def new_iso_axes(figsize=(7, 7)):
    """Create a 3D axes set up for an isometric view with correct aspect ratio."""
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection='3d')
    ax.view_init(elev=ISO_ELEV, azim=ISO_AZIM)
    ax.set_box_aspect((Lx, Ly, Lz))
    ax.set_xlabel('$x$ (m)', labelpad=12)
    ax.set_ylabel('$y$ (m)', labelpad=12)#[][]
    ax.set_zlabel('$z$ (m)', labelpad=12)
    ax.set_xlim(Lx, 0)
    ax.set_ylim(Ly, 0)
    ax.set_zlim(0, Lz)

    # Fill the figure almost entirely with the axes
    ax.set_position([0.02, 0.02, 0.96, 0.96])

    def format_axis_tick(value, position):
        return '0' if np.isclose(value, 0) else f'{value:.1f}'

    # 5 ticks spaced by 0.4 m
    x_ticks = np.arange(0, Lx + 0.001, 0.4)
    y_ticks = np.arange(0, Ly + 0.001, 0.4)
    z_ticks = np.arange(0, Lz + 0.001, 0.4)
    ax.set_xticks(x_ticks)
    ax.set_yticks(y_ticks)
    ax.set_zticks(z_ticks)
    ax.xaxis.set_major_formatter(FuncFormatter(format_axis_tick))
    ax.yaxis.set_major_formatter(FuncFormatter(format_axis_tick))
    ax.zaxis.set_major_formatter(FuncFormatter(format_axis_tick))

    return fig, ax

def plot_face(ax, X, Y, Z, data, alpha=1.0):
    """Draw one flat, data-colored patch at the given 3D coordinates."""
    facecolors = cmap(norm(data))
    ax.plot_surface(X, Y, Z, facecolors=facecolors, rstride=1, cstride=1,
                    shade=False, linewidth=0, antialiased=False, alpha=alpha)

#Coordinate grids reused for every face
Xxy, Yxy = np.meshgrid(x_coords, y_coords)
Xxz, Zxz = np.meshgrid(x_coords, z_coords)
Yyz, Zyz = np.meshgrid(y_coords, z_coords)

#Image 1 - the cube itself, all 6 outer faces
fig, ax = new_iso_axes()
plot_face(ax, Xxy, Yxy, np.full_like(Xxy, 0.0), points[0, :, :]) # bottom (z=0)
plot_face(ax, Xxy, Yxy, np.full_like(Xxy, Lz), points[-1, :, :]) # top (z=Lz)
plot_face(ax, Xxz, np.full_like(Xxz, 0.0), Zxz, points[:, 0, :]) # front (y=0)
plot_face(ax, Xxz, np.full_like(Xxz, Ly), Zxz, points[:, -1, :]) # back (y=Ly)
plot_face(ax, np.full_like(Yyz, 0.0), Yyz, Zyz, points[:, :, 0]) # left (x=0)
plot_face(ax, np.full_like(Yyz, Lx), Yyz, Zyz, points[:, :, -1]) # right (x=Lx)

# Slice location lines
slice_y_values = [0.8, 0.96, 1.12]
slice_r_values = [0.24, 0.1788854382, 0]#only valid for curent y values
slice_labels = ['C', 'B', 'A']
slice_colours = ['lightgreen', 'yellow', '#F4B183']

line_offset = 0.01 * Ly#############

for y_slice, label, colour in zip(slice_y_values, slice_labels, slice_colours):
    y_line = y_slice - line_offset if np.isclose(y_slice, Ly) else y_slice

    # Dotted line along the left edge
    ax.plot([0, 0], [y_line, y_line], [0, Lz], color=colour, linestyle=':', linewidth=SLICE_LINE, zorder=1000)
    # Dotted line along the top edge, reaching x = 1.6
    ax.plot([0, Lx], [y_line, y_line], [Lz, Lz], color=colour, linestyle=':', linewidth=SLICE_LINE, zorder=1000)
    # Dotted line along the right edge
    ax.plot([Lx, Lx], [y_line, y_line], [0, Lz], color=colour, linestyle=':', linewidth=SLICE_LINE, zorder=1000)
    # Dotted line along the bottom edge
    ax.plot([0, Lx], [y_line, y_line], [0, 0], color=colour, linestyle=':', linewidth=SLICE_LINE, zorder=1000)
    # Label at the far end of the line (x = 1.6)
    ax.text(Lx, y_line, Lz, label, fontsize=19, fontfamily='serif', color='black', ha='left', va='bottom', zorder=1001)

ax.set_title('3D case 1', pad=30)#@@@@@@@@@@@@@@@@@@@@@@@@@@@@
fig.savefig(r'.\3d_1.png', dpi=600, bbox_inches='tight')
plt.close(fig)


def draw_slice_axes(ax, y_val, n_ticks=5, tick_len_frac=0.02, fontsize=14):
    """Draw x and z axis lines, ticks and labels attached to the slice plane at y = y_val,
    instead of relying on mplot3d's default axes (which sit on the outer box, not the slice)."""
    def format_tick(value):
        return '0' if np.isclose(value, 0) else f'{value:.1f}'
    dz = tick_len_frac * Lz
    dx = tick_len_frac * Lx
    dy = tick_len_frac * Ly
    x_ticks = np.linspace(0, Lx, n_ticks)
    z_ticks = np.linspace(0, Lz, n_ticks)
    fontnewer=19

    # x-axis - runs along the bottom edge of the slice
    ax.plot([0, Lx], [y_val, y_val], [0, 0], color='black', linewidth=1, zorder=10)
    for tx in x_ticks:
        ax.plot([tx, tx], [y_val, y_val], [0, -dz], color='black', linewidth=1, zorder=10)
        ax.text(tx, y_val, -2 * dz, format_tick(tx), ha='center', va='top', fontsize=fontnewer)
    ax.text(Lx/2 - 4*dx, y_val - 8*dy, 0, '$x$ (m)', fontsize=fontnewer)

    # z-axis - runs along the left edge of the slice
    ax.plot([0, 0], [y_val, y_val], [0, Lz], color='black', linewidth=1, zorder=10)
    for tz in z_ticks:
        ax.plot([0, -dx], [y_val, y_val], [tz, tz], color='black', linewidth=1, zorder=10)
        ax.text(-2 * dx, y_val, tz, format_tick(tz), ha='right', va='center', fontsize=fontnewer)
    ax.text(-17*dx, y_val, Lz/2 + 4*dz, '$z$ (m)', fontsize=fontnewer)

def place_title_above_slice(fig, ax, y_val, text, gap=0.03, fontsize=14):
    """Put a title horizontally centred on the slice and just above its highest projected point.
    `gap` is in axes-fraction units."""
    fig.canvas.draw()
    M = ax.get_proj()
    corners = [(0, y_val, 0), (Lx, y_val, 0), (0, y_val, Lz), (Lx, y_val, Lz)]
    proj = [proj3d.proj_transform(cx, cy, cz, M)[:2] for cx, cy, cz in corners]
    disp = ax.transData.transform(proj)
    frac = ax.transAxes.inverted().transform(disp)
    x_centre = 0.5 * (frac[:, 0].min() + frac[:, 0].max())
    y_top = frac[:, 1].max()
    ax.text2D(x_centre, y_top + gap, text, transform=ax.transAxes, ha='center', va='bottom', fontsize=fontsize)

#Images 2-4 - xz-plane slices, all drawn at the same position
y_slice_indices = [round(slice_y_values[0]/(side/(divisions-1))), round(slice_y_values[1]/(side/(divisions-1))), round(slice_y_values[2]/(side/(divisions-1)))]######
labelcount=0
for n, yi in enumerate(y_slice_indices, start=2):
    fig, ax = new_iso_axes()
    ax.computed_zorder = False
    ax.set_axis_off()

    fig.patch.set_alpha(0.0)
    ax.patch.set_alpha(0.0)

    y_sampled = y_coords[yi]
    slice_data = points[:, yi, :]

    plot_face(ax, Xxz, np.full_like(Xxz, Y_DRAW), Zxz, slice_data)
    draw_slice_axes(ax, Y_DRAW)

    slice_colour = slice_colours[n - 2]
    ax.plot([0, 0], [Y_DRAW, Y_DRAW], [0, Lz], color=slice_colour, linestyle=':', linewidth=2*SLICE_LINE, zorder=1000)
    ax.plot([0, Lx], [Y_DRAW, Y_DRAW], [Lz, Lz], color=slice_colour, linestyle=':', linewidth=2*SLICE_LINE, zorder=1000)
    ax.plot([Lx, Lx], [Y_DRAW, Y_DRAW], [0, Lz], color=slice_colour, linestyle=':', linewidth=2*SLICE_LINE, zorder=1000)
    ax.plot([0, Lx], [Y_DRAW, Y_DRAW], [0, 0], color=slice_colour, linestyle=':', linewidth=2*SLICE_LINE, zorder=1000)

    theta = np.linspace(0, 2 * np.pi, 200)
    r = slice_r_values[labelcount]
    x = 0.8 + r * np.cos(theta)
    y = np.full_like(theta, Y_DRAW)
    z = 0.8 + r * np.sin(theta)
    ax.plot(x, y, z, color="black", linewidth=2)

    place_title_above_slice(fig, ax, Y_DRAW, f'Slice {slice_labels[labelcount]} ($y$ = {y_sampled:.2f} m)')
    labelcount += 1

    fig.savefig(f'3d_{n}.png', dpi=600, transparent=True)
    plt.close(fig)




