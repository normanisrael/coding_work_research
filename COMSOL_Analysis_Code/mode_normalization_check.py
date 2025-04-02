import numpy as np
import matplotlib.pyplot as plt
import hamiltonian_functions as hf
from scipy.interpolate import griddata

'''
modeA = '/home/nanophotonics/Research/Project_Top/EparT.txt' #nanorodmode.txt
modeB = '/home/nanophotonics/Research/Project_Top/EparT.txt' #nanorodmode.txt
mesh = '/home/nanophotonics/Research/Project_Top/ER5nmmesh_1par.mphtxt' #nanorodmesh.mphtxt
QNAf =  '/home/nanophotonics/Research/Project_Top/5nmeigF_1par.txt' #nanorodeigenf.txt
QNBf =  '/home/nanophotonics/Research/Project_Top/5nmeigF_1par.txt' #nanorodeigenf.txt

modeA = '/home/nanophotonics/Research/Project_Top/E20nmmodes_1par.txt'#cubemode.txt'
modeB = '/home/nanophotonics/Research/Project_Top/E20nmmodes_1par.txt'#cubemode.txt'
mesh = '/home/nanophotonics/Research/Project_Top/ER20nmmesh_1par.mphtxt'#Esqnmmesh_1par.mphtxt'
QNAf = '/home/nanophotonics/Research/Project_Top/E20nmeigenf.txt'#cubeeigenf.txt'
QNBf = '/home/nanophotonics/Research/Project_Top/E20nmeigenf.txt'#cubeeigenf.txt'

modeA = '/home/nanophotonics/Research/Project_Top/cubemode.txt'
modeB = '/home/nanophotonics/Research/Project_Top/cubemode.txt'
mesh = '/home/nanophotonics/Research/Project_Top/Esqnmmesh_1par.mphtxt'
QNAf = '/home/nanophotonics/Research/Project_Top/cubeeigenf.txt'
QNBf = '/home/nanophotonics/Research/Project_Top/cubeeigenf.txt'
'''

modeA = '/home/nanophotonics/Research/Project_Top/nanorodmode3.txt'
modeB = '/home/nanophotonics/Research/Project_Top/nanorodmode3.txt'
mesh = '/home/nanophotonics/Research/Project_Top/nanorodmesh2.mphtxt'
QNAf =  '/home/nanophotonics/Research/Project_Top/nrodeigf.txt'
QNBf =  '/home/nanophotonics/Research/Project_Top/nrodeigf.txt'

mesh_coordinates = np.loadtxt(mesh,dtype=float,skiprows=1)
EA = np.loadtxt(modeA,dtype=str,skiprows=9)
EB = np.loadtxt(modeB,dtype=str,skiprows=9)
QNA = np.loadtxt(QNAf,dtype=str,skiprows=5)
QNB = np.loadtxt(QNBf,dtype=str,skiprows=5)
eps0 = 8.85E-12

xmesh = mesh_coordinates[:,0]
ymesh = mesh_coordinates[:,1]
zmesh = mesh_coordinates[:,2]

xA_mode_coordinate = EA[:,0].astype(float)
yA_mode_coordinate = EA[:,1].astype(float)
zA_mode_coordinate = EA[:,2].astype(float)

xB_mode_coordinate = EB[:,0].astype(float)
yB_mode_coordinate = EB[:,1].astype(float)
zB_mode_coordinate = EB[:,2].astype(float)

modeA_no = 1
indxA = 3*modeA_no
ExA = np.array(hf.i2j(EA[:,indxA])).astype(complex)
EyA = np.array(hf.i2j(EA[:,indxA + 1])).astype(complex)
EzA = np.array(hf.i2j(EA[:,indxA + 2])).astype(complex)

modeB_no = 1
indxB = 3*modeB_no
ExB = np.array(hf.i2j(EB[:,indxB])).astype(complex)
EyB = np.array(hf.i2j(EB[:,indxB + 1])).astype(complex)
EzB = np.array(hf.i2j(EB[:,indxB + 2])).astype(complex)

fA = np.array(hf.i2j(QNA[:,0])).astype(complex)
fB = np.array(hf.i2j(QNB[:,0])).astype(complex)

QNNA = np.array(hf.i2j(QNA[:,3])).astype(complex)/eps0
QNNB = np.array(hf.i2j(QNB[:,3])).astype(complex)/eps0

eps_inf = 1.0#9.1
omp = 1.26E16#1.38E16
gam = 1.41E14#1.08E14
om0 = 0.0
OmegaA = 2*np.pi*fA[modeA_no-1] 
OmegaB = 2*np.pi*fB[modeB_no-1]

dpdw = hf.depsdw(OmegaA,eps_inf*eps0,omp,om0,gam)

nm = 1e-9

no_points = 4

center_of_region = np.array([0,0,0])*nm
radius_of_region = 100.0*nm 

region_x_coordinate, region_y_coordinate, region_z_coordinate = hf.region(xmesh,ymesh,zmesh,center_of_region,radius_of_region)

n, w = np.polynomial.legendre.leggauss(no_points)
x, y, z = np.meshgrid(n,n,n)
wx, wy, wz = np.meshgrid(w,w,w)
x, y, z = x.flatten(), y.flatten(), z.flatten()
wx, wy, wz = wx.flatten(), wy.flatten(), wz.flatten()

#All space
ux = xA_mode_coordinate.max()
lx = xA_mode_coordinate.min()
uy = yA_mode_coordinate.max()
ly = yA_mode_coordinate.min()
uz = zA_mode_coordinate.max()
lz = zA_mode_coordinate.min()

x_sample = (ux - lx)/2 * np.array(x) + (ux + lx)/2 
y_sample = (uy - ly)/2 * np.array(y) + (uy + ly)/2
z_sample = (uz - lz)/2 * np.array(z) + (uz + lz)/2

interp_method = 'nearest'

xxa = x_sample
yya = y_sample
zza = z_sample

xxb = x_sample
yyb = y_sample
zzb = z_sample 

phiEAx = griddata((xA_mode_coordinate, yA_mode_coordinate, zA_mode_coordinate),ExA,(xxa,yya,zza),method=interp_method)
phiEAy = griddata((xA_mode_coordinate, yA_mode_coordinate, zA_mode_coordinate),EyA,(xxa,yya,zza),method=interp_method)
phiEAz = griddata((xA_mode_coordinate, yA_mode_coordinate, zA_mode_coordinate),EzA,(xxa,yya,zza),method=interp_method)

phiEBx = griddata((xB_mode_coordinate, yB_mode_coordinate, zB_mode_coordinate),ExB,(xxb,yyb,zzb),method=interp_method)
phiEBy = griddata((xB_mode_coordinate, yB_mode_coordinate, zB_mode_coordinate),EyB,(xxb,yyb,zzb),method=interp_method)
phiEBz = griddata((xB_mode_coordinate, yB_mode_coordinate, zB_mode_coordinate),EzB,(xxb,yyb,zzb),method=interp_method)

inXA = hf.gaussian_int3d(np.conjugate(phiEAx)*phiEAx,wx,wy,wz,ux,lx,uy,ly,uz,lz)
inYA = hf.gaussian_int3d(np.conjugate(phiEAy)*phiEAy,wx,wy,wz,ux,lx,uy,ly,uz,lz)
inZA = hf.gaussian_int3d(np.conjugate(phiEAz)*phiEAz,wx,wy,wz,ux,lx,uy,ly,uz,lz)

inXB = hf.gaussian_int3d(np.conjugate(phiEBx)*phiEBx,wx,wy,wz,ux,lx,uy,ly,uz,lz)
inYB = hf.gaussian_int3d(np.conjugate(phiEBy)*phiEBy,wx,wy,wz,ux,lx,uy,ly,uz,lz)
inZB = hf.gaussian_int3d(np.conjugate(phiEBz)*phiEBz,wx,wy,wz,ux,lx,uy,ly,uz,lz)

NA = (inXA + inYA + inZA)*dpdw
NB = (inXB + inYB + inZB)*dpdw

exA = phiEAx/np.sqrt(QNNA[modeA_no-1])
eyA = phiEAy/np.sqrt(QNNA[modeA_no-1])
ezA = phiEAz/np.sqrt(QNNA[modeA_no-1])

exB = phiEBx/np.sqrt(QNNA[modeA_no-1])
eyB = phiEBy/np.sqrt(QNNA[modeA_no-1])
ezB = phiEBz/np.sqrt(QNNA[modeA_no-1])

innerProdX = hf.gaussian_int3d(np.conjugate(exA)*exB,wx,wy,wz,ux,lx,uy,ly,uz,lz)
innerProdY = hf.gaussian_int3d(np.conjugate(eyA)*eyB,wx,wy,wz,ux,lx,uy,ly,uz,lz)
innerProdZ = hf.gaussian_int3d(np.conjugate(ezA)*ezB,wx,wy,wz,ux,lx,uy,ly,uz,lz)

innerProd = (innerProdX + innerProdY + innerProdZ)#*dpdw
print(innerProd)
