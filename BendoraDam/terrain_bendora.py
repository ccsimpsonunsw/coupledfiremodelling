import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from numpy import savetxt
from math import pow
from scipy import interpolate

we = 9000.
sn = 9000.
sr = 8
t_dx = 30.
t_dy = 30.
f_dxy = 30.
h_dxy = 30.
z_dxy = 3.75
width = 2000.

t_nx = np.int32(we/t_dx)
t_ny = np.int32(sn/t_dy)
t = np.zeros([t_nx,t_ny])

h_nx = np.int32(we/h_dxy)
h_ny = np.int32(sn/h_dxy)
hgt = np.zeros([h_nx,h_ny])

z_nx = np.int32(we/z_dxy)
z_ny = np.int32(sn/z_dxy)
zsf = np.zeros([z_nx,z_ny])
ncat = np.zeros([z_nx,z_ny])

f_nx = np.int32(we/f_dxy)
f_ny = np.int32(sn/f_dxy)
fcat = np.zeros([f_nx,f_ny])

t_x = np.arange(0,t_dx*t_nx,t_dx)
t_y = np.arange(0,t_dy*t_ny,t_dy)
h_x = np.arange(0,we,h_dxy)
h_y = np.arange(0,sn,h_dxy)
z_x = np.arange(0,we,z_dxy)
z_y = np.arange(0,sn,z_dxy)

t = np.loadtxt('bendora_height.asc')
t = np.flipud(t)
#t[t<5] = 10000
base = np.amin(t)
#print (np.amax(t)-np.amin(t))
#t[t>9999] = base

splines1 = interpolate.RectBivariateSpline(t_x,t_y,t)
hgt = interpolate.RectBivariateSpline.__call__(splines1,h_x,h_y)
for i in range(0,h_nx):
  for j in range(0,h_ny):
    if hgt[j,i] < 1.:
      hgt[j,i] = 0.

ga = 40.
zone = np.rint(width/h_dxy)
xe1 = 55
xe2 = 220
ye1 = 80
ye2 = 245
for i in range(0,h_ny):
  for j in range(0,h_nx):
    if i <= ye1:
      if j <= xe1:
        tmp = np.sqrt((i-ye1)*(i-ye1)+(j-xe1)*(j-xe1))
      if j > xe1 and j < xe2:
        tmp = i-ye1
      if j >= xe2:
        tmp = np.sqrt((i-ye1)*(i-ye1)+(j-xe2)*(j-xe2))
    if i >= ye2:
      if j <= xe1:
        tmp = np.sqrt((i-ye2)*(i-ye2)+(j-xe1)*(j-xe1))
      if j > xe1 and j < xe2:
        tmp = i-ye2
      if j >= xe2:
        tmp = np.sqrt((i-ye2)*(i-ye2)+(j-xe2)*(j-xe2))           
    if j <= xe1 and i > ye1 and i < ye2:
      tmp = xe1-j
    if j >= xe2 and i > ye1 and i < ye2:
      tmp = j-xe2
    hgt[j,i] = hgt[j,i]*np.exp(tmp*tmp/(-2.*ga*ga))-1.0 
    if hgt[j,i] < base:
      hgt[j,i] = base

splines2 = interpolate.RectBivariateSpline(h_x,h_y,hgt)
zsf = interpolate.RectBivariateSpline.__call__(splines2,z_x,z_y)
for i in range(0,z_nx):
  for j in range(0,z_ny):
    if zsf[j,i] < base:
      zsf[j,i] = base

f_x = np.arange(0,we,f_dxy)
f_y = np.arange(0,sn,f_dxy)

n = np.loadtxt('bendora_landwater.asc')
n = np.flipud(n)
n[n>5] = 10
n[n<5] = 14
fcat = n

for i in range(0,f_nx):
  for j in range(0,f_ny):
    ncat[i*sr:i*sr+sr,j*sr:j*sr+sr] = fcat[i,j]

np.savetxt('terrain2d_hgt.bendora',hgt,fmt='%5u',delimiter=' ')
np.savetxt('terrain2d_zsf.bendora',zsf,fmt='%5u',delimiter=' ')
np.savetxt('input_fc.bendora',ncat,fmt='%4u',delimiter=' ')

#h_lev = np.arange(0.,3000.,20.)
#f_lev = np.arange(12.5,14.5,1.)
#plt.figure()
#plt.colorbar()
#plt.contourf(h_x,h_y,hgt,levels=h_lev)
#plt.contourf(f_x,f_y,fcat,levels=f_lev,alpha=0.3)
#plt.show()
