import os
import glob
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('classic')
plt.rc('figure',facecolor='w')



file_cloud_end = '/Users/oanavesa/Desktop/GradSchool/FirstYear/ASTRO506/data_headon/data_t00000300.txt'
xc,yc,zc,vxc,vyc,vzc = np.genfromtxt(file_cloud_end, usecols=(0,1,2,3,4,5), unpack=True, skip_header=2)


rv = []
vv = []
angxl = []
mass_of_particles = 1/40000
with open('ang2.txt', 'a') as h:
    for i in range(0,len(xc)):
        Lx = yc[i]*mass_of_particles*vzc[i] - zc[i]*mass_of_particles*vyc[i]
        Ly = zc[i]*mass_of_particles*vxc[i] - xc[i]*mass_of_particles*vzc[i]
        Lz = xc[i]*mass_of_particles*vyc[i] - yc[i]*mass_of_particles*vxc[i]


            #r_vectory = yc[j]-yc[i]
            #r_vectorz = zc[j]-zc[i]
            #velocity_vectory = vyc[j]-vyc[i]
            #velocity_vectorz = vzc[j]-vzc[i]
            #angular_momentumx = np.cross(r_vectorx,mass_of_particles*velocity_vectorx)
            #summa = np.sum(angular_momentumx)
            #agx.append(summa)
            #rv.append(Lz)
        h.writelines(np.transpose([str(Lx) + "  " + str(Ly) + "  " + str(Lz) + '\n']))
h.close()
