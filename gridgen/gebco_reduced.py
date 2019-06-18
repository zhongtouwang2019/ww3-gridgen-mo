#! /usr/bin/env python

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

pltchk = False

scalefac = 6

hiresfile = 'GEBCO_2014_2D.nc'
loresfile = 'gebco_reduced_%d' %scalefac + '.nc'

d = nc.Dataset(hiresfile)

dimy = 21600
dimx = 43200
offsx = 0
offsy = 0

#dimy = 1200
#dimx = 2400
#offsx = 3200
#offsy = 8400

deepmin = -702.

depout = np.zeros([dimy/scalefac, dimx/scalefac])
mskout = np.zeros([dimy/scalefac, dimx/scalefac])
latout = np.zeros(dimy/scalefac)
lonout = np.zeros(dimx/scalefac)

lptot = np.float(scalefac**2)
for lpx in range(scalefac):
    for lpy in range(scalefac):
        lpiter = lpx * scalefac + lpy
        print('Working on iteration %d' %lpiter + ' out of %d' %lptot)
        dep = d.variables['elevation'][lpy+offsy:dimy+offsy:scalefac,lpx+offsx:dimx+offsx:scalefac]
        #dep[dep < deepmin] = deepmin
        msk = np.zeros(np.shape(dep))
        msk[dep>0.] = 1.
        lat = d.variables['lat'][lpy+offsy:dimy+offsy:scalefac]
        lon = d.variables['lon'][lpx+offsx:dimx+offsx:scalefac]
        depout = depout + dep / lptot
        mskout = mskout + msk / lptot
        lonout = lonout + lon / lptot
        latout = latout + lat / lptot

if pltchk:
    plt.subplot(2,1,1)
    depouttmp = np.ma.masked_greater(depout, 5.0)
    plt.pcolormesh(lonout,latout,depouttmp,vmin=-500.)
    plt.colorbar()
    plt.subplot(2,1,2)
    plt.pcolormesh(lonout,latout,mskout)
    plt.colorbar()
    plt.show()

d.close()

# write out to a new netCDF file

with nc.Dataset(loresfile, 'w') as nbg:
    ndimx = nbg.createDimension('lon',size=dimx/scalefac)
    ndimy = nbg.createDimension('lat',size=dimy/scalefac)

    ndep = nbg.createVariable('lat','f8',dimensions=('lat'))
    ndep.units = 'degrees_east'
    ndep[:] = latout[:]

    ndep = nbg.createVariable('lon','f8',dimensions=('lon'))
    ndep.units = 'degrees_north'
    ndep[:] = lonout[:]

    ndep = nbg.createVariable('elevation','f8',dimensions=('lat','lon'))
    ndep.units = 'm'
    ndep[:,:] = depout[:,:]

    ndep = nbg.createVariable('landmask','f8',dimensions=('lat','lon'))
    ndep.units = '1'
    ndep[:,:] = mskout[:,:]

