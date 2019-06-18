#! /usr/bin/env python

# gridgen functions for atlantic 25km base grid

import gridgen as grd
import numpy as np
import matplotlib.pyplot as plt

# switches for steps in process
gensmc = False
basegridcrop = False
basegridmk12 = True
loadbase = False

# set up a dummy bathymtery for testing
bathyfile = '/project/ofrd/bathymetry/gebco_reduced_6.nc'
writedir = '/project/ofrd/waves/wavegrids/atlantic/grid_netCDF'
project = 'atlantic'

#---grid extents ---
# set up a regular grid and extents - this provides the SMC base
dx = 360.0 / 1024.0
dy = dx * 2.0 / 3.0

# ne corner should be spec'd 1 higher than anticipated regular
# array index, e.g. 0.0-1024.0 multipliers should return indices 0-1023, ie. 1024 cells 
swlat = -104.0 * dy
swlon = -280.0 * dx
nelat = 350.0 * dy
nelon = 191.0 * dx
extents = [swlon, swlat, nelon, nelat]

#--- smc grid ---
# uses the regular grid above as its base (coarse)grid
if gensmc:
    print('')
    print('*** Creating SMC base grid')

    basesmc = grd.createBasesmc(bathyfile, extents, dx, dy, 
                   name='atlantic', label='AS512base',
                   mindepth=None, dlim=3.0, drymin=None, drymax=0.7, 
                   bathytype='gebco', getpland=True)
				   
    # save the grid
    basefile = basesmc.writeNC(writedir=writedir)
    print('Data written to %s' %basefile)

    # visualise the grid
    grd.plotGridsmc(basesmc, latlon=False)

##-- load basegrid and remove cells from regions we don't want to compute
if basegridcrop:
    print('')
    print('*** Loading SMC base grid')
    basesmcrm = grd.loadNCsmc(writedir+'/atlantic_AS512base.nc')
    # visualise the grid
    grd.plotGridsmc(basesmcrm)

    print('')
    print('*** Cropping the base grid to remove unwanted seas')
    basesmcrm.label = 'AS512crop'
    basesmcrm.markRegion(18., -24.4, 67.4, 50.47, marker='dry') # remove e-med, w-indian ocean and gulf
    basesmcrm.markRegion(1.5, 30.3, 18., 46.1, marker='dry') # remove w-med (1st section) 
    basesmcrm.markRegion(-5.2, 34.75, 1.5, 41.33, marker='dry') # remove w-med (2nd section) 
    basesmcrm.markRegion(-100., -24.4, -63.75, 7.12, marker='dry') # remove e-pacific (1st section) 
    basesmcrm.markRegion(-99., 6.85, -89.3, 16.75, marker='dry') # remove e-pacific (2nd section) 
    basesmcrm.markRegion(-90.2, 6.85, -84.32, 14.13, marker='dry') # remove e-pacific (3rd section) 
    basesmcrm.markRegion(-84.38, 7., -77.70, 8.65, marker='dry') # remove e-pacific (4th section) 
    basesmcrm.markRegion(-92.2, 42.85, -76.1, 49.0, marker='dry') # remove great lakes 
    basesmcrm.markRegion(-99.2, 50.52, -70.0, 69.9, marker='dry') # remove Hudson Bay 
    basesmcrm.markRegion(-99.2, 65.80, -80.0, 83.0, marker='dry') # remove NW passage 
    basesmcrm.markRegion(-57.4, 80.40, -30.0, 83.0, marker='dry') # remove coastal sea N of Greenland
    basesmcrm.delDryCells()
    basefile = basesmcrm.writeNC(writedir=writedir)
    print('Data written to %s' %basefile)
    grd.plotGridsmc(basesmcrm)

##-- load basegrid and mark additional tier cells in shallow water
if basegridmk12:
    print('')
    print('*** Loading SMC base grid')
    basesmcrm = grd.loadNCsmc(writedir+'/atlantic_AS512crop.nc')
    # visualise the grid
    grd.plotGridsmc(basesmcrm)

    print('')
    print('*** Adding tier cells in shallow water')
    basesmcrm.label = 'AS512L1EUK'
    basesmcrm.markRegion(-25.0, 27.0, 12.5, 68.0, marker='tier') # wide UK region - west corresponds to global Euro zone
    basefile = basesmcrm.writeNC(writedir=writedir)
    print('Data written to %s' %basefile)
    grd.plotGridsmc(basesmcrm)

#----
# load the base SMC tier
if loadbase:
    print('')
    print('*** Loading SMC cropped grid')
    basesmcrm = grd.loadNCsmc(writedir+'/atlantic_AS512L1EUK.nc')
    # visualise the grid
    grd.plotGridsmc(basesmcrm,latlon=False)


