#! /usr/bin/env python

# gridgen functions for global 25km base grid

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
writedir = '/project/ofrd/waves/wavegrids/global/grid_netCDF'
project = 'global'

#---grid extents ---
# set up a regular grid and extents - this provides the SMC base
dx = 360.0 / 1024.0
dy = dx * 2.0 / 3.0

swlat = -342.0 * dy
swlon = 0.0 * dx
# ne corner should be spec'd 1 higher than anticipated regular
# array index, e.g. 0.0-1024.0 multipliers should return indices 0-1023, ie. 1024 cells 
nelat = 367.0 * dy
nelon = 1024.0 * dx
extents = [swlon, swlat, nelon, nelat]

#--- smc grid ---
# uses the regular grid above as its base (coarse)grid
if gensmc:
    print('')
    print('*** Creating SMC base grid')

    basesmc = grd.createBasesmc(bathyfile, extents, dx, dy, 
                   name='global', label='GS512base',
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
    basesmcrm = grd.loadNCsmc(writedir+'/global_GS512base.nc')
    # visualise the grid
    grd.plotGridsmc(basesmcrm)

    print('')
    print('*** Cropping the base grid to remove unwanted seas')
    basesmcrm.label = 'GS512crop'
    basesmcrm.markRegion(87.5, 41.5, 90.7, 43.9, marker='dry') # inland lake in central asia
    basesmcrm.markRegion(62.3, 41.2, 63.8, 42.8, marker='dry') # inland lake next to Caspian
    basesmcrm.markRegion(136.0, -31.9, 140.9, -26.7, marker='dry') # inland lake in Australia
    basesmcrm.markRegion(267.0, 42.4, 284.6, 49.4, marker='dry') # Great Lakes
    basesmcrm.markRegion(5.2, 33.4, 7.8, 35.0, marker='dry') # North African Lake1 (Libya??)
    basesmcrm.markRegion(23.90, 28.33, 31.58, 30.50, marker='dry') # North African Lake2 (Egypt??)
    basesmcrm.markRegion(35.2, 30.5, 36.1, 32.6, marker='dry') # North African Lake3 (Israel??)
    basesmcrm.delDryCells()
    basefile = basesmcrm.writeNC(writedir=writedir)
    print('Data written to %s' %basefile)
    grd.plotGridsmc(basesmcrm)

##-- load basegrid and mark additional tier cells in shallow water
if basegridmk12:
    print('')
    print('*** Loading SMC base grid')
    basesmcrm = grd.loadNCsmc(writedir+'/global_GS512crop.nc')
    # visualise the grid
    grd.plotGridsmc(basesmcrm)

    print('')
    print('*** Adding tier cells in shallow water')
    basesmcrm.label = 'GS512L1EUK'
    basesmcrm.markRegion(335.0, 27.0, 360.0, 68.0, marker='tier') # Euro region west of Greenwich
    basesmcrm.markRegion(0.0, 27.0, 42.0, 68.0, marker='tier') # Euro region east of Greenwich
    basefile = basesmcrm.writeNC(writedir=writedir)
    print('Data written to %s' %basefile)
    grd.plotGridsmc(basesmcrm)

#----
# load the base SMC tier
if loadbase:
    print('')
    print('*** Loading SMC cropped grid')
    basesmcrm = grd.loadNCsmc(writedir+'/global_GS512L1EUK.nc')
    # visualise the grid
    grd.plotGridsmc(basesmcrm,latlon=False)


