#! /usr/bin/env python

# gridgen functions for global 25km base grid

import gridgen as grd
import numpy as np
import matplotlib.pyplot as plt

# switches for steps in process
gentier = False
combtier = True
gridmk6 = False
chkgrd = False

# set up a dummy bathymtery for testing
bathyfile = '/project/ofrd/bathymetry/gebco_reduced_6.nc'
writedir = '/project/ofrd/waves/wavegrids/global/grid_netCDF'
project = 'global'

#----
# generate the first SMC tier
if gentier:
    print('')
    print('*** Loading SMC base grid')
    basesmcrm = grd.loadNCsmc(writedir+'/global_GS512L1.nc')
    # visualise the grid
    grd.plotGridsmc(basesmcrm,latlon=False)

    print('')
    print('*** Creating tier for SMC grid')
    tier1smc = grd.createTiersmc(bathyfile, basesmcrm,
                                  label='GS512L2-sub',
                                  mindepth=None, dlim=3.0, drymin=None, drymax=0.7, 
                                  bathytype='gebco', getpland=True)
    # save the grid
    tier1file = tier1smc.writeNC(writedir=writedir)
    print('Data written to %s' %tier1file)
    # visualise the grid
    grd.plotGridsmc(tier1smc)

#---
# combined the (cropped) basegrid and tier1
if combtier:
    print('')
    print('*** Loading SMC base grid')
    basesmcrm = grd.loadNCsmc(writedir+'/global_GS512L1.nc')
    # visualise the grid
    grd.plotGridsmc(basesmcrm,latlon=False)
    print('')
    print('*** Loading SMC tier')
    tier1smc = grd.loadNCsmc(writedir+'/global_GS512L2-sub.nc')
    # visualise the grid
    grd.plotGridsmc(tier1smc,latlon=False)

    print('')
    print('*** Combining tier and base')
    combbaset1 = grd.joinTiersmc(basesmcrm,tier1smc,bathyfile)
    combbaset1.label = 'GS512L2-comb'

    # save the grid
    combbaset1file = combbaset1.writeNC(writedir=writedir)
    print('Data written to %s' %combbaset1file)
    # visualise the grid
    grd.plotGridsmc(combbaset1,latlon=False)

##-- load grid and mark additional tier cells in shallow water
if gridmk6:
    print('')
    print('*** Loading SMC grid')
    basesmcrm = grd.loadNCsmc(writedir+'/global_GS512L2-comb.nc')
    # visualise the grid
    grd.plotGridsmc(basesmcrm)

    print('')
    print('*** Adding tier cells in specified regions')
    basesmcrm.label = 'GS512L2'
    #basesmcrm.markDepths(20.0, marker='tier') # global 6km cells active at <10m
    #basesmcrm.markRegion(330.0, 25.0, 360.0, 68.0, marker='tier', depthlim=80.0) # Euro region west of Greenwich
    #basesmcrm.markRegion(0.0, 25.0, 42.0, 68.0, marker='tier', depthlim=80.0) # Euro region east of Greenwich
    basefile = basesmcrm.writeNC(writedir=writedir)
    print('Data written to %s' %basefile)
    grd.plotGridsmc(basesmcrm)

##-- check grid
if chkgrd:
    print('')
    print('*** Loading SMC grid')
    basesmcrm = grd.loadNCsmc(writedir+'/global_GS512L2-comb.nc')
    # visualise the grid
    grd.plotGridsmc(basesmcrm, latlon=False)
