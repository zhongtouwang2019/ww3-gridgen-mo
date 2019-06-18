#! /usr/bin/env python

# gridgen functions for atlantic 25km base grid

import gridgen as grd
import numpy as np
import matplotlib.pyplot as plt

# switches for steps in process
gentier = False
combtier = False
gridmk6 = True
chkgrd = False

# set up a dummy bathymtery for testing
bathyfile = '/project/ofrd/bathymetry/gebco_reduced_6.nc'
writedir = '/project/ofrd/waves/wavegrids/atlantic/grid_netCDF'
project = 'atlantic'

#----
# generate the first SMC tier
if gentier:
    print('')
    print('*** Loading SMC base grid')
    basesmcrm = grd.loadNCsmc(writedir+'/atlantic_AS512L1EUK.nc')
    # visualise the grid
    grd.plotGridsmc(basesmcrm,latlon=False)

    print('')
    print('*** Creating tier for SMC grid')
    tier1smc = grd.createTiersmc(bathyfile, basesmcrm,
                                  label='AS512L2EUK-sub',
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
    basesmcrm = grd.loadNCsmc(writedir+'/atlantic_AS512L1EUK.nc')
    # visualise the grid
    grd.plotGridsmc(basesmcrm,latlon=False)
    print('')
    print('*** Loading SMC tier')
    tier1smc = grd.loadNCsmc(writedir+'/atlantic_AS512L2EUK-sub.nc')
    # visualise the grid
    grd.plotGridsmc(tier1smc,latlon=False)

    print('')
    print('*** Combining tier and base')
    combbaset1 = grd.joinTiersmc(basesmcrm,tier1smc,bathyfile)
    combbaset1.label = 'AS512L2EUK-comb'

    # save the grid
    combbaset1file = combbaset1.writeNC(writedir=writedir)
    print('Data written to %s' %combbaset1file)
    # visualise the grid
    grd.plotGridsmc(combbaset1,latlon=False)

##-- load grid and mark additional tier cells in shallow water
if gridmk6:
    print('')
    print('*** Loading SMC grid')
    basesmcrm = grd.loadNCsmc(writedir+'/atlantic_AS512L2EUK-comb.nc')
    # visualise the grid
    grd.plotGridsmc(basesmcrm)

    print('')
    print('*** Adding tier cells in specified regions')
    basesmcrm.label = 'AS512L2EUK'
    basesmcrm.markRegion(-18.0, 30.0, 14.8, 63.0, marker='tier', depthlim=320.0) # wide UK region - west corresponds to global Euro zone
    basefile = basesmcrm.writeNC(writedir=writedir)
    print('Data written to %s' %basefile)
    grd.plotGridsmc(basesmcrm)

##-- check grid
if chkgrd:
    print('')
    print('*** Loading SMC grid')
    basesmcrm = grd.loadNCsmc(writedir+'/atlantic_AS512L2EUK.nc')
    # visualise the grid
    grd.plotGridsmc(basesmcrm, latlon=False)
