#! /usr/bin/env python

# gridgen functions for global 25km base grid

import gridgen as grd
import numpy as np
import matplotlib.pyplot as plt

# switches for steps in process
gentier = False
combtier = False
gridmk3 = True

# set up a dummy bathymtery for testing
bathyfile = '/project/ofrd/bathymetry/gebco_reduced_6.nc'
writedir = '/project/ofrd/waves/wavegrids/global/grid_netCDF'
project = 'global'

#----
# generate the SMC tier
if gentier:
    print('')
    print('*** Loading SMC base grid')
    basesmcrm = grd.loadNCsmc(writedir+'/global_GS512L2EUK.nc')
    # visualise the grid
    grd.plotGridsmc(basesmcrm,latlon=False)

    print('')
    print('*** Creating tier for SMC grid')
    # set limit for dry cells below zero - to remove very shallow coastal cells from model 
    tier1smc = grd.createTiersmc(bathyfile, basesmcrm,
                                  label='GS512L3EUK-sub',
                                  mindepth=None, dlim=3.0, drymin=None, drymax=0.7, 
                                  bathytype='gebco', getpland=True)
    # save the grid
    tier1file = tier1smc.writeNC(writedir=writedir)
    print('Data written to %s' %tier1file)
    # visualise the grid
    grd.plotGridsmc(tier1smc)

#---
# combined the basegrid and tier
if combtier:
    print('')
    print('*** Loading SMC base grid')
    basesmcrm = grd.loadNCsmc(writedir+'/global_GS512L2EUK.nc')
    # visualise the grid
    grd.plotGridsmc(basesmcrm,latlon=False)
    print('')
    print('*** Loading SMC tier')
    tier1smc = grd.loadNCsmc(writedir+'/global_GS512L3EUK-sub.nc')
    # visualise the grid
    grd.plotGridsmc(tier1smc,latlon=False)

    print('')
    print('*** Combining tier and base')
    combbaset1 = grd.joinTiersmc(basesmcrm,tier1smc,bathyfile)
    combbaset1.label = 'GS512L3EUK-comb'

    # save the grid
    combbaset1file = combbaset1.writeNC(writedir=writedir)
    print('Data written to %s' %combbaset1file)
    # visualise the grid
    grd.plotGridsmc(combbaset1,latlon=False)

##-- load grid and mark additional tier cells in regions of special interest
if gridmk3:
    print('')
    print('*** Loading SMC grid')
    basesmcrm = grd.loadNCsmc(writedir+'/global_GS512L3EUK-comb.nc')
    # visualise the grid
    grd.plotGridsmc(basesmcrm)

    print('')
    print('*** Deselecting tier cells and removing dry cells in specified regions')
    basesmcrm.label = 'GS512L3EUK'
    basesmcrm.unmarkCells(marker='tier', box=[344.0,45.0,9.4,61.15], osbox=True, thruzero=True)
    basesmcrm.delDryCells(celltype='dry')
    basefile = basesmcrm.writeNC(writedir=writedir)
    print('Data written to %s' %basefile)
    grd.plotGridsmc(basesmcrm)

