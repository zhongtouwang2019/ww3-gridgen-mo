#! /usr/bin/env python

# gridgen functions for global 25km base grid

import gridgen as grd
import numpy as np
import matplotlib.pyplot as plt

# switches for steps in process
gentier = False
combtier = False
writeww3 = True

# set up a dummy bathymtery for testing
bathyfile = '/project/ofrd/bathymetry/gebco_reduced_6.nc'
writedir = '/project/ofrd/waves/wavegrids/global/grid_netCDF'
project = 'global'

#----
# generate the SMC tier
if gentier:
    print('')
    print('*** Loading SMC base grid')
    basesmcrm = grd.loadNCsmc(writedir+'/global_GS512L2-comb.nc')
    # visualise the grid
    grd.plotGridsmc(basesmcrm,latlon=False)

    print('')
    print('*** Creating tier for SMC grid')
    # set limit for dry cells below zero - to remove very shallow coastal cells from model 
    tier1smc = grd.createTiersmc(bathyfile, basesmcrm,
                                  label='GS512L3-sub',
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
    basesmcrm = grd.loadNCsmc(writedir+'/global_GS512L2-comb.nc')
    # visualise the grid
    grd.plotGridsmc(basesmcrm,latlon=False)
    print('')
    print('*** Loading SMC tier')
    tier1smc = grd.loadNCsmc(writedir+'/global_GS512L3-sub.nc')
    # visualise the grid
    grd.plotGridsmc(tier1smc,latlon=False)

    print('')
    print('*** Combining tier and base')
    combbaset1 = grd.joinTiersmc(basesmcrm,tier1smc,bathyfile)
    combbaset1.label = 'GS512L3-comb'

    # save the grid
    combbaset1file = combbaset1.writeNC(writedir=writedir)
    print('Data written to %s' %combbaset1file)
    # visualise the grid
    grd.plotGridsmc(combbaset1,latlon=False)

if writeww3:
    print('')
    print('*** Loading combined grid')
    combsmc = grd.loadNCsmc(writedir+'/global_GS512L3-comb.nc')
    # visualize the grid
    grd.plotGridsmc(combsmc, latlon=False)

    print('')
    print('*** Removing dry cells')
    # cell removal here runs some risk of putting higher tier cells by coast??
    #pcdry = 0.5
    #combsmc.markDrypc(pcdry)
    #combsmc.markDepths(5.0,marker='dry')
    #combsmc.unmarkCells(marker='tier')
    combsmc.delDryCells(celltype='alldry')
    # visualize the grid
    grd.plotGridsmc(combsmc, latlon=False)

    print('')
    print('*** Writing smc cells to text file')
    combsmc.depthmin = 15.
    combsmc.sortCells()
    #combsmc.writeWW3(writedir=writedir+'/'+project, mindepth=combsmc.depthmin, writemindepth=True)
    combsmc.writeWW3(writedir='/project/ofrd/waves/wavegrids/global/GS512L3')
