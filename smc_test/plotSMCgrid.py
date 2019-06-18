"""
##
##  It reads cell files and uses own projection to draw 
##  the SMC grid. Projected polygons are collected into 
##  vert variables for grid and subsequent swh plots. 
##
#;  First created: For 4 views     J G Li   26 Nov 2008
#;  Modified for Color grids       J G Li   26 Aug 2009
#;  Adapted for 25km SMC grids     J G Li    3 Feb 2010
#;  Updated G25SMC-12c grids.      J G Li   25 May 2011
#;  Sterographic projection nR.    J G Li   30 Jun 2011
#;  Adapted for G50SMC grid.       J G Li   19 Aug 2011
#;  Extended to fill the Arctic.   J G Li    5 Oct 2011
#;  Rectify polar cell position.   J G Li   25 Oct 2011
#;  Simplify with readcell and steromap.  JGLi12Nov2014
##
##  Converted into a Python function.     JGLi05Dec2018
##  Save ELat/Lon and sx/yc in file.      JGLi11Dec2018
##  Add color map and draw color bar.     JGLi14Dec2018
##  Adapted for SMC36125 grid plot.       JGLi03Jan2019
##  Import resource and set stacksize.    JGLi07Jan2019
##  Use polycollections for two plots.    JGLi30Jan2019
##  Adapted for SMC61250 global grid.     JGLi18Feb2019
##  Adapted for SMC61250 global grid.     JGLi18Feb2019
##  Adapted for Andy's SMC6125  grid.     JGLi12Mar2019
##
"""

##  Import relevant modules and functions

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from datetime import datetime
import sys

# Jian-Guo's local codes for plotting
MCodes='/home/h05/frjl/Python/MyCodes/'
sys.path.append(MCodes)
from readcell import readcell   
from steromap import steromap
from rgbcolor import rgbcolor

def pltSMCgrid(ModlName, nlevs, dx, dy, x0, y0, nx, ny,
                Wrkdir, Cel_file='ww3Cels.dat', Arc_file=None, DatGMC=None,
                figtype='.ps'):

    print(" Program started at %s " % datetime.now().strftime('%F %H:%M:%S'))

##  Read global and Arctic part cells. 
    if DatGMC is None:
        DatGMC=Wrkdir
    Cel_file = DatGMC+'/'+Cel_file
    if Arc_file is not None:
        Arc_file = DatGMC+'/'+Cel_file

    if Arc_file is not None:
        headrs, cel = readcell( [Cel_file, Arc_file] ) 
        na = int( headrs[1].split()[0] )
        nb = int( headrs[1].split()[1] )
    else:
        headrs, cel = readcell([Cel_file]) 
        na = 0
        nb = 0
    ng = int( headrs[0].split()[0] )
    nc = ng + na
    print(' Merged total cel number = %d' %nc)

##  Size-1 cell increments
    dxlon = dx / 2.0**(nlevs-1)
    dylat = dy / 2.0**(nlevs-1)

##  Origin and cells numbers for size-1 cell i j indices.
    #xlon0 = x0 - 2.0*dxlon
    #ylat0 = y0 - 2.0*dylat
    xlon0 = x0 - dx/2.0
    ylat0 = y0 - dy/2.0
    nx = (nx+1) * 2**(nlevs-1) + 1
    ny = (ny+1) * 2**(nlevs-1) + 1

##  Half the latitude grid to set j=0 on the Equator.
#   neqt=(ny-1)/2
#   xlon=(np.arange(nx))*dxlon
#   ylat=(np.arange(ny)-neqt)*dylat
    neqt=0
    xlon=(np.arange(nx))*dxlon + xlon0
    ylat=(np.arange(ny))*dylat+ylat0
    print(' Full and Half lat grid = %d %d' %(ny-1, neqt))

##  Adjust cel[:,1] by neqt for easy indexing in ylat
#   cel[:,1]=cel[:,1]+neqt

##  Maximum j row number in Global part
    imx = cel[:,0].max()
    jmx = cel[:,1].max()
    dimx= cel[:,2].max()
    djmx= cel[:,3].max()
    khmx= cel[:,4].max()
    print(' Maximum i j di dj k = %d %d %d %d %d' %(imx, jmx, dimx, djmx, khmx))
    khmn= cel[:,3].min()
    print(' Minimum depth in cel = %.1f' %khmn)

##  Extra array in config for Arctic part
#   Arctic = True
    ngabjm = [ng, na, nb, jmx]
    if Arc_file is None:
        Arctic = False

##  Use own color map and defined depth colors 
    colrfile = MCodes+'rgbspectrum.dat'
    colrs = rgbcolor( colrfile )

##  Maximum mapping radius.
    radius=10.0

##  Possible selection of your plot types. 
    gorloc={0:'Global',1:'EuroArc',2:'Pacific',3:'Atlantic'}

##  Prompt selection choices and ask for one input
    print(gorloc)
    instr = input(' *** Please enter your selected number here > ')
    m = int(instr)
    pltype=gorloc.get(m, 'Invalid_selection')
    if( pltype == 'Invalid_selection' ): 
        print("Invalid selection, program terminated.")
        exit()

    print(" Draw SMC grid "+pltype)

    if( pltype == 'Global'):
##  Whole global projection angle from N Pole to be 90.0
        pangle=90.0
        plon= 0.0 
        plat= 23.5 
        clrbxy=[ -9.6,-12.6, 19.0,  1.0]
        sztpxy=[ 16.0, 10.0,-10.2,-10.3]
        rngsxy=[-10.0, 10.0,-13.0, 10.0]
        from smcglobl import smcglobl

    if( pltype == 'EuroArc'):
##  Euro-Arctic regional plot
        pangle=27.5 
        plon=  4.0
        plat= 69.0
        clrbxy=[  7.5, -4.5,   0.8,  7.0]
        sztpxy=[ 11.0, 15.0,   4.0, -5.0]
        rngsxy=[-10.0, 10.0, -14.0, 14.0]
        papror='portrait'

    if( pltype == 'Pacific'):
##  West Pacific regional plot
        plon= 138.0
        plat=  18.0
        pangle=33.0 
        clrbxy=[-13.6,  7.5,   7.5,  0.8]
        sztpxy=[ 15.0, 10.0, -10.0,  7.0]
        rngsxy=[-15.0, 15.0, -10.0, 10.0]
        papror='landscape'

    if( pltype == 'Atlantic'):
##  Whole global projection angle from N Pole to be 90.0
        pangle=64.0
        plon= 320.0 
        plat= 23.6 
        clrbxy=[ -9.6, -9.9, 19.0,  0.8]
        sztpxy=[ 10.0, 12.0, -6.4, -4.6]
        rngsxy=[-10.0, 10.0,-10.0, 11.0]
        papror='portrait'

#   print " Start loop over cells ... "
    print( " Start loop over cells at %s " % datetime.now().strftime('%F %H:%M:%S') )

##  Initial verts and ncels variable for polycollections.
    nvrts = []
    ncels = []
    svrts = []
    scels = []
##  Exclude last cell, the North Polar Cell, and Arctic boundary cells.
    for i in range(nc): 
        if( (i < ng-nb) or (i >= ng+nb) ):
            xc=[cel[i,0],cel[i,0]+cel[i,2],cel[i,0]+cel[i,2],cel[i,0]]
            yc=[cel[i,1],cel[i,1],cel[i,3]+cel[i,1],cel[i,3]+cel[i,1]]
            slat=ylat[yc]
            slon=xlon[xc]

##  Convert slat slon to elat elon with given new pole
            elat,elon,sxc,syc = steromap(slat,slon,plat,plon,Pangl=pangle,Onecl=True)

            if( (elat[0] >= 0.0) and (rngsxy[0] < sxc[0] < rngsxy[1])
                                 and (rngsxy[2] < syc[0] < rngsxy[3]) ):
                nvrts.append(list(zip(sxc,syc)))
                ncels.append(i)

            if( (elat[0] <  0.0) and (pltype == 'Global') ):
                svrts.append(list(zip(sxc,syc)))
                scels.append(i)

##  End of cell i loop excluding polar cell
    print( " Processing polar cell at %s " % datetime.now().strftime('%F %H:%M:%S') )
#;  Polar cell as a octagon, note polar cell size set as size-64 for flux calculation
#;  As it mergers 8 size-64 cells, each side of the octagon should be size-64 
#;  10 apexes are calculated to conform with other cells.  To plot it, drop the last apex.
##  Use a square box for polar cell to fit the new array size.  JGLi30Jan2019
#   i=nc-1
#   xc=cel[i,0]+np.arange(4)*(nx-1)/4
#   yc=xc*0+cel[i,1]
#   slat=ylat[yc]
#   slon=xlon[xc]

#;  Convert slat slon to elat elon with given new pole
#   elat, elon, sxc, syc = steromap( slat, slon, plat, plon, Onecl=True )

#   if( (elat[0] >= 0.0) and (rngsxy[0] < sxc[0] < rngsxy[1])
#                        and (rngsxy[2] < syc[0] < rngsxy[3]) ):
#       nvrts.append(list(zip(sxc,syc)))
#       ncels.append(i)

#   if( (elat[0] <  0.0) and (pltype == 'Global') ):
#       svrts.append(list(zip(sxc,syc)))
#       scels.append(i)

##  Set plot size and limits and message out anchor point.
    rdpols=[radius, pangle, plon, plat]
    config=np.array([rdpols, sztpxy, rngsxy, clrbxy, ngabjm])
    pzfile=DatGMC+'/' + ModlName + '_Vrts'+pltype[0:4]+'.npz'

##  Store selected north and south verts and cell numbers for swh plots.
##  Use the np.savez to save 3/5 variables in one file.  JGLi22Feb2019 
    if( pltype == 'Global' ):
        np.savez( pzfile, nvrt=nvrts, ncel=ncels, cnfg=config, 
                          svrt=svrts, scel=scels)
    else:
        np.savez( pzfile, nvrt=nvrts, ncel=ncels, cnfg=config) 

##  These variables could be loaded back by
#   vrtcls = np.load(DatGMC+'S625Vrts'+pltype[0:4]+'.npz')
#   nvrts = vrtcls['nvrt'] ; ncels = vrtcls['ncel']; config=vrtcls['cnfg']; 
#   svrts = vrtcls['svrt'] ; scels = vrtcls['scel']
##

##  Draw your selected grid plot.
    psfile=Wrkdir + '/' + ModlName + '_' + pltype[0:4] + 'grd' + figtype 
    bufile=MCodes + '/ECBuoys.dat'

    if( pltype == 'Global'):
        from smcglobl import smcglobl
        smcglobl( cel, nvrts,ncels,svrts,scels,colrs, config,
             mdlname= ModlName, Arctic=Arctic, psfile=psfile)

    else:
        from smclocal import smclocal
        smclocal( cel, nvrts,ncels,colrs,config, Arctic=Arctic, 
              mdlname= ModlName,  buoys=bufile, psfile=psfile,
              paprorn=papror)

    print( " Program finished at %s " % datetime.now().strftime('%F %H:%M:%S') )

## End of smc6125grids program ##
    return


#---
# main program
#Wrkdir='/project/ofrd/waves/wavegrids/atlsmc/A36125'
#ModlName='A36125'
#nlevs=4
## regular grid info at coarsest level
#dx = 360.0 / 1024.0
#dy = dx * 2.0 / 3.0
#x0 = -98.261719
#y0 = -24.257812
#nx = 471
#ny = 454

#Wrkdir='/project/ofrd/waves/wavegrids/global'
#ModlName='S36125UK'
#nlevs=4
## regular grid info at coarsest level
#dx = 360.0 / 1024.0
#dy = dx * 2.0 / 3.0
#x0 = 0.175781
#y0 = -80.039062
#nx = 1024
#ny = 709
#arcfile = None

#Wrkdir='/project/ofrd/waves/wavegrids/global/GS512L3'
#ModlName='GS512L3'
#nlevs=3
## regular grid info at coarsest level
#dx = 360.0 / 1024.0
#dy = dx * 2.0 / 3.0
#x0 = 0.17578125
#y0 = -80.03906250
#nx = 1024
#ny = 709
#arcfile = None

Wrkdir='/project/ofrd/waves/wavegrids/global/GS512L4EUK'
ModlName='GS512L4EUK'
nlevs=4
# regular grid info at coarsest level
dx = 360.0 / 1024.0
dy = dx * 2.0 / 3.0
x0 = 0.17578125
y0 = -80.03906250
nx = 1024
ny = 709
arcfile = None

## main program
#Wrkdir='/project/ofrd/waves/wavegrids/atlantic/AS512L4EUK'
#ModlName='AS512L4EUK'
#nlevs=4
## regular grid info at coarsest level
#dx = 360.0 / 1024.0
#dy = dx * 2.0 / 3.0
#x0 = -98.26171875
#y0 = -24.25781250
#nx = 471
#ny = 454
#arcfile = None


pltSMCgrid(ModlName, nlevs, dx, dy, x0, y0, nx, ny, Wrkdir, Arc_file=arcfile, figtype='.ps')
