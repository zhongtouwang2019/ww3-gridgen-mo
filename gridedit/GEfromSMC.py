#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import iris
from subprocess import call

# write the google earth file header info
def WriteGEHeader( outp, collection_name ):

    outp.write('<?xml version="1.0" encoding="UTF-8"?>\n')
    outp.write('<kml xmlns="http://www.opengis.net/kml/2.2">\n')
    outp.write('<Document>\n')
    outp.write('  <name>%s' %collection_name +'</name>\n')
    outp.write('  <open>0</open>\n')
    outp.write('  <Style id="PolyStyle12">\n')
    outp.write('    <IconStyle>\n')
    outp.write('      <color>ffff9933</color>\n')
    outp.write('      <colorMode>normal</colorMode>\n')
    outp.write('      <scale>1.0</scale>\n')
    outp.write('      <icon>\n')
    outp.write('        <href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>\n')
    outp.write('      </icon>\n')
    outp.write('    </IconStyle>\n')
    outp.write('  </Style>\n')
    outp.write('  <Style id="PolyStyle23">\n')
    outp.write('    <IconStyle>\n')
    outp.write('      <color>ffffff99</color>\n')
    outp.write('      <colorMode>normal</colorMode>\n')
    outp.write('      <scale>1.0</scale>\n')
    outp.write('      <icon>\n')
    outp.write('        <href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>\n')
    outp.write('      </icon>\n')
    outp.write('    </IconStyle>\n')
    outp.write('  </Style>\n')
    outp.write('  <Style id="PolyStyle34">\n')
    outp.write('    <IconStyle>\n')
    outp.write('      <color>ff669933</color>\n')
    outp.write('      <colorMode>normal</colorMode>\n')
    outp.write('      <scale>1.0</scale>\n')
    outp.write('      <icon>\n')
    outp.write('        <href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>\n')
    outp.write('      </icon>\n')
    outp.write('    </IconStyle>\n')
    outp.write('  </Style>\n')
    outp.write('  <Style id="PolyStyle45">\n')
    outp.write('    <IconStyle>\n')
    outp.write('      <color>ff00ff00</color>\n')
    outp.write('      <colorMode>normal</colorMode>\n')
    outp.write('      <scale>1.0</scale>\n')
    outp.write('      <icon>\n')
    outp.write('        <href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>\n')
    outp.write('      </icon>\n')
    outp.write('    </IconStyle>\n')
    outp.write('  </Style>\n')
    outp.write('  <Style id="PolyStyle56">\n')
    outp.write('    <IconStyle>\n')
    outp.write('      <color>ff99ffcc</color>\n')
    outp.write('      <colorMode>normal</colorMode>\n')
    outp.write('      <scale>1.0</scale>\n')
    outp.write('      <icon>\n')
    outp.write('        <href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>\n')
    outp.write('      </icon>\n')
    outp.write('    </IconStyle>\n')
    outp.write('  </Style>\n')
    outp.write('  <Style id="PolyStyle67">\n')
    outp.write('    <IconStyle>\n')
    outp.write('      <color>ff66ffff</color>\n')
    outp.write('      <colorMode>normal</colorMode>\n')
    outp.write('      <scale>1.0</scale>\n')
    outp.write('      <icon>\n')
    outp.write('        <href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>\n')
    outp.write('      </icon>\n')
    outp.write('    </IconStyle>\n')
    outp.write('  </Style>\n')
    outp.write('  <Style id="PolyStyle78">\n')
    outp.write('    <IconStyle>\n')
    outp.write('      <color>ff3399ff</color>\n')
    outp.write('      <colorMode>normal</colorMode>\n')
    outp.write('      <scale>1.0</scale>\n')
    outp.write('      <icon>\n')
    outp.write('        <href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>\n')
    outp.write('      </icon>\n')
    outp.write('    </IconStyle>\n')
    outp.write('  </Style>\n')
    outp.write('  <Style id="PolyStyle89">\n')
    outp.write('    <IconStyle>\n')
    outp.write('      <color>ff6666cc</color>\n')
    outp.write('      <colorMode>normal</colorMode>\n')
    outp.write('      <scale>1.0</scale>\n')
    outp.write('      <icon>\n')
    outp.write('        <href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>\n')
    outp.write('      </icon>\n')
    outp.write('    </IconStyle>\n')
    outp.write('  </Style>\n')
    outp.write('  <Style id="PolyStyle910">\n')
    outp.write('    <IconStyle>\n')
    outp.write('      <color>ff000099</color>\n')
    outp.write('      <colorMode>normal</colorMode>\n')
    outp.write('      <scale>1.0</scale>\n')
    outp.write('      <icon>\n')
    outp.write('        <href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>\n')
    outp.write('      </icon>\n')
    outp.write('    </IconStyle>\n')
    outp.write('  </Style>\n')
    outp.write('  <Style id="PolyStyle10p">\n')
    outp.write('    <IconStyle>\n')
    outp.write('      <color>ff0000cc</color>\n')
    outp.write('      <colorMode>normal</colorMode>\n')
    outp.write('      <scale>1.0</scale>\n')
    outp.write('      <icon>\n')
    outp.write('        <href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>\n')
    outp.write('      </icon>\n')
    outp.write('    </IconStyle>\n')
    outp.write('  </Style>\n')
    outp.write('  <Style id="PolyStyleMinVl">\n')
    outp.write('    <IconStyle>\n')
    outp.write('      <color>ff5c1237</color>\n')
    outp.write('      <colorMode>normal</colorMode>\n')
    outp.write('      <scale>1.0</scale>\n')
    outp.write('      <icon>\n')
    outp.write('        <href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>\n')
    outp.write('      </icon>\n')
    outp.write('    </IconStyle>\n')
    outp.write('  </Style>\n')
    outp.write('  <Style id="PolyStyleBC">\n')
    outp.write('    <IconStyle>\n')
    outp.write('      <color>ff000000</color>\n')
    outp.write('      <colorMode>normal</colorMode>\n')
    outp.write('      <scale>1.0</scale>\n')
    outp.write('      <icon>\n')
    outp.write('        <href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>\n')
    outp.write('      </icon>\n')
    outp.write('    </IconStyle>\n')
    outp.write('    <ListStyle>\n')
    outp.write('      <listItemType>checkHideChildren</listItemType>\n')
    outp.write('    </ListStyle>\n')
    outp.write('  </Style>\n')

    return    

# write the google earth closing statement
def WriteGEEnd( outp ):

    outp.write('</Document>\n')
    outp.write('</kml>\n')

    return

# write the google earth placemark statement
def WriteGEPlace( modelname,lppt,x,y,lat,lon,depth,ylevel,obstrpc=None,minval=10.0 ):

    markerheight = 30.0

    outp.write('  <Placemark>\n')

    # select style for pushpin
    if depth >= 200.0:
        outp.write('      <styleUrl>#PolyStyle12</styleUrl>\n')
    elif depth >= 150.0:
        outp.write('      <styleUrl>#PolyStyle23</styleUrl>\n')
    elif depth >= 120.0:
        outp.write('      <styleUrl>#PolyStyle34</styleUrl>\n')
    elif depth >= 100.0:
        outp.write('      <styleUrl>#PolyStyle45</styleUrl>\n')
    elif depth >= 80.0:
        outp.write('      <styleUrl>#PolyStyle56</styleUrl>\n')
    elif depth >= 60.0:
        outp.write('      <styleUrl>#PolyStyle67</styleUrl>\n')
    elif depth >= 40.0:
        outp.write('      <styleUrl>#PolyStyle78</styleUrl>\n')
    elif depth >= 20.0:
        outp.write('      <styleUrl>#PolyStyle89</styleUrl>\n')
    elif depth >= 10.0:
        outp.write('      <styleUrl>#PolyStyle910</styleUrl>\n')
    elif minval != 10.0:
        if depth >= minval:
            outp.write('      <styleUrl>#PolyStyle10p</styleUrl>\n')
        else:
            outp.write('      <styleUrl>#PolyStyleMinVl</styleUrl>\n')
    else:
      outp.write('      <styleUrl>#PolyStyle10p</styleUrl>\n')

    outp.write('    <description>\n')
    outp.write('            %s' %modelname + '\n')
    outp.write('            Cell index number %d' %lppt + ' (first index=1)\n')
    outp.write('            X-Y cell number: %d' %x + ', %d' %y + '\n')
    outp.write('            Ylevel: %d' %ylevel + ' m\n')
    outp.write('            Position (lat,lon): %8.3f' %lat + ', %8.3f' %lon + ' deg.dec\n')
    if depth < minval:
        outp.write('            Depth: %8.2f' %minval + ' m\n')
    else:
        outp.write('            Depth: %8.2f' %depth + ' m\n')
    if obstrpc is not None:
        outp.write('            Obstr: %4.1f' %obstrpc + ' percent\n')
    outp.write('    </description>\n')
    outp.write('    <Point>\n')
    outp.write('      <coordinates>\n')
    outp.write('            %8.3f' %lon + ', %8.3f' %lat + ',%4.1f' %markerheight +'\n')
    outp.write('      </coordinates>\n')
    outp.write('    </Point>\n')
    outp.write('  </Placemark>\n')

    return

def ReadSMCTxt(smccellsfile,nrlv,jshift,dlon,dlat,slon,slat,rotated=False,rlon=0.0,rlat=0.0,obstrfile=None):

    levscl = 2.0 ** (nrlv-1.0)

    print('Reading '+smccellsfile)
    with open(smccellsfile,'r') as inp:
        txtdata = inp.readlines()
        inp.close()
    obstr = False
    if obstrfile is not None:
        obstr = True
        print('Reading '+obstrfile)
        with open(obstrfile,'r') as inp:
            obsdata = inp.readlines()
            inp.close()

    xy = np.empty( [len(txtdata)-1,2] )
    latlon = np.empty( [len(txtdata)-1,2] )
    depth = np.empty( len(txtdata)-1 )
    xlevel = np.empty( len(txtdata)-1 )
    ylevel = np.empty( len(txtdata)-1 )
    if obstr:
        obstrpc = np.empty( len(txtdata)-1 )

    for lp in range(1,len(txtdata)):
        xy[lp-1,0] = np.float(txtdata[lp].split()[0])
        xy[lp-1,1] = np.float(txtdata[lp].split()[1])
        xlevel[lp-1] = np.float(txtdata[lp].split()[2])
        ylevel[lp-1] = np.float(txtdata[lp].split()[3])
        depth[lp-1] = np.float(txtdata[lp].split()[4])
        if obstr:
           obstrpc[lp-1] = np.float(obsdata[lp])

        # lat lon calculations are based on the 0,0 cell being at slat, slon
        # dlon and dlat are assumed to be prescribed at the coarsest cell resolution
        #latlon[lp-1,0] = (slon - 0.5 * dlon / levscl) + xy[lp-1,0] * dlon / levscl # cell lower left corner
        latlon[lp-1,0] = (slon - 0.5 * dlon) + xy[lp-1,0] * dlon / levscl # cell lower left corner
        latlon[lp-1,0] = latlon[lp-1,0] + (xlevel[lp-1]/2.0) * dlon / levscl # cell centre
        if latlon[lp-1,0] > 180.0:
            latlon[lp-1,0] = latlon[lp-1,0] - 360.0 # put on -180 to 180 grid
        #latlon[lp-1,1] = (slat - 0.5 * dlon / levscl) + xy[lp-1,1] * dlat / levscl # cell lower left corner
        latlon[lp-1,1] = (slat - 0.5 * dlat) + xy[lp-1,1] * dlat / levscl # cell lower left corner
        latlon[lp-1,1] = latlon[lp-1,1] + (ylevel[lp-1]/2.0) * dlat / levscl # cell centre

    if rotated:
        lontmp = latlon[:,0]
        lattmp = latlon[:,1]
        rlonstmp, rlatstmp = iris.analysis.cartography.unrotate_pole(lontmp,lattmp,rlon,rlat)
        latlon[:,0] = rlonstmp
        latlon[:,1] = rlatstmp

    print(np.min(latlon[:,0]))
    print(np.max(latlon[:,0]))
    print(np.min(latlon[:,1]))
    print(np.max(latlon[:,1]))

    print('Read of '+smccellsfile+' completed')

    if obstr:
        return xy, latlon, depth, ylevel, obstrpc
    else:
        return xy, latlon, depth, ylevel, None

#-- main program

# set model
#modelname = 'S36125'
#modelname = 'A36125'
#modelname = 'AMM15SMC'
#modelname = 'GS512L4EUK'
modelname = 'AS512L4EUK'

# set levels constraints for points to output
#levslist = [1,2]
levslist = None

# set region
#UK - standard grid domain
#latlims = [46.0, 62.75]
#lonlims = [-16.0, 13.0]
#regname = 'UK'

# set region
#Euro - standard grid domain
#latlims = [30.2, 65.9]
#lonlims = [-19.8, 41.8]
#regname = 'Euro'

## set region
#latlims = [62, 68]
#lonlims = [-42, -20]
#regname = 'GreenIceChk'

# set region
latlims = [40, 50]
lonlims = [-30, -20]
regname = 'MidAtlChk'

#UK coastal domain
#latlims = [48.5, 62.00]
#lonlims = [-8.0, 3.0]
#regname = 'UKcoastal'

#Norway
#latlims = [57.0, 63.0]
#lonlims = [0.0, 9.0]
#regname = 'Norway'

#Baltic settings
#latlims = [53.0, 66.0]
#lonlims = [9.0, 30.0]
#regname = 'Baltic'

# Carib settings
#latlims = [7.0, 26.0]
#lonlims = [-89.0, -58.0]
#regname = 'Carib'

# KoJap settings
#latlims = [30.0, 41.0]
#lonlims = [121.0, 137.0]
#regname = 'KoJap'

# Med settings
#latlims = [30.0, 47.0]
#lonlims = [-10.0, 42.0]
#regname = 'Med'

# Caspian settings
#latlims = [36.5, 47.0]
#lonlims = [46.0, 54.5]
#regname = 'Caspian'

# Arabia
#latlims = [5.0, 30.0]
#lonlims = [31.0, 78.0]
#regname = 'Arabia'

# Brazil
#latlims = [-30.0, -20.0]
#lonlims = [-55.0, -35.0]
#regname = 'Brazil'

# Ascension
#latlims = [-9.0, -7.0]
#lonlims = [-15.0, -13.0]
#regname = 'Ascension'

#UK SW settings
#latlims = [49.6, 51.8]
#lonlims = [-6.6, 0.5]
#regname = 'UK-SW'

#UK NE settings
#latlims = [54.3, 55.4]
#lonlims = [-1.65, -0.25]
#regname = 'UK-NE'


# set the input file data
if modelname == 'S36125':
    smccellsfile = '/hpc/home/d01/frxs/FCM_WW3CONFIG/trunk/GblSMC_361225/grid/S36125MCels.dat'
    nrlv = 4
    jshift = 2816
    dlon = 0.35156250
    dlat = 0.23437500
    slon = 0.0
    slat = 0.0
    rotated = False
    rlon = 0.0
    rlat = 0.0
    depthmin = 15.0
    if levslist == None:
        levslist = [1,2,4,8]
elif modelname == 'A36125':
    smccellsfile = '/hpc/home/d01/frxs/FCM_WW3CONFIG/trunk/AtlSMC_361225/grid/Atlan36125.dat'
    nrlv = 4
    jshift = 2816
    dlon = 0.35156250
    dlat = 0.23437500
    slon = 0.0
    slat = 0.0
    rotated = False
    rlon = 0.0
    rlat = 0.0
    depthmin = 15.0
    if levslist == None:
        levslist = [1,2,4,8]
elif modelname == 'AMM15SMC':
    smccellsfile = '/hpc/home/d01/frxs/FCM_WW3CONFIG/trunk/amm15_smc/grid/amm15s.ww3Cels.dat'
    nrlv = 2
    jshift = 2816
    dlon = 0.0270 # value for coarsest cell size
    dlat = 0.0270 # value for coarsest cell size
    slon = -10.8895 # cell centre
    slat = -7.2942 # cell centre
    rotated = True
    rlon = 177.5
    rlat = 37.5
    depthmin = 10.0
    if levslist == None:
        levslist = [1,2]
elif modelname == 'GS512L4EUK':
    smccellsfile = '/hpc/home/d01/frxs/FCM_WW3CONFIG/r2214_PS43/GW1.0/configs/'+modelname+'/ww3Cels.dat'
    smcobstrfile = '/hpc/home/d01/frxs/FCM_WW3CONFIG/r2214_PS43/GW1.0/configs/'+modelname+'/ww3Obstr.dat'
    nrlv = 4
    jshift = 0
    dlon = 0.35156250
    dlat = 0.23437500
    slon = 0.17578125
    slat = -80.03906250
    rotated = False
    rlon = 0.0
    rlat = 0.0
    depthmin = 15.0
    if levslist == None:
        levslist = [1,2,4,8]
elif modelname == 'AS512L4EUK':
    smccellsfile = '/hpc/home/d01/frxs/FCM_WW3CONFIG/r2214_PS43/GW1.0/configs/'+modelname+'/ww3Cels.dat'
    smcobstrfile = '/hpc/home/d01/frxs/FCM_WW3CONFIG/r2214_PS43/GW1.0/configs/'+modelname+'/ww3Obstr.dat'
    nrlv = 4
    jshift = 0
    dlon = 0.35156250
    dlat = 0.23437500
    slon = -98.26171875 
    slat = -24.25781250
    rotated = False
    rlon = 0.0
    rlat = 0.0
    depthmin = 15.0
    if levslist == None:
        levslist = [1,2,4,8]


# set working file names and variables
workdir = '/data/cr1/frxs/WW3Grids'
GEBathyFile = workdir + '/' + modelname+'_'+regname+'.kml'

# read in the SMC text file data
xy, latlon, depth, ylevel, obstrpc = ReadSMCTxt(smccellsfile,nrlv,jshift,dlon,dlat,
                                                 slon,slat,rotated,rlon,rlat,
                                                 obstrfile=smcobstrfile)

# write out to ge file
with open(GEBathyFile,'w') as outp:

    WriteGEHeader( outp, modelname+': '+regname )

    for lppt in range(np.shape(xy)[0]):
        # level constraint
        if ylevel[lppt] in levslist:
            # lat constraint
            if latlon[lppt,1] >= latlims[0] and latlon[lppt,1] <= latlims[1]:
                # lon constraint
                if latlon[lppt,0] >= lonlims[0] and latlon[lppt,0] <= lonlims[1]:
                    if obstrpc is not None:
                       WriteGEPlace(modelname,lppt+1,xy[lppt,0],xy[lppt,1],latlon[lppt,1],latlon[lppt,0],
                                     depth[lppt],ylevel[lppt],obstrpc=obstrpc[lppt],minval=depthmin)
                    else:
                        WriteGEPlace(modelname,lppt+1,xy[lppt,0],xy[lppt,1],latlon[lppt,1],latlon[lppt,0],
                                      depth[lppt],ylevel[lppt],minval=depthmin)
        if np.mod(lppt,100000) == 0:
            print('Searched %d' %lppt + ' points out of %d' %np.shape(xy)[0])


    WriteGEEnd( outp )

    outp.close()


# zip and copy to kmz file
call(["zip","-r",workdir+'/'+modelname+"_"+regname+".kmz",workdir+'/'+modelname+"_"+regname+".kml"])
#call(["mv",workdir+modelname+".zip",workdir+modelname+".kmz"])
call(["cp",workdir+'/'+modelname+"_"+regname+".kmz","./"])


