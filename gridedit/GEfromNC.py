#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from subprocess import call

# write the google earth file header info
def WriteGEHeader( outp ):

    outp.write('<?xml version="1.0" encoding="UTF-8"?>\n')
    outp.write('<kml xmlns="http://www.opengis.net/kml/2.2">\n')
    outp.write('<Document>\n')
    outp.write('  <name>name</name>\n')
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
def WriteGEPlace( modelname,x,y,lat,lon,depth,minval=10.0 ):

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
    outp.write('            X-Y cell number: %d' %x + ', %d' %y + '\n')
    outp.write('            Position (lat,lon): %8.3f' %lat + ', %8.3f' %lon + ' deg.dec\n')
    outp.write('            Depth: %8.2f' %depth + ' m\n')
    outp.write('    </description>\n')
    outp.write('    <Point>\n')
    outp.write('      <coordinates>\n')
    outp.write('            %8.3f' %lon + ', %8.3f' %lat + ',%4.1f' %markerheight +'\n')
    outp.write('      </coordinates>\n')
    outp.write('    </Point>\n')
    outp.write('  </Placemark>\n')

    return



# set the input file data
workdir = '/data/cr1/frxs/WW3Grids/'

#WW3file = '/hpc/data/d01/frxs/AS20W/outputs/as20w_20160101.nc'
#modelname = 'AS20W'

WW3file = '/hpc/data/d01/frxs/UAE30S/outputs/uae30s_20160101.nc'
modelname = 'UAE30S'

#WW3file = '/hpc/projects/ocean/nemo/AMM7/CO6_INPUT_DATA/Nemo_Inputs/CO6_MOMME_Corrected_bathy_meter.nc'
#modelname = 'AMM7_CO6'

# set a minimum depth value used by the model
depthmin = 5.0

# set working file names and variables
GEBathyFile = workdir + modelname+'.kml'

# read in the NEMO file data
print('Reading depth data from ' + WW3file)
d = nc.Dataset(WW3file)

depths = d.variables['dpt']
lats   = d.variables['latitude']
lons   = d.variables['longitude']

#depths = d.variables['Bathymetry']
#lats   = d.variables['lat']
#lons   = d.variables['lon']


# write out to ge file
with open(GEBathyFile,'w') as outp:

    WriteGEHeader( outp )

    for x in range(len(lons)):
        for y in range(len(lats)):
            if depths[0,y,x] > 0.0:
                WriteGEPlace(modelname,x,y,lats[y],lons[x],depths[0,y,x],minval=depthmin)
            if y==len(lats)-1 and np.mod(x,100) == 0:
                print 'Written data from %d' %x + ' columns out of %d' %len(lons)


    WriteGEEnd( outp )

    outp.close()


# zip and copy to kmz file
call(["zip","-r",workdir+modelname+".kmz",workdir+modelname+".kml"])
#call(["mv",workdir+modelname+".zip",workdir+modelname+".kmz"])
call(["cp",workdir+modelname+".kmz","./"])


