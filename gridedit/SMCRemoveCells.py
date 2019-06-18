#!/usr/bin env python

import numpy as np

def readEditKML(kmlfile):
    '''Read in a kml file containing grid cells to be deleted'''

    with open(kmlfile,'r') as inp:
        txtdata = inp.readlines()
        inp.close()

    xyvals = []
    for lp in range(len(txtdata)):
        chktxt = txtdata[lp].split(':')
        if len(chktxt)>1 and chktxt[0].strip()=='X-Y cell number':        
            xystrs = chktxt[1].split(',') 
            xyvals.append([np.int(xystrs[0]),np.int(xystrs[1])])

    return xyvals


def editSMCTxt(smccellsfile, xyvals, obstrfile=None):
    '''Read SMC cells and obstruction file, then remove points
       at locations read from the Edit file'''

    cellsedit = smccellsfile + '.edit'
    obstredit = smcobstrfile + '.edit'

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

    tmptxt = txtdata[0].split()
    gridnums = []
    for lp in range(len(tmptxt)):
        gridnums.append(np.int(tmptxt[lp]))
    cellsvals  = []
    obstrvals = []
    for lp in range(1,len(txtdata)):
        tmptxt = txtdata[lp].split()
        chkxy = [np.int(tmptxt[0]),np.int(tmptxt[1])]
        if chkxy not in xyvals:
            cellsvals.append(txtdata[lp])
            if obstrfile is not None:
                obstrvals.append(obsdata[lp])
        else:
            gridnums[0] = gridnums[0] - 1
            sclind = np.int(np.log2(np.int(tmptxt[3]))) + 1
            gridnums[sclind] = gridnums[sclind] - 1

    print('Writing to '+cellsedit)
    with open(cellsedit,'w') as outp:
        for lp in range(len(gridnums)):
            outp.write(' %i' %gridnums[lp])
        outp.write('\n')
        for lp in range(len(cellsvals)):
            outp.write('%s' %cellsvals[lp])
        outp.close()

    if obstrfile is not None:
        print('Writing to '+obstredit)
        with open(obstredit,'w') as outp:
            outp.write(' %i' %gridnums[0] + ' 1\n')
            for lp in range(len(obstrvals)):
                outp.write('%s' %obstrvals[lp])
            outp.close()

    return


## main program
#editfile = '/project/ofrd/waves/wavegrids/global/GS512L4EUK/EditGS512L4EUK.kml'
#smccellsfile = '/project/ofrd/waves/wavegrids/global/GS512L4EUK/ww3Cels.dat'
#smcobstrfile = '/project/ofrd/waves/wavegrids/global/GS512L4EUK/ww3Obstr.dat'

editfile = '/project/ofrd/waves/wavegrids/atlantic/AS512L4EUK/EditAS512L4EUK.kml'
smccellsfile = '/project/ofrd/waves/wavegrids/atlantic/AS512L4EUK/ww3Cels.dat'
smcobstrfile = '/project/ofrd/waves/wavegrids/atlantic/AS512L4EUK/ww3Obstr.dat'

xyvals = readEditKML(editfile)
print(xyvals)

editSMCTxt(smccellsfile, xyvals, obstrfile=smcobstrfile)





