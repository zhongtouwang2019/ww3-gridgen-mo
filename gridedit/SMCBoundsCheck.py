#!/usr/bin env python

import numpy as np

def createBdycln(smccellsfile, bpfile, outdir='.'):
    '''Generate a boundary point index file from cell and boundary point lists'''

    bpindfile = outdir + '/ww3Bdycln.dat'

    print('Reading '+smccellsfile)
    with open(smccellsfile,'r') as inp:
        smcdata = inp.readlines()
        inp.close()

    print('Reading '+bpfile)
    with open(bpfile,'r') as inp:
        bpdata = inp.readlines()
        inp.close()
    bplist = []
    for lp in range(len(bpdata)):
        bpsplit = bpdata[lp].split()
        bplist.append([np.int(bpsplit[0]), np.int(bpsplit[1]),
                        np.int(bpsplit[2]), np.int(bpsplit[3])])

    print('Checking for matching points and writing boundary cell indices')
    print('Writing results out to %s' %bpindfile)
    with open(bpindfile, 'w') as outp:
        for lp in range(1,len(smcdata)):
            clsplit = smcdata[lp].split()
            chkcell = [np.int(clsplit[0]), np.int(clsplit[1]),
                        np.int(clsplit[2]), np.int(clsplit[3])]
            if chkcell in bplist:
                print('Found boundary point at index %d' %lp)
                outp.write(' %d\n' %lp)
        outp.close()

    return

## main program
smccellsfile = '/project/ofrd/waves/wavegrids/atlantic/AS512L4EUK/ww3Cels.dat'
bpfile = '/project/ofrd/waves/wavegrids/atlantic/AS512L4EUK/ww3BoundCels.dat'
outdir = '/project/ofrd/waves/wavegrids/atlantic/AS512L4EUK'
createBdycln(smccellsfile, bpfile, outdir=outdir)
