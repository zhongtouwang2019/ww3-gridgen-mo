# library for gridgen functions
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import copy

class regGrid:

    def __init__(self, name='regular_grid'):
        self.name = name
        self.rtd = False     # rotated grid logical
        self.plat = 90.0
        self.plon = 0.0
        self.depthmin = 10.0 # minimum depth for model grid
        self.dlim = 0.0      # dry depth threshold
        self.drymin = 0.01   # proportion cells to define wet point
        self.drymax = 0.99   # proportion cells to define dry point
        self.llx = None
        self.lly = None
        self.dx = None
        self.dy = None
        self.nx = None
        self.ny = None
        self.midcellsx = None
        self.midcellsy = None
        self.cellbounds = None
        self.celldepths = None

    def setDryLimits(self, mindepth=10.0, dlim=0.0, drymin=0.0, drymax=0.99):
        '''Update the limits for defining model minimum depth and wet and dry cells'''
        self.mindepth = mindepth
        self.dlim   = dlim
        self.drymin = drymin
        self.drymax = drymax

    def setExtents(self, llx, lly, urx, ury, dx, dy):
        '''Sets the grid; for regular defines llx,lly as cell midpoint'''
        nx, ny = setXYdims(llx, lly, urx, ury, dx, dy)
        self.llx = llx
        self.lly = lly
        self.dx = dx
        self.dy = dy
        self.nx = nx
        self.ny = ny

    def setGrid(self):
        '''Sets the grid using cell midpoints'''
        self.midcellsx, self.midcellsy = setGridRegular(self.llx, self.lly,
                                                        self.dx, self.dy,
                                                        self.nx, self.ny)
    def setBounds(self):
        '''Set cell boundaries for fill depth calculations'''
        self.cellbounds = setCellBoundsRegular(self.midcellsx, self.midcellsy,
                                               self.dx, self.dy)

    def setDepths(self, rdlats, rdlons, rdbathy, median_depth=False):
        '''Set cell depths using fill depth calculations'''
        self.celldepths = fillCells(self.cellbounds, rdlats, rdlons,
                                    rdbathy, dlim=self.dlim,
                                    drymin=self.drymin, drymax=self.drymax,
                                    rotated=self.rtd, median_depth=False)

    def writeWW3(self, writedir='.'):
        '''Write the grid data out to ww3_grid compatible files'''
        writedepths = cells2grid(self.midcellsx, self.midcellsy,
                                 self.celldepths[:,0])
        writemask = cells2grid(self.midcellsx, self.midcellsy,
                                self.celldepths[:,2])
        writeblock = cells2grid(self.midcellsx, self.midcellsy,
                                self.celldepths[:,1])
        writeWW3regular(writedepths, writemask, writeblock, 
                        self.llx, self.lly, self.dx, self.dy, self.nx, self.ny, 
                        llscale=1.0, ldscale=1.0, 
                        depthlim=self.dlim, mindepth=self.mindepth, 
                        depscale=1.0, blkscale=100.0, 
                        rtd=self.rtd, plat=self.plat, plon=self.plon,
                        writedir='.')        

class smcGrid:

    def __init__(self, name='smc_grid', label='basegrid'):
        self.name = name
        self.label = label
        self.ntiers = 1
        self.rtd = False     # rotated grid logical
        self.plat = 90.0
        self.plon = 0.0
        self.depthmin = 10.0 # model minimum depth
        self.dlim = 0.0      # dry depth threshold
        self.drymin = 0.01   # proportion dry cells to define wet cell
        self.drymax = 0.99   # proportion dry cells to define dry cell
        self.llx = None # lower left hand corner x
        self.lly = None # lower left hand corner y
        self.dx = None  # dx for highest tier resolution
        self.dy = None  # dy for highest tier resolution
        self.nx = None  # number x cells for base resolution (regular grid)
        self.ny = None  # number y cells for base resolution (regular grid)
        self.midcellsx = None
        self.midcellsy = None
        self.smccells  = None
        self.cellbounds = None
        self.celldepths = None

    def setDryLimits(self, mindepth=None, dlim=None, drymin=None, drymax=None):
        '''Update the limits for defining model minimum depth and wet and dry cells'''
        if mindepth is not None:
            self.mindepth = mindepth
        if dlim is not None:
            self.dlim   = dlim
        if drymin is not None:
            self.drymin = drymin
        if drymax is not None:
            self.drymax = drymax

    def setExtents(self, llx, lly, urx, ury, dx, dy):
        '''Sets the grid; for SMC defines llx,lly as cell lower left corner'''
        nx, ny = setXYdims(llx, lly, urx, ury, dx, dy)
        self.llx = llx
        self.lly = lly
        self.dx = dx
        self.dy = dy
        self.nx = nx
        self.ny = ny

    def setGridFromRegular(self):
        '''Sets the grid using cell lower left corners'''
        self.midcellsx, self.midcellsy, self.smccells = setGridRegular(
                                                           self.llx, self.lly,
                                                           self.dx, self.dy,
                                                           self.nx, self.ny,
                                                           midpoint=False,
                                                           smc=True)

    def setCFLCells(self):
        '''Applies CFL relaxation criteria to high latitude cells'''
        print('[INFO] Applying CFL relaxation for high latitude cells')
        clatmin = np.cos(np.min(np.abs(self.midcellsy)) * np.pi / 180.0)
        xmax = np.max(self.smccells[:,0])
        rtier = 1
        # rtier <=10 should allow us to get to 89.95N from equatorial cell
        while rtier <= 10:
            # test for cells at factor 2 reduction from largest dx
            clats = np.cos(self.midcellsy * np.pi / 180.0)
            inds = np.where(clats <= clatmin/2.0**rtier)[0]
            if np.size(inds) == 0:
                rtier = 99
            else:
                # double cell sizes in x dimension and mark every other
                # cell for removal
                print('[INFO] Identified %d cells for CFL relaxation step %d'
                      %tuple([np.size(inds),rtier]))
                evenx = np.where(np.mod(self.smccells[inds,0],2**rtier)==0)
                self.smccells[inds[evenx],2] = self.smccells[inds[evenx],2] * 2
                self.midcellsx[inds[evenx]] = self.midcellsx[inds[evenx]] + \
                                              self.dx * 2.0**(rtier-2)
                oddx = np.where(np.mod(self.smccells[inds,0],2**rtier)!=0)
                self.smccells[inds[oddx],2] = -1
                # if last cell in row is even, mark for removal
                # this prevents the grid overlapping the preset extents
                #if np.mod(xmax,2**rtier) == 0:
                maxx = np.where((self.smccells[inds,0]+self.smccells[inds,2]/2) > xmax)
                self.smccells[inds[maxx],2] = -1
                # remove unused cells
                retain = np.where(self.smccells[:,2] != -1)[0]
                self.smccells = self.smccells[retain,:]
                self.midcellsx = self.midcellsx[retain]
                self.midcellsy = self.midcellsy[retain]
            rtier = rtier + 1

    def setBounds(self):
        '''Set cell boundaries for fill depth calculations'''
        self.cellbounds = setCellBoundsSMC(self.smccells, self.llx,
                                           self.lly, self.dx, self.dy)

    def setMids(self):
        '''Set cell centres'''
        self.midcellsx, self.midcellsy = setCellMidsSMC(self.smccells, 
                                                        self.llx, self.lly,
                                                        self.dx, self.dy)

    def setDepths(self, rdlats, rdlons, rdbathy, pland=None, median_depth=False, setadj=True):
        '''Set cell depths using fill depth calculations'''
        self.celldepths = fillCells(self.cellbounds, rdlats, rdlons,
                                    rdbathy, dlim=self.dlim,
                                    drymin=self.drymin, drymax=self.drymax,
                                    pland=pland, rotated=self.rtd, 
                                    median_depth=median_depth, smc=True, setadj=setadj)

    def markDepths(self, depthlim, marker='tier'):
        '''Marks cells shallower than given depth for tier'''
        print('[INFO] Marking cells above depth limit %.2f m' %depthlim)
        if depthlim > 0.0:
            print('[WARN] Depth is set greater than zero, changing sign for depth negative convention')
            depthlim = depthlim * -1.0
        if marker.lower() == 'tier':
            self.celldepths[(self.celldepths[:,0]>=depthlim) & (self.celldepths[:,2]==1), 2] = -1
        elif marker.lower() == 'dry':
            self.celldepths[self.celldepths[:,0]>=depthlim, 2] = 0

    def markDrypc(self, pcdry):
        '''Marks cells with dry percentage greater than threshold as dry'''
        self.celldepths[self.celldepths[:,1]>=pcdry, 2] = 0

    def markRegion(self, xsw, ysw, xne, yne, marker='tier', depthlim=None):
        '''Marks cells in a given region for tier'''
        print('[INFO] Setting cells in region %.2fE, %.2fN, %.2fE, %.2fN'
                %tuple([xsw, ysw, xne, yne]) + ' to ' + marker)
        inds = np.where(((self.midcellsx[:] >= xsw) & 
                          (self.midcellsx[:] < xne) &
                          (self.midcellsy[:] >= ysw) & 
                          (self.midcellsy[:] < yne)))
        if marker.lower() == 'tier':
            if depthlim is not None:
                print('[INFO] Only marking cells above depth limit %.2f m' %depthlim)
                if depthlim > 0.0:
                    print('[WARN] Depth is set greater than zero, changing sign for depth negative convention')
                    depthlim = depthlim * -1.0
                for ind in inds[0]:
                    if self.celldepths[ind,2] == 1 and \
                        self.celldepths[ind,0] > depthlim:
                        self.celldepths[ind,2] = -1
            else:
                for ind in inds[0]:
                    if self.celldepths[ind,2] == 1:
                        self.celldepths[ind,2] = -1
        elif marker.lower() == 'dry':
            self.celldepths[inds, 2] = 0

    def unmarkCells(self, marker='tier', box=None, osbox=True, thruzero=False):
        '''Sets marked cells to wet setting'''
        if box is not None:
            xsw = box[0]
            ysw = box[1]
            xne = box[2]
            yne = box[3]
            if osbox:
                print('[INFO] Setting cells outside region %.2fE, %.2fN, %.2fE, %.2fN'
                       %tuple(box[:]) + ' to wet')
                if thruzero:
                    inds = np.where((((self.midcellsx[:] < xsw) & 
                                       (self.midcellsx[:] > xne)) |
                                       (self.midcellsy[:] < ysw) | 
                                       (self.midcellsy[:] > yne)))
                else:
                    inds = np.where(((self.midcellsx[:] < xsw) | 
                                      (self.midcellsx[:] > xne) |
                                      (self.midcellsy[:] < ysw) | 
                                      (self.midcellsy[:] > yne)))
            else:
                print('[INFO] Setting cells in region %.2fE, %.2fN, %.2fE, %.2fN'
                       %tuple(box[:]) + ' to wet')
                if thruzero:
                    inds = np.where((((self.midcellsx[:] >= xsw) | 
                                       (self.midcellsx[:] < xne)) &
                                       (self.midcellsy[:] >= ysw) & 
                                       (self.midcellsy[:] < yne)))
                else:
                    inds = np.where(((self.midcellsx[:] >= xsw) & 
                                      (self.midcellsx[:] < xne) &
                                      (self.midcellsy[:] >= ysw) & 
                                      (self.midcellsy[:] < yne)))
            if marker.lower() == 'tier':
                for ind in inds[0]:
                    if self.celldepths[ind,2] < 0:
                        self.celldepths[ind,2] = self.celldepths[ind,2]+2
            elif marker.lower() == 'dry':
                for ind in inds[0]:
                    if self.celldepths[ind,2] == 0 or self.celldepths[ind,2] == -2:
                        self.celldepths[ind,2] = 1
        else:
            if marker.lower() == 'tier':
                self.celldepths[self.celldepths[:,2] == -1, 2] = 1
                self.celldepths[self.celldepths[:,2] == -2, 2] = 0
            elif marker.lower() == 'dry':
                self.celldepths[self.celldepths[:,2] == 0, 2] = 1

    def delDryCells(self,celltype='dry'):
        '''Remove dry cells from SMC grid arrays'''
        self.smccells, self.midcellsx, self.midcellsy, \
        self.cellbounds, self.celldepths = removeCells(self.smccells,
                                                       self.midcellsx,
                                                       self.midcellsy,
                                                       self.cellbounds,
                                                       self.celldepths,
                                                       celltype=celltype)

    def delTierCells(self,celltype='tier'):
        '''Remove tier cells from SMC grid arrays'''
        self.smccells, self.midcellsx, self.midcellsy, \
        self.cellbounds, self.celldepths = removeCells(self.smccells,
                                                       self.midcellsx,
                                                       self.midcellsy,
                                                       self.cellbounds,
                                                       self.celldepths,
                                                       celltype=celltype)

    def sortCells(self):
        self.smccells, self.midcellsx, self.midcellsy, \
        self.cellbounds, self.celldepths = sortSMC(self.smccells,
                                                   self.midcellsx,
                                                   self.midcellsy,
                                                   self.cellbounds,
                                                   self.celldepths)

    def writeWW3(self, writedir='.', mindepth=None, writemindepth=False):
        '''Write the grid data out to ww3_grid compatible files'''
        if mindepth is None:
            mindepth = self.depthmin
        writeWW3smc(self.smccells, self.celldepths, self.ntiers,  
                    self.llx, self.lly, self.dx, self.dy, self.nx, self.ny, 
                    llscale=1.0, ldscale=1.0, 
                    depthlim=self.dlim, mindepth=mindepth,
                    writemindepth = writemindepth, 
                    depscale=1.0, blkscale=100.0, 
                    rtd=self.rtd, plat=self.plat, plon=self.plon,
                    writedir=writedir)        

    def writeBounds(self, writedir='.', lon360=True,
                     north=False, east=False, south=False, west=False):
        '''Writes boundary cell centre values to ww3_grid compatible files'''
        writeBoundsmc(self, writedir=writedir, lon360=lon360, 
                       north=north, east=east, south=south, west=west)

    def writeNC(self, writedir='.', filename=None):
        '''Write core SMC grid and metadata out to a netCDF file'''
        writeNCsmc(self, writedir=writedir, filename=filename)


#--- helper functions ---        

def chkAdj(ind, smcbounds, altbounds=None, nexttier=False):
    '''Finds cells intersecting indexed cell corners'''

    # corners are ordered as SW, NW, NE, SE
    # 1st corners for box we want to test, 2nd corner for intersect to test
    if nexttier:
        corners = [ [[0,1],[0,3]], [[0,1],[2,1]], 
                    [[0,3],[0,1]], [[0,3],[2,3]],
                    [[2,3],[0,3]], [[2,3],[2,1]],
                    [[2,1],[2,3]], [[2,1],[0,1]] ]
        intersects = [None, None, None, None, None, None, None, None]
    else:
        corners = [ [[0,1],[2,1]], [[0,3],[0,1]], [[2,3],[0,3]], [[2,1],[2,3]] ]
        intersects = [None, None, None, None]

    if altbounds is not None:
        testbounds = altbounds
    else:
        testbounds = smcbounds

    for lp, crnrs in enumerate(corners):
        xcrnr = smcbounds[ind, crnrs[0][0]]
        ycrnr = smcbounds[ind, crnrs[0][1]]
        # find cell with interescting corner
        distsw = np.sqrt((testbounds[:,crnrs[1][0]]-xcrnr)**2. + \
                                  (testbounds[:,crnrs[1][1]]-ycrnr)**2.)
        mdind = np.argmin(distsw)
        if distsw[mdind] < 0.0001:
            intersects[lp] = mdind

    return intersects


def setXYdims(llx, lly, urx, ury, dx, dy):
    '''Define regular grid nx,ny based on extents and dx,dy'''

    nx = np.int(np.ceil((urx - llx) / dx))
    ny = np.int(np.ceil((ury - lly) / dy))
    print('[INFO] Regular grid dimensions nx,ny = %d,%d' %tuple([nx, ny]))

    return nx, ny


def setGridRegular(llx, lly, dx, dy, nx, ny, midpoint=True,
                   smc=False):
    '''Set up a regular grid mesh based on extents and dx,dy'''

    if smc:
        print('[INFO] Returning grid with %d seapoints' %(nx*ny))
    else:
        print('[INFO] Returning grid with nx,ny = %d,%d' %tuple([nx, ny]))
    if midpoint:
        midoffsx = 0.0
        midoffsy = 0.0
    else:
        midoffsx = dx / 2.0
        midoffsy = dy / 2.0

    midcellsx = llx + np.arange(nx) * dx + midoffsx
    midcellsy = lly + np.arange(ny) * dy + midoffsy
    if smc:
        smccellsx = np.tile(midcellsx,ny)
        smccellsy = np.repeat(midcellsy,nx)
        smccells  = np.empty([nx*ny,4],dtype=int)
        smccells[:,0] = np.tile(np.arange(nx),ny)
        smccells[:,1] = np.repeat(np.arange(ny),nx)
        smccells[:,2] = np.ones(nx*ny)
        smccells[:,3] = np.ones(nx*ny)
        return smccellsx, smccellsy, smccells
    else:
        return midcellsx, midcellsy


def setCellBoundsRegular(midcellsx, midcellsy, dx, dy):
    '''Returns a cells boundary list for a regular grid,
       using mid cell values and dx,dy'''

    print('[INFO] Setting cell boundaries array for regular grid')
    nx = len(midcellsx)
    ny = len(midcellsy)
    cell_bounds = np.zeros([ny*nx, 4])
    for lpy, y in enumerate(midcellsy):
        for lpx, x in enumerate(midcellsx):
            lpt = lpy * nx + lpx
            cell_bounds[lpt,0] = midcellsx[lpx] - 0.5 * dx # sw corner x
            cell_bounds[lpt,1] = midcellsy[lpy] - 0.5 * dy # sw corner y
            cell_bounds[lpt,2] = midcellsx[lpx] + 0.5 * dx # ne corner x
            cell_bounds[lpt,3] = midcellsy[lpy] + 0.5 * dy # ne corner y

    return cell_bounds


def setCellBoundsSMC(smccells, llx, lly, dx, dy):
    '''Returns a cells boundary list for a SMC grid,
       using SMC cells, lower left corners and dx,dy'''

    print('[INFO] Setting cell boundaries array for SMC grid')
    cell_bounds = np.zeros([np.shape(smccells)[0], 4])
    for lp in range(np.shape(smccells)[0]):
        cell_bounds[lp,0] = llx + np.float(smccells[lp,0]) * dx  # sw corner x
        cell_bounds[lp,1] = lly + np.float(smccells[lp,1]) * dy  # sw corner y
        cell_bounds[lp,2] = llx + np.float((smccells[lp,0] + \
                             smccells[lp,2])) * dx               # ne corner x
        cell_bounds[lp,3] = lly + np.float((smccells[lp,1] + \
                             smccells[lp,3])) * dy               # ne corner y

    return cell_bounds


def setCellMidsSMC(smccells, llx, lly, dx, dy):
    '''Returns cells centre lists for a SMC grid,
       using SMC cells, lower left corners and dx,dy'''

    print('[INFO] Setting cell centres array for SMC grid')
    midcellsx = np.zeros(np.shape(smccells)[0])
    midcellsy = np.zeros(np.shape(smccells)[0])
    for lp in range(np.shape(smccells)[0]):
        midcellsx[lp] = llx + (np.float(smccells[lp,0]) + \
                               np.float(smccells[lp,2]) / 2.0) * dx
        midcellsy[lp] = lly + (np.float(smccells[lp,1]) + \
                               np.float(smccells[lp,3]) / 2.0) * dy

    return midcellsx, midcellsy


def fillCells(cell_bounds, rdx, rdy, rdbathy, dlim=0.0, drymin=0.0, 
               drymax=0.99, pland=None, rotated=False, 
               median_depth=False, smc=False, setadj=False):
    '''Returns a list of depth and land-sea data to correspond
    with cell bounds list'''

    print('[INFO] Calculating cell depths')
    ncells = np.shape(cell_bounds)[0]
    # cell depths array as depth, proportion of dry cells and cell type
    cell_depths = np.zeros([ncells,3])
    cell_depths[:,2] = 1 # set to default ww3 wet cell value

    if dlim > 0.0:
        print('[WARN] Dry depth limit is set greater than zero, changing sign for depth negative convention')
        dlim = dlim * -1.0

    # if rdx and rdy are 1D arrays, combine to form 2d arrays
    #if len(np.shape(rdx)) == 1:
    #    chkx, chky = np.meshgrid(rdx, rdy)
    #else:
    chkx = rdx
    chky = rdy

    for lp in range(np.shape(cell_bounds)[0]):
        if np.mod(lp, 2500) == 0:
            print('[INFO] ... done %d points out of %d' %tuple([lp, ncells]))
        xsw = cell_bounds[lp,0]
        ysw = cell_bounds[lp,1]
        xne = cell_bounds[lp,2]
        yne = cell_bounds[lp,3]
        if len(np.shape(rdx)) == 1:
            # regular bathy
            indsx = np.where((chkx >= xsw) & (chkx < xne))
            indsy = np.where((chky >= ysw) & (chky < yne))
            ndepths = np.size(indsx) * np.size(indsy)
        else:
            # rotated pole bathy
            inds = np.where(((chkx >= xsw) & (chkx < xne) &
                              (chky >= ysw) & (chky < yne)))
            ndepths = np.size(inds) / 2
        if ndepths > 0:
            if len(np.shape(rdx)) == 1:
                # regular bathy
                bathytmp = rdbathy[np.min(indsy):np.max(indsy)+1,
                                    np.min(indsx):np.max(indsx)+1].flatten()
            else:
                # rotated pole bathy
                bathytmp = rdbathy[inds]
            # only use wet depths in calculations
            if np.size(bathytmp[bathytmp<dlim]) > 0:
                if median_depth:
                    depth = np.median(bathytmp[bathytmp<dlim])
                else:
                    depth = np.mean(bathytmp[bathytmp<dlim])
            else:
                depth = 99.99
            # use all depths for dry percentage calculation
            pcdry = np.size(np.where(bathytmp >= dlim)[0])
            # add wet cell land percentages if this info has been loaded in
            if pland is not None:
                if len(np.shape(rdx)) == 1:
                    # regular bathy
                    plandtmp = pland[np.min(indsy):np.max(indsy)+1,
                                      np.min(indsx):np.max(indsx)+1].flatten()
                else:
                    # rotated pole bathy
                    plandtmp = pland[inds]
                if np.size(bathytmp[bathytmp < dlim]) > 0:
                    plandsum = np.sum(plandtmp[bathytmp < dlim])
                    pcdry = np.float(pcdry) + plandsum
            pcdry = np.float(pcdry) / np.float(ndepths)
            cell_depths[lp,0] = depth
            cell_depths[lp,1] = pcdry
            # mark cells for removal/tiering based on percentage dry
            if pcdry >= drymax:
                # reject dry cells
                cell_depths[lp,2] = 0
            elif pcdry > drymin:
                # set partially dry points for tiering
                cell_depths[lp,2] = -1
        else:
            print('[WARNING] No source data found in cell, returning zero value')

    # second pass through cells to switch cells adjacent to coast to type -2
    # sets required additional tiering in next step
    if smc and setadj:
        print('[INFO] Checking for points adjacent to dry cells')
        inds = np.where(cell_depths[:,2] == 0)
        adjdry = []
        for cnt, lp in enumerate(inds[0]):
            if np.mod(cnt, 2500) == 0:
                print('[INFO] ... done %d points out of %d' %tuple([cnt, np.size(inds)]))
            intersects = chkAdj(lp, cell_bounds, altbounds=None)
            switch_drytype = False
            if np.any(intersects is not None):
                for chkcell in intersects:
                    if chkcell is not None:
                        if cell_depths[chkcell,2] != 0:
                            cell_depths[chkcell,2] = -1
                            switch_drytype = True
            if switch_drytype:
                adjdry.append(lp)
        if np.size(np.array(adjdry)) > 0:
            cell_depths[adjdry,2] = -2

    # for non-smc grids set cells marked -1 to 1 (wet)
    if not smc:
        print('[INFO] Not SMC grid - switching tier values to wet cells')
        cell_depths[cell_depths[:,2] == -2, 2] = -1
        cell_depths[:,2] = np.abs(cell_depths[:,2])
        print(np.min(cell_depths[lp,2]))

    return cell_depths


def cells2grid(midx, midy, cellsarr):
    '''Map 1D array for regular grid back to 2D xy array'''
    nx = len(midx)
    ny = len(midy)
    cells_grid = cellsarr.reshape([ny, nx])

    return cells_grid


def smcTier(llx, lly, dx, dy, smccells, cellbounds, celldepths=None):
    '''Split cell bound box into 4 cells.
       If celldepths are provided split based on percentage
       of dry points'''

    if celldepths is not None:
        bounds_tmp = cellbounds[np.where(celldepths[:,2] < 0)[0],:]
        smc_tmp    = smccells[np.where(celldepths[:,2] < 0)[0],:]
    else:
        bounds_tmp = np.copy(cellbounds)
        smc_tmp    = np.copy(smccells)
    ncells = np.shape(bounds_tmp)[0]

    ndx = dx / 2.0
    ndy = dy / 2.0

    print('[INFO] Call to SMC Tier will return %d new cells' %(ncells*4))
    new_smc    = np.zeros([ncells*4,4],dtype=int)
    for lp in range(ncells):
        for lpc in range(4):
            new_smc[lp*4+lpc,0] = smc_tmp[lp,0] * 2 + np.int(np.floor(lpc / 2)) * \
                                  smc_tmp[lp,2]
            new_smc[lp*4+lpc,1] = smc_tmp[lp,1] * 2 + np.int(np.mod(lpc, 2)) * \
                                  smc_tmp[lp,3]
            new_smc[lp*4+lpc,2] = smc_tmp[lp,2]
            new_smc[lp*4+lpc,3] = smc_tmp[lp,3]
    new_bounds = setCellBoundsSMC(new_smc, llx, lly, ndx, ndy)

    rept_checks = False
    if rept_checks:
        print('')
        print('*** Checking new_smc for repeats')
        bdchk = []
        for lp in range(np.shape(new_smc)[0]):
            if tuple(new_smc[lp,:]) not in bdchk:
                bdchk.append(tuple(new_smc[lp,:]))
            else:
                print('Repeat: ',lp,new_smc[lp,:])
            if np.mod(lp,1000) == 0:
                print('Checked %d' %lp)

        print('')
        print('*** Checking new_bounds for repeats')
        bdchk = []
        for lp in range(np.shape(new_bounds)[0]):
            if tuple(new_bounds[lp,:]) not in bdchk:
                bdchk.append(tuple(new_bounds[lp,:]))
            else:
                print('Repeat: ',lp,new_bounds[lp,:])
            if np.mod(lp,1000) == 0:
                print('Checked %d' %lp)

    return ndx, ndy, new_smc, new_bounds
            

def removeCells(smccells, midcellsx, midcellsy, cellbounds, celldepths,
                celltype='dry'):
    '''Use the cell depths markers to remove cells from the SMC arrays'''

    if celltype.lower() == 'dry':
        typeval = 0
        print('[INFO] Removing %s cells from SMC array' %celltype.lower())
    elif celltype.lower() == 'alldry':
        typeval = 0
        celldepths[celldepths[:,2] == -2, 2] = 0
        print('[INFO] Removing dry-tier and dry cells from SMC array')
    elif celltype.lower() == 'tier':
        typeval = -1
        print('[INFO] Removing %s cells from SMC array' %celltype.lower())
    elif celltype.lower() == 'alltier':
        typeval = -1
        celldepths[celldepths[:,2] == -2, 2] = -1
        print('[INFO] Removing dry-tier and tier cells from SMC array')

    inds = np.where(celldepths[:,2] != typeval)
    print('[INFO] %d of %d cells kept'
          %tuple([np.size(inds[0]),np.size(smccells[:,0])]))
    new_smccells = smccells[inds[0],:]
    new_midcellsx = midcellsx[inds[0]]
    new_midcellsy = midcellsy[inds[0]]
    new_cellbounds = cellbounds[inds[0],:]
    new_celldepths = celldepths[inds[0],:]

    return new_smccells, new_midcellsx, new_midcellsy, \
           new_cellbounds, new_celldepths

   
def sortSMC(smccells, midcellsx, midcellsy, cellbounds, celldepths):
    '''Reorder SMC cell arrays; needed before write out'''

    print('[INFO] Sorting SMC cells for write out')
    smctmp = np.array([smccells[:,0], smccells[:,1],
                       smccells[:,2], smccells[:,3]])
    sortinds = np.lexsort(smctmp)
    smccells = smccells[sortinds,:]
    midcellsx = midcellsx[sortinds]
    midcellsy = midcellsy[sortinds]
    cellbounds = cellbounds[sortinds,:]
    celldepths = celldepths[sortinds,:]

    return smccells, midcellsx, midcellsy, cellbounds, celldepths

#--- primary functions ---

def createBasesmc(bathyfile, extents, dx, dy, 
                   name='smc_grid', label='basegrid',
                   mindepth=None, dlim=None, drymin=None, drymax=None, 
                   bathytype='gebco', getpland=True, setadj=True):
    '''Create a basic SMC grid object from a bathymetry file'''

    # generate the smc grid
    basesmc = smcGrid(name=name, label=label)
    if mindepth is not None:    
        basesmc.setDryLimits(mindepth=mindepth)
    if dlim is not None:    
        basesmc.setDryLimits(dlim=dlim)
    if drymin is not None:    
        basesmc.setDryLimits(drymin=drymin)
    if drymax is not None:    
        basesmc.setDryLimits(drymax=drymax)
    basesmc.setExtents(extents[0], extents[1], extents[2], extents[3], dx, dy)
    basesmc.setGridFromRegular()
    basesmc.setCFLCells()
    basesmc.setBounds()

    # read in bathymetry data from standard dataset and set model depths
    lon360 = False
    if np.max(basesmc.cellbounds[:,0]) > 180.0: lon360=True
    if bathytype == 'gebco':
        read_lats, read_lons, read_bathy, elevscale, read_pland = readGEBCO(bathyfile, 
                                                                             getpland=True, lon360=lon360)
    basesmc.setDepths(read_lons, read_lats, read_bathy, pland=read_pland, setadj=setadj)

    # run check on cells to add additional tiered cells upstream of 
    # high latitude cells where tiering is needed for next level
    print('[INFO] Checking for high latitude cells that need splitting at next tier')
    inds = np.where(basesmc.celldepths[:,2] == -1)
    tcells = 0
    for cnt, lp in enumerate(inds[0]):
        if np.mod(cnt, 2500) == 0:
            print('[INFO] ... done %d points out of %d' %(cnt, len(inds[0])))
        intersects = chkAdj(lp, basesmc.cellbounds, altbounds=None, nexttier=True)
        if np.any(intersects is not None):
           for chkcell in intersects:
                if chkcell is not None:
                    if basesmc.celldepths[chkcell,2] == 1:
                        if basesmc.smccells[chkcell,2] >=  basesmc.smccells[lp,2]*2:
                            basesmc.celldepths[chkcell,2] = -1
                            tcells = tcells + 1
    print('[INFO] %d cells found' %tcells)

    basesmc.delDryCells()

    return basesmc


def createTiersmc(bathyfile, smcobj,
                   name=None, label='newtier',
                   mindepth=None, dlim=None, drymin=None, drymax=None, 
                   bathytype='gebco', getpland=True, setadj=True, deldry=False):
    '''Create a basic SMC grid object from a bathymetry file'''

    # read in bathymetry data from standard dataset
    lon360 = False
    if np.max(smcobj.cellbounds[:,0]) > 180.0: lon360=True
    if bathytype == 'gebco':
        read_lats, read_lons, read_bathy, elevscale, read_pland = readGEBCO(bathyfile, 
                                                                             getpland=True, lon360=lon360)

    # generate the new tier smc grid object
    if name is None:
        name = smcobj.name
    smctier = smcGrid(name=name, label=label)
    if mindepth is not None:    
        smctier.setDryLimits(mindepth=mindepth)
    if dlim is not None:    
        smctier.setDryLimits(dlim=dlim)
    if drymin is not None:    
        smctier.setDryLimits(drymin=drymin)
    if drymax is not None:    
        smctier.setDryLimits(drymax=drymax)
    smctier.rtd = smcobj.rtd
    smctier.plat = smcobj.plat
    smctier.plon = smcobj.plon
    smctier.llx = smcobj.llx
    smctier.lly = smcobj.lly
    smctier.nx = smcobj.nx
    smctier.ny = smcobj.ny
    smctier.ntiers = smcobj.ntiers + 1

    # define the new tier
    smctier.dx, smctier.dy, smctier.smccells, \
     smctier.cellbounds = smcTier(smcobj.llx, smcobj.lly,
                                   smcobj.dx, smcobj.dy,
                                   smcobj.smccells, smcobj.cellbounds,
                                   celldepths=smcobj.celldepths)

    # generate new cell centres
    smctier.setMids()
 
    # get depth values for the new tier
    smctier.setDepths(read_lons, read_lats, read_bathy, pland=read_pland, setadj=setadj)

    # remove dry cells from the new arrays
    # dont do this for now as want to run cell checking in the join function
    if deldry:
        smctier.delDryCells()
 
    return smctier


def joinTiersmc(basesmc, smctier, bathyfile, tiernext=True):
    '''Combines two smc grids, where second object is next tier'''

    # set up new object 
    # initially this is a copy of the base grid
    newsmc = copy.deepcopy(basesmc)
    print('[INFO] %d cells in base file' %np.shape(newsmc.smccells)[0])
    print('[INFO] %d cells in tier file' %np.shape(smctier.smccells)[0])

    # remove cells marked for tiering from the base arrays
    newsmc.delTierCells(celltype='alltier')
    print('[INFO] %d cells in new base file' %np.shape(newsmc.smccells)[0])

    # run checks on cells to add additional tiered cells next to dry land
    print('[INFO] Checking for points adjacent to tier dry cells in base')
    inds = np.where(smctier.celldepths[:,2] == -2)
    tcells = 0
    for cnt, lp in enumerate(inds[0]):
        if np.mod(cnt, 2500) == 0:
            print('[INFO] ... done %d points out of %d' %(cnt, len(inds[0])))
        intersects = chkAdj(lp, smctier.cellbounds, altbounds=newsmc.cellbounds, nexttier=True)
        if np.any(intersects is not None):
            for chkcell in intersects:
                if chkcell is not None:
                    if newsmc.celldepths[chkcell,2] == 1:
                        newsmc.celldepths[chkcell,2] = -3
                        tcells = tcells + 1
    print('[INFO] %d cells found' %tcells)
    # then, for higher level tiers check how this impacts further cells
    niter = smctier.ntiers - 2
    if niter >0:
        print('[INFO] Checking for higher tier cells that need splitting as a result of dry cell test')
        print('[INFO] Performing %d iterations' %niter)
        for tierlp in range(niter):
            print('[INFO] Performing iteration %d of %d' %(tierlp+1, niter))
            inds = np.where(newsmc.celldepths[:,2] == -3)
            tcells = 0
            for cnt, lp in enumerate(inds[0]):
                if np.mod(cnt, 2500) == 0:
                    print('[INFO] ... done %d points out of %d' %(cnt, len(inds[0])))
                if tierlp > 0:
                    newsmc.celldepths[lp,2] = -1 
                intersects = chkAdj(lp, newsmc.cellbounds, altbounds=None, nexttier=True)
                if np.any(intersects is not None):
                    for chkcell in intersects:
                        if chkcell is not None:
                            if newsmc.celldepths[chkcell,2] == 1:
                                if (newsmc.smccells[chkcell,3] >  newsmc.smccells[lp,3]) or \
                                    (newsmc.smccells[chkcell,2] >=  newsmc.smccells[lp,2]*2):
                                    newsmc.celldepths[chkcell,2] = -3
                                    tcells = tcells + 1
            print('[INFO] %d cells found' %tcells)
    # used -3 above to ensure only newly marked cells are checked against during iterations
    # now set any leftovers to the usual -1 for tiering
    newsmc.celldepths[newsmc.celldepths[:,2] == -3, 2] = -1

    # create a subtier and calculate bathy for additional tiered cells
    print('[INFO] Creating additional tier cells and revised depths')
    subsmc = createTiersmc(bathyfile, newsmc,
                            label='subtier',
                            mindepth=newsmc.depthmin, dlim=newsmc.dlim, 
                            drymin=newsmc.drymin, drymax=newsmc.drymax, 
                            bathytype='gebco', getpland=True, setadj=True)
    # remove any cells marked for tiering, these occur from assessing
    # the bathymetry at a refined level, but should not be further tiered in
    # as already marked as not for tiering in previous level
    #plotGridsmc(subsmc)
    #subsmc.unmarkCells(marker='tier')
    #subsmc.delDryCells()
    #plotGridsmc(subsmc)

    # now check for any jumps this may have caused???

    # remove new cells marked for subtiering from the base arrays
    newsmc.delTierCells(celltype='alltier')
    print('[INFO] %d cells in new base file' %np.shape(newsmc.smccells)[0])
    
    # update new object metadata using smctier
    # no update to name, llx, lly and rotated pole as should be same across all tiers
    print('[INFO] Updating combined grid metadata to use lower tier')
    newsmc.label = basesmc.label + '-' + smctier.label
    newsmc.ntiers = smctier.ntiers
    newsmc.depthmin = smctier.depthmin
    newsmc.dlim = smctier.dlim
    newsmc.drymin = newsmc.drymin
    newsmc.drymax = newsmc.drymax
    newsmc.dx = smctier.dx
    newsmc.dy = smctier.dy
    newsmc.nx = smctier.nx
    newsmc.ny = smctier.ny
    # adjust smc cell index values for base grid prior to combination
    newsmc.smccells = newsmc.smccells * 2

    print('[INFO] Combining grid cell and depth arrays')
    # append the subtier data to the original arrays
    newsmc.smccells = np.concatenate([newsmc.smccells, subsmc.smccells])
    newsmc.midcellsx = np.concatenate([newsmc.midcellsx, subsmc.midcellsx])
    newsmc.midcellsy = np.concatenate([newsmc.midcellsy, subsmc.midcellsy])
    newsmc.cellbounds = np.concatenate([newsmc.cellbounds, subsmc.cellbounds])
    newsmc.celldepths = np.concatenate([newsmc.celldepths, subsmc.celldepths])

    # remove dry cells from the tier array - no longer needed
    smctier.delDryCells()

    # append the new data to the original arrays
    newsmc.smccells = np.concatenate([newsmc.smccells, smctier.smccells])
    newsmc.midcellsx = np.concatenate([newsmc.midcellsx, smctier.midcellsx])
    newsmc.midcellsy = np.concatenate([newsmc.midcellsy, smctier.midcellsy])
    newsmc.cellbounds = np.concatenate([newsmc.cellbounds, smctier.cellbounds])
    newsmc.celldepths = np.concatenate([newsmc.celldepths, smctier.celldepths])
    print('[INFO] %d cells in new base file' %np.shape(newsmc.smccells)[0])

    # run second check on cells to add additional tiered cells upstream of 
    # cells where tiering is maintained for next level...
    if tiernext:
        print('[INFO] Checking for higher tier cells that need splitting at next tier')
        niter = newsmc.ntiers - 1
        print('[INFO] Performing %d iterations' %niter)
        for tierlp in range(niter):
            print('[INFO] Performing iteration %d of %d' %(tierlp+1, niter))
            if tierlp == 0:
                inds = np.where(newsmc.celldepths[:,2] < 0)
            else:
                inds = np.where(newsmc.celldepths[:,2] == -3)
            tcells = 0
            for cnt, lp in enumerate(inds[0]):
                if np.mod(cnt, 2500) == 0:
                    print('[INFO] ... done %d points out of %d' %(cnt, len(inds[0])))
                if tierlp > 0:
                    newsmc.celldepths[lp,2] = -1 
                intersects = chkAdj(lp, newsmc.cellbounds, altbounds=None, nexttier=True)
                if np.any(intersects is not None):
                    for chkcell in intersects:
                        if chkcell is not None:
                            if newsmc.celldepths[chkcell,2] == 1:
                                if (newsmc.smccells[chkcell,3] >  newsmc.smccells[lp,3]) or \
                                    (newsmc.smccells[chkcell,2] >=  newsmc.smccells[lp,2]*2):
                                    newsmc.celldepths[chkcell,2] = -3
                                    tcells = tcells + 1
            print('[INFO] %d cells found' %tcells)
        # used -3 above to ensure only newly marked cells are checked against during iterations
        # now set any leftovers to the usual -1 for tiering
        newsmc.celldepths[newsmc.celldepths[:,2] == -3, 2] = -1

    return newsmc


def loadNCsmc(smcfile):
    '''Load an existing SMC grid from netCDF file'''

    print('[INFO] Loading existing SMC grid from %s' %smcfile) 

    smcobj = smcGrid()

    d = nc.Dataset(smcfile)
    smcobj.name = d.variables['grid_name'][0]
    smcobj.label = d.variables['label'][0]
    smcobj.ntiers = d.variables['ntiers'][0]
    if d.variables['rtd'][0] == 1:
        smcobj.rtd = True
    else:
        smcobj.rtd = False
    smcobj.plat = d.variables['plat'][0]
    smcobj.plon = d.variables['plon'][0]
    smcobj.depthmin = d.variables['depthmin'][0]
    smcobj.dlim = d.variables['dlim'][0]
    smcobj.drymin = d.variables['drymin'][0]
    smcobj.drymax = d.variables['drymax'][0]
    smcobj.llx = d.variables['llx'][0]
    smcobj.lly = d.variables['lly'][0]
    smcobj.dx = d.variables['dx'][0]
    smcobj.dy = d.variables['dy'][0]
    smcobj.nx = d.variables['nx'][0]
    smcobj.ny = d.variables['ny'][0]
    smcobj.smccells = d.variables['smccells'][:,:]
    smcobj.celldepths = d.variables['celldepths'][:,:]
    smcobj.setBounds()
    smcobj.setMids()

    return smcobj


#--- read-write functions ---

def readGEBCO(filein, getpland=True, lon360=False):

    elevscale = -1.0 #GEBCO water depths are -ve

    print('[INFO] Reading file derived from GEBCO dataset')
    print('[INFO] Reading baseline bathy data from: %s' %filein)

    d = nc.Dataset(filein)
    lats  = d.variables['lat'][:]
    print('[INFO] Reading %d latitude points' %len(lats))
    lons  = d.variables['lon'][:]
    print('[INFO] Reading %d longitude points' %len(lons))
    elev  = d.variables['elevation'][:]
    if getpland:
        print('[INFO] Getting existing land mask information')
        pland = d.variables['landmask'][:]
    else:
        pland = None
    d.close()

    if elevscale < 0.:
        print('[INFO] Bathy convention is depth negative')
    else:
        print('[INFO] Bathy convention is depth positive - resetting to depth negative')
        elev = elev * elevscale * -1.        
        elevscale = elevscale * -1.

    if lon360:
        print('[INFO] Shifting arrays to 0->360 degree reference frame')
        lons[lons<0.0] = 360+lons[lons<0.0]
        rollind = np.where(lons == np.min(lons))[0][0]
        lons = np.roll(lons,-1*rollind)
        elev = np.roll(elev,-1*rollind,axis=1)
        if pland is not None:
            pland = np.roll(pland,-1*rollind,axis=1)

    return lats, lons, elev, elevscale, pland


def writeWW3regular(writedepths, writemask, writeblock,
                    llx, lly, dx, dy, nx, ny, 
                    llscale=1.0, ldscale=1.0, 
                    depthlim=0.0, mindepth=10.0,
                    depscale=1.0, blkscale=100.0, 
                    rtd=False, plat=90.0, plon=0.0,
                    writedir='.'):
    '''Write out a regular grid to WAVEWATCH III grid arrays and
       metadata files'''

    WW3Bathy  = writedir+'/ww3depths.txt'
    unitbathy = 30
    idlabathy = 3
    idfmbathy = 1

    WW3Mask  = writedir+'/ww3mask.txt'
    unitmask = 33
    idlamask = 3
    idfmmask = 1

    WW3Block  = writedir+'/ww3block.txt'
    unitblock = 31
    idlablock = 3
    idfmblock = 1

    WW3Meta  = writedir+'/ww3meta.txt'
    WW3GDef  = writedir+'/ww3.grid_def'

    # write out the depth file and scaling value
    # use flipud assuming we read in from bottom to
    # top and want to write out the same way
    print('[INFO] Writing depth info to '+WW3Bathy)
    np.savetxt(WW3Bathy, np.flipud(writedepths)*depscale,
               delimiter=" ", fmt="%i")

    # write out the depth file and scaling value
    print('[INFO] Writing mask info to '+WW3Mask)
    np.savetxt(WW3Mask, np.flipud(writemask),
               delimiter=" ", fmt="%i")

    # write out the blocking file
    print('[INFO] Writing blocking info to '+WW3Block)
    np.savetxt(WW3Block, np.flipud(writeblock)*blkscale,
               delimiter=" ", fmt="%2i")

    # write info to metadata file
    mdx  = dx * ldscale
    mdy  = dy * ldscale
    mllx = llx * llscale
    mlly = lly * llscale

    # estimate limits on CFL and 2nd order swell age
    maxlat  = llx + np.float(ny) * dy
    minlon  = 1853.0 * 60.0 * dx * np.cos(np.pi*maxlat/180.0)
    maxcg   = 1.4 * 9.81 * 25.0 / (4.0 * np.pi)
    cflstep = minlon / maxcg
    sagemax = 0.5 * minlon**2.0 * 12.0 / \
              ((2.0*np.pi*maxcg/24.0)**2.0 * cflstep)

    print('[INFO] Writing WW3 metadata to '+WW3Meta)
    with open(WW3Meta,'w') as inp:
        inp.write('$ Grid minimum cell dx: %.2f' %minlon +
                  'm at latitude %.3f' %maxlat +' degrees\n')
        inp.write('$ CFL minimum timestep (needs rounding down): %i' %cflstep + \
                  ' seconds\n')
        inp.write('$ Estimated maximum swell age for 24 direction ' + \
                  'spectrum: %i' %sagemax +' seconds\n')
        inp.write('$ Minimum depth set for model at %.1f' %mindepth + 'm\n')
        if rtd:
            inp.write('$ Grid specified on rotated pole with plat = %.1f' %plat + \
                      ', plot = %.1f\n' %plon)
        inp.write('$\n')
        inp.write('$ Define grid rules -------------------------------------------------- $\n')
        inp.write('$ Four records containing :\n')
        inp.write('$  1 NX, NY. As the outer grid lines are always defined ' + \
                  'as land\n')
        inp.write('$    points, the minimum size is 3x3.\n')
        inp.write('$  2 Grid increments SX, SY (degr.or m) and scaling ' + \
                  '(division) factor.\n')
        inp.write('$    If NX*SX = 360., latitudinal closure is applied.\n')
        inp.write('$  3 Coordinates of (1,1) (degr.) and scaling (division) ' + \
                  'factor.\n')
        inp.write('$  4 Limiting bottom depth (m) to discriminate between ' + \
                  'land and sea\n')
        inp.write('$    points, minimum water depth (m) as allowed in model, ' + \
                  'unit number\n')
        inp.write('$    of file with bottom depths, scale factor for bottom ' + \
                  'depths (mult.),\n')
        inp.write('$    IDLA, IDFM, format for formatted read, FROM and ' + \
                  'filename.\n')
        inp.write('$\n')
        inp.write('$ Define grid -------------------------------------------------------- $\n')
        inp.write('$\n')
        inp.write(' %i' %nx + ' %i' %ny +'\n')
        inp.write(' %5.3f' %mdx +' %5.3f' %mdy +' %5.1f' %ldscale +'\n')
        inp.write(' %5.3f' %mllx +' %5.3f' %mlly + ' %5.1f' %llscale +'\n')
        inp.write(' %5.2f' %depthlim +' %5.1f' %mindepth +' %i' %unitbathy + \
                  ' %5.1f' %depscale +' %i' %idlabathy +' %i' %idfmbathy + \
                  " '(....)' 'NAME' '%s" %WW3Bathy +"'\n")
        inp.write('$\n')
        inp.write('$ Sub-grid blocking input file\n') 
        inp.write(' %i' %unitblock +' %5.1f' %blkscale + \
                  ' %i' %idlablock +' %i' %idfmblock + \
                  " '(....)' 'NAME' '%s" %WW3Block +"'\n")  
        inp.write('$\n')
        inp.write('$ Input boundary points and excluded points -------------------------- $\n') 
        inp.write('$\n')
        inp.write(' %i' %unitmask +' %i' %idlamask +' %i' %idfmmask + \
                  " '(....)' 'NAME' '%s" %WW3Mask +"'\n")  
        inp.close()

    # write grid data to grid_def file
    print('[INFO] Writing grid_def metadata to '+WW3GDef)
    with open(WW3GDef,'w') as inp:
        inp.write(' %i' %nx + ' %i' %ny +'\n')
        inp.write(' %5.3f' %llx +' %5.3f' %lly + ' %8.6f' %dx +' %8.6f' %dy +'\n')
        inp.write(' %.1f' %plon + ' %.1f\n' %plat)
        inp.close()
    

def writeWW3smc(smccells, celldepths, ntiers,
                llx, lly, dx, dy, nx, ny, 
                llscale=1.0, ldscale=1.0, 
                depthlim=0.0, mindepth=10.0, 
                writemindepth=False,
                depscale=1.0, blkscale=100.0, 
                rtd=False, plat=90.0, plon=0.0,
                writedir='.'):
    '''Write out a regular grid to WAVEWATCH III grid arrays and
       metadata files'''

    WW3Cels = writedir+'/ww3Cels.dat'
    WW3Obst = writedir+'/ww3Obstr.dat'
    unitbathy = 30
    idlabathy = 3
    idfmbathy = 1
    WW3Meta  = writedir+'/smc.ww3meta.txt'
    WW3GDef  = writedir+'/smc.ww3.grid_def'

    # write out the cells file
    print('[INFO] Writing cell info to '+WW3Cels)
    with open(WW3Cels,'w') as inp:
        ncells = np.shape(smccells)[0]
        outcells = ' %d' %ncells
        for lp in range(ntiers):
            ntcells = np.size(np.where(smccells[:,3]==np.int(2**lp)))
            outcells = outcells + ' %d' %ntcells
        inp.write(outcells+'\n')
        for lp in range(ncells):
            # always write out positive cell depths for now...
            if writemindepth:
                inp.write(' %5d %5d %2d %2d ' %tuple(smccells[lp,:]) + \
                           '%5d\n' %np.max([np.abs(celldepths[lp,0]),mindepth]))
            else:
                inp.write(' %5d %5d %2d %2d ' %tuple(smccells[lp,:]) + '%5d\n' %np.abs(celldepths[lp,0]))
        inp.close()

    # write out the obstruction file
    print('[INFO] Writing obstruction info to '+WW3Obst)
    with open(WW3Obst,'w') as inp:
        ncells = np.shape(smccells)[0]
        outcells = ' %d' %ncells
        inp.write(outcells+'   1\n')
        for lp in range(ncells):
            # write out integer cell land percentages...
            landpc = np.int(celldepths[lp,1]*100.0)
            if landpc < 1:
                landpc = 0
            inp.write(' %2d\n' %landpc)
        inp.close()

    # write info to metadata file

    # calculating output metadata here
    mdx = 2.0**(ntiers-1) * dx / ldscale # values for grid_def file - based on largest smc cell size
    mdy = 2.0**(ntiers-1) * dy / ldscale # values for grid_def file - based on largest smc cell size
    mllx = llx * llscale + mdx * ldscale / 2.0 # grid centres for ll cell - based on largest smc cell
    mlly = lly * llscale + mdy * ldscale / 2.0 # grid centres for ll cell - based on largest smc cell

    # calculate limits on CFL and 2nd order swell age - based on largest smc cell size
    maxlat  = (lly / llscale) + np.float(ny) * (dy / ldscale)
    minlon  = 1853.0 * 60.0 * (2.0**(ntiers-1) * dx / ldscale) * np.cos(np.pi*maxlat/180.0)
    maxcg   = 1.4 * 9.81 * 25.0 / (4.0 * np.pi)
    cflstep = minlon / maxcg
    sagemax = 0.5 * minlon**2.0 * 12.0 / ((2.0*np.pi*maxcg/24.0)**2.0 * cflstep)

    # write grid data to grid.inp metadata file
    # note grid parameters are defined by the samllest cell size
    # and use the small cell centre for the sw corner
    print('[INFO] Writing WW3 metadata to '+WW3Meta)
    with open(WW3Meta,'w') as inp:
        inp.write('$ Grid minimum cell dx: %.2f' %minlon + \
                  'm at latitude %.3f' %maxlat +' degrees\n')
        inp.write('$ CFL minimum timestep (needs rounding down): %i' %cflstep + \
                  ' seconds\n')
        inp.write('$ Estimated maximum swell age for 24 direction spectrum: %i' %sagemax \
                  +' seconds\n')
        inp.write('$ Minimum depth set for model at %f' %mindepth + 'm\n')
        if rtd:
            inp.write('$ Grid specified on rotated pole with plat = %.1f' %plat + \
                      ', plot = %.1f\n' %plon)
        inp.write('$\n')
        inp.write('$ Define grid rules -------------------------------------------------- $\n')
        inp.write('$ Four records containing :\n')
        inp.write('$  1 NX, NY. As the outer grid lines are always defined' + \
                  ' as land\n')
        inp.write('$    points, the minimum size is 3x3.\n')
        inp.write('$  2 Grid increments SX, SY (degr.or m) and scaling ' + \
                  '(division) factor.\n')
        inp.write('$    If NX*SX = 360., latitudinal closure is applied.\n')
        inp.write('$  3 Coordinates of (1,1) (degr.) and scaling (division) ' + \
                  'factor.\n')
        inp.write('$  4 Limiting bottom depth (m) to discriminate between ' + \
                  'land and sea\n')
        inp.write('$    points, minimum water depth (m) as allowed in model, ' + \
                  'unit number\n')
        inp.write('$    of file with bottom depths, scale factor for bottom ' + \
                  'depths (mult.),\n')
        inp.write('$    IDLA, IDFM, format for formatted read, FROM and filename.\n')
        inp.write('$\n')
        inp.write('$ Define grid -------------------------------------------------------- $\n')
        inp.write('$\n')
        inp.write(' %i' %nx + ' %i' %ny +'\n')
        # line below adds vn4.18 cell details, ntiers, jshift, nbdys, ishift
        inp.write(' %d' %ntiers +' 0 0 0\n')
        inp.write(' %10.8f' %mdx +' %10.8f' %mdy +' %5.1f' %ldscale +'\n')
        inp.write(' %12.8f' %mllx +' %12.8f' %mlly +' %5.1f' %llscale +'\n')
        inp.write(' %5.2f' %depthlim +' %5.1f' %mindepth +' %i' %unitbathy + \
                  ' %5.1f' %depscale +' %i' %idlabathy +' %i' %idfmbathy + \
                  " '(....)' 'NAME' '%s" %WW3Cels +"'\n")
        inp.write('$\n')
        inp.close()

    # write grid data to grid_def file
    # note grid_def parameters for pre-procesing are defined by the largest cell size
    # and use a largest cell centre for the sw corner
    #nxdef = np.int(nx / 2.0**(ntiers-1))
    #nydef = np.int(ny / 2.0**(ntiers-1))
    print('[INFO] Writing grid_def metadata to '+WW3GDef)
    with open(WW3GDef,'w') as inp:
        inp.write(' %i' %nx + ' %i' %ny +'\n')
        inp.write(' %12.8f' %mllx +' %12.8f' %mlly +' %10.8f' %mdx +' %10.8f' %mdy +'\n')
        inp.write(' %6.2f' %plon +' %6.2f' %plat +'\n')
        inp.close()


def writeBoundsmc(smcobj, writedir='.', lon360=True, 
                    north=False, east=False, south=False, west=False):
    '''Writes boundary cell centre and index values'''

    # not sure if the max values in below are correct, may need to be largest cells only??
    chkbounds = [[north, 1, np.max(smcobj.smccells[:,1])],
                 [east, 0, np.max(smcobj.smccells[:,0])],
                 [south, 1, 0], [west, 0, 0]]

    outfile = writedir + '/ww3Bounds.dat'
    print('[INFO] Writing boundary data to: '+outfile)
    with open(outfile, 'w') as inp:
        for chk in chkbounds:
            if chk[0]: 
                inds = np.where(smcobj.smccells[:,chk[1]] == chk[2])
                blonlat = np.transpose(np.array([smcobj.midcellsx[inds], smcobj.midcellsy[inds]]))
                if chk[1] == 0:
                    sortinds = np.lexsort(np.array([blonlat[:,0],blonlat[:,1]]))
                else:
                    sortinds = np.lexsort(np.array([blonlat[:,1],blonlat[:,0]]))
                blonlat = blonlat[sortinds,:]
                if lon360:
                    blonlat[blonlat[:,0]<0.0,0] = blonlat[blonlat[:,0]<0.0,0] + 360.0
                for ind in range(np.shape(blonlat)[0]):     
                    inp.write('%12.7f %12.7f 0.0 0.0 1\n' %(blonlat[ind,0],blonlat[ind,1]))
        inp.close()

    outfile = writedir + '/ww3BoundCels.dat'
    print('[INFO] Writing boundary cell values to: '+outfile)
    with open(outfile, 'w') as inp:
        for chk in chkbounds:
            if chk[0]: 
                inds = np.where(smcobj.smccells[:,chk[1]] == chk[2])
                bcels = np.array(smcobj.smccells[inds[0],:])
                for lp in range(np.shape(bcels)[0]):     
                    inp.write(' %5d %5d %2d %2d\n' %tuple([bcels[lp,0],bcels[lp,1],bcels[lp,2],bcels[lp,3]]))
        inp.close()

    return


def writeNCsmc(smcobj, writedir='.', filename=None):
    '''Writes core SMC grid and metadata to a netCDF file'''

    if filename is None:
        filename = smcobj.name.replace(' ','') + '_' + smcobj.label.replace(' ','') + \
                    '.nc'
    ncfile = writedir+'/'+filename
    with nc.Dataset(ncfile,'w') as outf:
        ncells = outf.createDimension('seapoints', size=np.shape(smcobj.smccells)[0])
        ninds = outf.createDimension('cellindex', size=np.shape(smcobj.smccells)[1])
        ndeps = outf.createDimension('depthtype', size=np.shape(smcobj.celldepths)[1])

        ncvar = outf.createVariable('smccells','i4',dimensions=('seapoints','cellindex'))
        ncvar.units = '1'
        ncvar.long_name = 'SMC cells layout'
        ncvar.comment = 'x index, y index, x cell size, y cell size'
        ncvar[:,:] = smcobj.smccells[:,:]

        ncvar = outf.createVariable('celldepths','f4',dimensions=('seapoints','depthtype'))
        ncvar.units = 'm'
        ncvar.long_name = 'depths and land'
        ncvar.comment = 'cell depths, normalised land factors and tier markers (1,-1)'
        ncvar[:,:] = smcobj.celldepths[:,:]

        ncvar = outf.createVariable('grid_name',np.str)
        ncvar.long_name = 'grid name'
        ncvar[0] = smcobj.name

        ncvar = outf.createVariable('label',np.str)
        ncvar.long_name = 'grid label'
        ncvar[0] = smcobj.label

        ncvar = outf.createVariable('ntiers','i2')
        ncvar.long_name = 'number of tiers'
        ncvar[:] = smcobj.ntiers

        ncvar = outf.createVariable('rtd','i2')
        ncvar.long_name = 'rotated pole logical'
        ncvar.comment = '1 for true, 0 for false'
        if smcobj.rtd:
            ncvar[:] = 1
        else:
            ncvar[:] = 0

        ncvar = outf.createVariable('plat','f4')
        ncvar.long_name = 'rotated pole latitude'
        ncvar[:] = smcobj.plat

        ncvar = outf.createVariable('plon','f4')
        ncvar.long_name = 'rotated pole longitude'
        ncvar[:] = smcobj.plon

        ncvar = outf.createVariable('depthmin','f4')
        ncvar.long_name = 'model minimum depth'
        ncvar[:] = smcobj.depthmin

        ncvar = outf.createVariable('dlim','f4')
        ncvar.long_name = 'land-sea depth threshold'
        ncvar[:] = smcobj.dlim

        ncvar = outf.createVariable('drymin','f4')
        ncvar.long_name = 'land coverage threshold for defining wet cells'
        ncvar[:] = smcobj.drymin

        ncvar = outf.createVariable('drymax','f4')
        ncvar.long_name = 'land coverage threshold for defining dry cells'
        ncvar[:] = smcobj.drymax

        ncvar = outf.createVariable('llx','f4')
        ncvar.long_name = 'longitude of lower left corner'
        ncvar[:] = smcobj.llx

        ncvar = outf.createVariable('lly','f4')
        ncvar.long_name = 'latitude of lower left corner'
        ncvar[:] = smcobj.lly

        ncvar = outf.createVariable('dx','f4')
        ncvar.long_name = 'grid longitude delta'
        ncvar[:] = smcobj.dx

        ncvar = outf.createVariable('dy','f4')
        ncvar.long_name = 'grid latitude delta'
        ncvar[:] = smcobj.dy

        ncvar = outf.createVariable('nx','i4')
        ncvar.long_name = 'number of x cells at base resolution'
        ncvar[:] = smcobj.nx

        ncvar = outf.createVariable('ny','i4')
        ncvar.long_name = 'number of y cells at base resolution'
        ncvar[:] = smcobj.ny

    ncfile = writedir+'/'+filename
    return ncfile


##--- visualization

def plotGridsmc(smcobj, latlon=True):
    '''Visualise the grid depths, dry cells and cell type'''

    if latlon:
        gx = smcobj.midcellsx
        gy = smcobj.midcellsy
    else:
        gx = smcobj.smccells[:,0]
        gy = smcobj.smccells[:,1]

    plt.subplot(2,2,1)
    plt.scatter(gx, gy, c=smcobj.celldepths[:,0],
                 vmin=-500., vmax=10., edgecolor='None', s=10)
    plt.colorbar()
    plt.title('Gridded depth')

    plt.subplot(2,2,2)
    plt.scatter(gx, gy, c=smcobj.celldepths[:,1],
                 vmin=0.0, vmax=0.5, edgecolor='None', s=10)
    plt.colorbar()
    plt.title('Percentage dry')

    plt.subplot(2,2,3)
    plt.scatter(gx, gy, c=smcobj.celldepths[:,2], edgecolor='None', s=10)
    plt.colorbar()
    plt.title('Cell type')

    plt.subplot(2,2,4)
    plt.scatter(gx, gy, c=np.log2(smcobj.smccells[:,3]), edgecolor='None', s=10)
    plt.colorbar()
    plt.title('Cell factor (y)')

    plt.show()
    plt.close()
    return
