#!/usr/bin/env python3

import sys
import numpy
from scipy.interpolate import griddata
from scipy.spatial import KDTree

################################################################################
# WRITE QMAP GRID TO .DX FORMAT                                                #


## READGRID FUNCTION ----------------------------------------------------------#

def readGrid(infname):
    ### parse qmap grid file ###
    xindex =  0
    yindex =  1
    zindex =  2
    dindex =  3
    vindex = -1

    orig_gridpts = []
    orig_values  = []

    incx   = incy   = incz = -1.0
    startx = starty = startz = -1.0
    endx   = endy   = endz   = -1.0

    handle = open(infname, "r")
    for line in handle:
        linearr = line.split()
        if linearr[0][0] == "#":
            if len(linearr) > 2:
                # parse comment to get starting x,y,z #
                if linearr[1] == "MIN" and linearr[2] == "X,Y,Z":
                    startx = float(linearr[3])
                    starty = float(linearr[4])
                    startz = float(linearr[5])
                # parse comment to get ending x,y,z   #
                elif linearr[1] == "MAX" and linearr[2] == "X,Y,Z":
                    endx = float(linearr[3])
                    endy = float(linearr[4])
                    endz = float(linearr[5])
                # parse comment to get grid increment #
                elif linearr[1] == "GRID" and linearr[2] == "SPACING":
                    incx = incy = incz = float(linearr[3])
        # parse grid data to get x,y,z and total potential energy difference #
        else:
            orig_gridpts.append([float(linearr[xindex]),
                                 float(linearr[yindex]),
                                 float(linearr[zindex])])
            orig_values.append(float(linearr[vindex]))
    handle.close()

    return (orig_gridpts, orig_values, startx,starty,startz, endx,endy,endz, incx,incy,incz)

## DONE READGRID FUNCTION -----------------------------------------------------#


## INTERPOLATEGRID FUNCTION ---------------------------------------------------#

def interpolateGrid(orig_gridpts, orig_values, startx,starty,startz, endx,endy,endz, incx,incy,incz):
    # generate new regular grid points #
    new_gridpts = []

    x  = startx
    nx = 0
    while x <= endx:
        y  = starty
        ny = 0
        while y <= endy:
            z  = startz
            nz = 0
            while z <= endz:
                new_gridpts.append([x,y,z])
                z  += incz
                nz += 1
            y  += incy
            ny += 1
        x  += incx
        nx += 1

    ### interpolate data on regular x,y,z grid ###
    val_grid = griddata(orig_gridpts, orig_values, new_gridpts, method="linear", fill_value=0.0)

    return (new_gridpts, val_grid, nx,ny,nz)

## DONE INTERPOLATEGRID FUNCTION ----------------------------------------------#


## WRITEGRID FUNCTION ---------------------------------------------------------#

def writeGrid(outfname, startx,starty,startz, incx,incy,incz, nx,ny,nz, values):

    outf = open(outfname, "w")

    outf.write("object 1 class gridpositions counts %d %d %d\n" % (nx,ny,nz))
    outf.write("origin %.2f %.2f %.2f\n" % (startx,starty,startz))
    outf.write("delta %.2f 0.00 0.00\n" % incx)
    outf.write("delta 0.00 %.2f 0.00\n" % incy)
    outf.write("delta 0.00 0.00 %.2f\n" % incz)
    outf.write("object 2 class gridconnections counts %d %d %d\n" % (nx,ny,nz))
    outf.write("object 3 class array type double rank 0 items n data follows\n")

    for v in values:
        outf.write("%.4e\n" % v)

    outf.write('attribute "dep" string "positions"\n')
    outf.write('object "regular positions regular connections" class field\n')
    outf.write('component "positions" value 1\n')
    outf.write('component "connections" value 2\n')
    outf.write('component "data" value 3\n')

    outf.close()

    return

## DONE WRITEGRID FUNCTION ----------------------------------------------------#



if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.stderr.write("usage %s infile.qmapout\n" % sys.argv[0])
        sys.stderr.write("  converts the qmapout file to .dx format for chimera\n")
        sys.exit(1)

    infname  = sys.argv[1]
    outfname = infname + ".dx"

    (orig_gridpts, orig_values, startx,starty,startz, endx,endy,endz, incx,incy,incz) = readGrid(infname)
    (new_gridpts, val_grid, nx,ny,nz) = interpolateGrid(orig_gridpts, orig_values, startx,starty,startz, endx,endy,endz, incx,incy,incz)

    # zero out the values of new grid points that don't have a 'close' original point #
    CLOSE_CUTOFF = 1.5  # angstrom cutoff distance for 'closeness' (H atom = 1.1)
    tr = KDTree(orig_gridpts)
    (dists,neighbors) = tr.query(new_gridpts, k=1, distance_upper_bound=CLOSE_CUTOFF)
    for i in range(len(dists)):
        if dists[i] == numpy.inf:
            val_grid[i] = 0.0   # set energy difference to 0.0 no close point exists
                                # these should be 'inside' the protein, and we don't
                                # want to look at them in chimera

    writeGrid(outfname, startx,starty,startz, incx,incy,incz, nx,ny,nz, val_grid)

    exit(0)
