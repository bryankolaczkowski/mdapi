#!/usr/bin/env python3

import sys
from optparse import OptionParser, OptionGroup
import math
import numpy

### set up command-line argument parsing ###
optparser = OptionParser(usage="usage: %prog [options] GRID1.csv [GRID2.csv ... ]",
                         version="%prog v0.1",
                         description="process mdapi grids for visualization or comparison")

optparser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                     help="print some runtime info to the screen [default: %default]")
optparser.add_option("-q", "--quiet", action="store_false", dest="verbose",
                     help="no news is (probably) good news")
optparser.add_option("-g", "--gridsize", action="store", type="float", dest="gridsize",
                     help="set the grid stepsize to NUM angstroms [default: %default]",
                     metavar="NUM")
optparser.add_option("--groups", action="store", type="string", dest="groups",
                     help="set up grid groups for averaging [default: creates a separate group for each input file]. Groups are identified using 'name1(file1,file2,...);name2(file1,file2,...);...nameN(file1,file2,...)' syntax", metavar="GRPSTR")

optparser.set_defaults(verbose=True, gridsize=1.0, groups="")

(options, args) = optparser.parse_args()

if len(args) < 1 and not options.groups:
    optparser.error("incorrect number of arguments: at least one GRID.csv is required, unless you specify grid volume groups using the --groups option")


### organize groups ###
if options.verbose:
    sys.stdout.write("organizing grid volume groups... ")
    sys.stdout.flush()

grid_groups = {}

if not options.groups:
    for infname in args:
        grid_groups[infname] = [infname]

else:
    grparr = options.groups.split(";")
    for grp in grparr:
        grpdefarr = grp.split("(")
        grpname   = grpdefarr[0]
        grpelems  = grpdefarr[1][:-1].split(",")
        grid_groups[grpname] = grpelems

if options.verbose:
    sys.stdout.write("done.\n")
    sys.stdout.write("found %d groups:\n" % len(grid_groups.keys()))
    for (gid,gelems) in grid_groups.items():
        sys.stdout.write("  %s: %s\n" % (gid,gelems))


### pass 1; we're just going to read the min and max coordinates ###
### across all input grid files. also calculate total sum for    ###
### each grid file.                                              ###
### do NOT allow files with "protein" in them to be included in  ###
### the calculation of minimum coordinates; these are the        ###
### locations of protein residues, not fragments!                ###
minx = -sys.float_info.max
miny = -sys.float_info.max
minz = -sys.float_info.max

maxx = sys.float_info.max
maxy = sys.float_info.max
maxz = sys.float_info.max

sums = {}

if options.verbose:
    sys.stdout.write("finding min and max coordinates across all input files...\n")
    donef = 0
    totlf = 0
    for gelems in grid_groups.values():
        totlf += len(gelems)

for elemarr in grid_groups.values():
    for infname in elemarr:

        ## get min and max dims for this file ##

        myminx =  sys.float_info.max
        mymaxx = -sys.float_info.max

        myminy =  sys.float_info.max
        mymaxy = -sys.float_info.max

        myminz =  sys.float_info.max
        mymaxz = -sys.float_info.max

        # also store sum for normalization (if needed) #
        mysum  = 0.0

        handle = open(infname, "r")
        handle.readline()
        for line in handle:
            xyz = [float(x) for x in line.strip().split(",")[1:]]
            mysum += 1.0

            if xyz[0] < myminx:
                myminx = xyz[0]
            if xyz[0] > mymaxx:
                mymaxx = xyz[0]

            if xyz[1] < myminy:
                myminy = xyz[1]
            if xyz[1] > mymaxy:
                mymaxy = xyz[1]

            if xyz[2] < myminz:
                myminz = xyz[2]
            if xyz[2] > mymaxz:
                mymaxz = xyz[2]

        handle.close()

        ## update global min and max dims                      ##
        # shrink box if needed to make sure we have data from   #
        # every input file                                      #
        # make sure we do NOT include protein residue locations #
        # in the global min/max dimensions calculations; for    #
        # now, we'll just make sure 'protein' is not in the     #
        # filename; this is compatible with the mdapi.py naming #
        # conventions for the x,y,z location maps               #
        if infname.find(".protein") == -1:
            if myminx > minx:
                minx = myminx
            if mymaxx < maxx:
                maxx = mymaxx

            if myminy > miny:
                miny = myminy
            if mymaxy < maxy:
                maxy = mymaxy

            if myminz > minz:
                minz = myminz
            if mymaxz < maxz:
                maxz = mymaxz

        sums[infname] = mysum

        if options.verbose:
            donef += 1
            sys.stdout.write("  done %s (%d/%d)\n" % (infname, donef, totlf))


nx = math.floor( (maxx-minx)/options.gridsize )
ny = math.floor( (maxy-miny)/options.gridsize )
nz = math.floor( (maxz-minz)/options.gridsize )

n = nx * ny * nz

startx = minx + (options.gridsize/2.0)
starty = miny + (options.gridsize/2.0)
startz = minz + (options.gridsize/2.0)

incx = options.gridsize
incy = options.gridsize
incz = options.gridsize

if options.verbose:
    sys.stdout.write("done.\n")
    sys.stdout.write("generated %.4f angstrom grid of %d points (x:%d,y:%d,z:%d)\n" % (options.gridsize,n,nx,ny,nz))
    sys.stdout.write("min: (%.4f,%.4f,%.4f); max: (%.4f,%.4f,%.4f)\n" % (minx,miny,minz, maxx,maxy,maxz))


### pass2; all grid variables set constant across all input ###
### files so they can be directly compared, assuming they   ###
### are in the same location in coordinate space.           ###
### next convert the files to .dx format                    ###
if options.verbose:
    sys.stdout.write("writing grid volume files...\n")
    donef = 0
    totlf = len(grid_groups.keys())

for (grpname,grpelems) in grid_groups.items():

    values = []

    for infname in grpelems:
        # read data file #
        data = []
        handle = open(infname, "r")
        handle.readline()
        for line in handle:
            data.append( [float(x) for x in line.strip().split(",")[1:]] )
        handle.close()

        # convert to numpy array and generate histogram #
        dataarray = numpy.array(data)
        (myvalues,binedges) = numpy.histogramdd(dataarray, bins=[nx,ny,nz], range=((minx,maxx),(miny,maxy),(minz,maxz)), density=False)

        # enter data into values array; keep a running sum for averaging #
        if len(values) == 0:
            for x in range(nx):
                for y in range(ny):
                    for z in range(nz):
                        values.append(myvalues[x][y][z])

        else:
            for x in range(nx):
                for y in range(ny):
                    for z in range(nz):
                        values[(x*nz*ny) + (y*nz) +z] += myvalues[x][y][z]

    # average over all fragmaps in the group #
    for i in range(len(values)):
        values[i] /= len(grpelems)

    # write .dx format #
    outfname = grpname + ".dx"
    outf     = open(outfname, "w")

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

    if options.verbose:
        donef += 1
        sys.stdout.write("  wrote %s (%d/%d)\n" % (outfname, donef, totlf))

if options.verbose:
    sys.stdout.write("done.\n")
    sys.stdout.write("processing finished.\n")
