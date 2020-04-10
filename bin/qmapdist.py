#!/usr/bin/env python3

import sys
import numpy as np
from optparse          import OptionParser, OptionGroup
from scipy.interpolate import griddata
from scipy.ndimage     import gaussian_filter
from scipy.spatial     import KDTree
from _version          import __version__

INTERP_CHOICES = ["nearest", "linear"]

### default values ###
DEF_MIN_DIST     = 1.5   # small distance in angstroms (H radius = 1.1A)
DEF_MAX_DIST     = 6.0   # 'far away' from a protein's surface
DEF_SMALL_ENER   = 0.01  # 'small' energy difference value
DEF_GRID_SIZE    = 1.0   # grid spacing in angstroms
DEF_SMOOTH_SIGMA = 1.0   # standard-deviation for gaussian smoothing
DEF_INTERPOLATE_METHOD = INTERP_CHOICES[0]  # grid interpolation method
DEF_WRITE_DX     = False # write .dx files for viewing with chimera

### function definitions ###

# returns ((x,y,z), potential_energy)
# if (x,y,z) is < min_dist, this is 'inside' the protein, and gets value 0.0
# if (x,y,z) is > max_dist, this is 'far away' from the protein, and also
# gets value 0.0
def get_data_from_file(filename, min_dist=DEF_MIN_DIST, max_dist=DEF_MAX_DIST):
    xyz      = []
    energy   = []

    handle = open(filename, "r")
    for line in handle:
        # skip over comment lines #
        if line[0] == "#":
            continue
        # parse line data #
        linedata = [float(x) for x in line.split()]
        xyzloc = linedata[0:3]
        dist   = linedata[3]
        ener   = linedata[6]
        if dist < min_dist or dist > max_dist:
            ener = 0.0
        # record line data #
        xyz.append(xyzloc)
        energy.append(ener)
    handle.close()

    return (np.array(xyz), np.array(energy))
# done get_data_from_file

# returns (x,y,z,energies)
def get_smoothed_data(orig_points, orig_values,
                      min_x,min_y,min_z, max_x,max_y,max_z,
                      interpolate_method=DEF_INTERPOLATE_METHOD,
                      grid_size=DEF_GRID_SIZE, sigma=DEF_SMOOTH_SIGMA,
                      small_value=DEF_SMALL_ENER, min_dist=DEF_MIN_DIST):

    # interpolate onto a regular grid #
    x = np.arange(min_x, max_x+grid_size, grid_size)
    y = np.arange(min_y, max_y+grid_size, grid_size)
    z = np.arange(min_z, max_z+grid_size, grid_size)

    new_points_lst = []
    for xv in x:
        for yv in y:
            for zv in z:
                new_points_lst.append([xv,yv,zv])

    new_points = np.array(new_points_lst)
    new_values = griddata(orig_points, orig_values, new_points, method=interpolate_method, fill_value=0.0)

    # smooth the regular grid data #
    nv_cube = new_values.reshape((x.size, y.size, z.size))
    smoothed_values_cube = gaussian_filter(nv_cube, sigma=sigma)
    nv_oned = smoothed_values_cube.ravel()

    # clean results #
    tr = KDTree(orig_points)
    dists,neighbors = tr.query(new_points, k=1, distance_upper_bound=min_dist)
    for i in range(len(dists)):
        # set energy to 0.0 if it is already 'very small', #
        # or if no 'close' point exists in the actual data #
        # These may be 'inside' the protein, or very far   #
        # away, and we don't want those values to dominate #
        if abs(nv_oned[i]) < small_value or dists[i] == np.inf:
            nv_oned[i] = 0.0

    return (x,y,z,nv_oned)
# done get_smoothed_data

# writes a .dx file for viewing the volume data in chimera
def write_dx(outfname, x,y,z, values):
    dp = 2

    startx = round(x[0], dp)
    starty = round(y[0], dp)
    startz = round(z[0], dp)
    nx     = x.size
    ny     = y.size
    nz     = z.size
    incx   = round(x[1] - x[0], dp)
    incy   = round(y[1] - y[0], dp)
    incz   = round(z[1] - z[0], dp)

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
# done write_dx

# returns euclidean distance between 2 potential energy distributions
# in 3D, assuming the energy distributions are on the same xyz grid
def qmap_distance(file1, file2,
                  min_dist=DEF_MIN_DIST,
                  max_dist=DEF_MAX_DIST,
                  small_energy=DEF_SMALL_ENER,
                  grid_size=DEF_GRID_SIZE,
                  smooth_sigma=DEF_SMOOTH_SIGMA,
                  interpolate_method=DEF_INTERPOLATE_METHOD,
                  write_dx_file=DEF_WRITE_DX):

    xyz1, energy1 = get_data_from_file(file1, min_dist=min_dist, max_dist=max_dist)
    xyz2, energy2 = get_data_from_file(file2, min_dist=min_dist, max_dist=max_dist)

    min_x = (np.amin(xyz1[:,0]) + np.amin(xyz2[:,0]))/2.0
    max_x = (np.amax(xyz1[:,0]) + np.amax(xyz2[:,0]))/2.0
    min_y = (np.amin(xyz1[:,1]) + np.amin(xyz2[:,1]))/2.0
    max_y = (np.amax(xyz1[:,1]) + np.amax(xyz2[:,1]))/2.0
    min_z = (np.amin(xyz1[:,2]) + np.amin(xyz2[:,2]))/2.0
    max_z = (np.amax(xyz1[:,2]) + np.amax(xyz2[:,2]))/2.0

    x1,y1,z1,e1 = get_smoothed_data(xyz1, energy1,
                                    min_x,min_y,min_z,
                                    max_x,max_y,max_z,
                                    interpolate_method=interpolate_method,
                                    grid_size=grid_size,
                                    sigma=smooth_sigma,
                                    small_value=small_energy,
                                    min_dist=min_dist)

    x2,y2,z2,e2 = get_smoothed_data(xyz2, energy2,
                                    min_x,min_y,min_z,
                                    max_x,max_y,max_z,
                                    interpolate_method=interpolate_method,
                                    grid_size=grid_size,
                                    sigma=smooth_sigma,
                                    small_value=small_energy,
                                    min_dist=min_dist)

    if write_dx_file:
        write_dx(file1+'.dx', x1, y1, z1, e1)
        write_dx(file2+'.dx', x2, y2, z2, e2)

    return np.linalg.norm(e1 - e2)
# done qmap_distance


### main ###

if __name__ == '__main__':

    ### set up command-line argument parsing ###
    optparser = OptionParser(usage="usage: %prog [options] GRID1.qmapout [GRID2.qmapout]",
                             version=__version__,
                             description="calculate distance between two qmap grid outputs")

    optparser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                         help="print some runtime info to the screen [default: %default]")
    optparser.add_option("-q", "--quiet", action="store_false", dest="verbose",
                         help="no news is (probably) good news")

    optparser.add_option("--min-dist", action="store", type="float", dest="min_dist",
                         help="set minimum distance to NUM angstroms [default: %default]",
                         metavar="NUM")
    optparser.add_option("--max-dist", action="store", type="float", dest="max_dist",
                         help="set maximum distance to NUM angstroms [default: %default]",
                         metavar="NUM")

    optparser.add_option("--small-energy", action="store", type="float", dest="small_energy",
                         help="set very-small energy difference to NUM [default: %default]",
                         metavar="NUM")

    optparser.add_option("--grid-size", action="store", type="float", dest="grid_size",
                         help="set grid size to NUM angstroms [default: %default]",
                         metavar="NUM")

    optparser.add_option("--sigma", action="store", type="float", dest="smooth_sigma",
                         help="set gaussian smoothing sigma to NUM [default: %default]",
                         metavar="NUM")

    optparser.add_option("--interpolate", action="store", type="choice", choices=INTERP_CHOICES, dest="interpolate_method",
                         help="set minimum distance to NUM angstroms [default: %default]",
                         metavar="NUM")

    optparser.add_option("-w", "--write", action="store_true", dest="write_dx_file",
                         help="write .dx files for chimera visualization [default: %default]")
    optparser.add_option("-n", "--nowrite", action="store_false", dest="write_dx_file",
                         help="don't write .dx files")


    optparser.set_defaults(verbose=True,
                           min_dist=DEF_MIN_DIST,
                           max_dist=DEF_MAX_DIST,
                           small_energy=DEF_SMALL_ENER,
                           grid_size=DEF_GRID_SIZE,
                           smooth_sigma=DEF_SMOOTH_SIGMA,
                           interpolate_method=DEF_INTERPOLATE_METHOD,
                           write_dx_file=DEF_WRITE_DX)

    (options, args) = optparser.parse_args()

    if len(args) < 1:
        optparser.error("incorrect number of arguments: at least one .qmapout file is required")

    if len(args) > 2:
        optparser.error("incorrect number of arguments: at most two .qmapout files are supported")

    # 2 qmapout files - calculate distance #
    if len(args) == 2:
        fn1 = args[0]
        fn2 = args[1]
        distance = qmap_distance(fn1, fn2,
                                 min_dist=options.min_dist,
                                 max_dist=options.max_dist,
                                 small_energy=options.small_energy,
                                 grid_size=options.grid_size,
                                 smooth_sigma=options.smooth_sigma,
                                 interpolate_method=options.interpolate_method,
                                 write_dx_file=options.write_dx_file)

        sys.stdout.write('%s,%s,%.4f\n' % (fn1,fn2,distance))

    # 1 qmapout file - write .dx for chimera viewing #
    else:
        fn1 = args[0]
        xyz1, energy1 = get_data_from_file(fn1, min_dist=options.min_dist, max_dist=options.max_dist)

        min_x = np.amin(xyz1[:,0])
        max_x = np.amax(xyz1[:,0])
        min_y = np.amin(xyz1[:,1])
        max_y = np.amax(xyz1[:,1])
        min_z = np.amin(xyz1[:,2])
        max_z = np.amax(xyz1[:,2])

        x1,y1,z1,e1 = get_smoothed_data(xyz1, energy1, min_x, min_y, min_z, max_x, max_y, max_z,
                                        interpolate_method=options.interpolate_method,
                                        grid_size=options.grid_size, sigma=options.smooth_sigma,
                                        small_value=options.small_energy, min_dist=options.min_dist)
        write_dx(fn1+'.dx', x1, y1, z1, e1)

    sys.exit(0)
