#!/usr/bin/env python3

import sys
import numpy as np
from optparse      import OptionParser, OptionGroup
from scipy.ndimage import gaussian_filter
from qmapdist      import write_dx
from _version      import __version__

### default values ###
DEF_GRID_SIZE    = 1.0    # grid spacing in angstroms
DEF_SMOOTH_SIGMA = 1.0    # standard-deviation for gaussian smoothing
DEF_NORMSUM      = 1000.0 # normalize counts to sum to this value
DEF_WRITE_DX     = False  # write .dx files for viewing with chimera

### function definitions ###

# returns (xyz) coordinates
def get_data_from_file(filename):
    xyz = []

    handle = open(filename, "r")
    handle.readline()
    for line in handle:
        xyz.append([float(x) for x in line.strip().split(',')[1:]])
    handle.close()

    return np.array(xyz)
# done get_data_from_file

# returns (x,y,z,values)
def get_smoothed_data(orig_points, min_x,min_y,min_z, max_x,max_y,max_z,
                      grid_size=DEF_GRID_SIZE, sigma=DEF_SMOOTH_SIGMA,
                      norm_sum=DEF_NORMSUM):

    # calculate bin count in each dimension #
    nx = round((max_x - min_x) / grid_size)
    ny = round((max_y - min_y) / grid_size)
    nz = round((max_z - min_z) / grid_size)

    # generate count histogram #
    counts,bins = np.histogramdd(orig_points, bins=[nx,ny,nz],
                                 range=((min_x,max_x),(min_y,max_y),(min_z,max_z)),
                                 density=False)

    # normalize counts #
    if norm_sum > 0.0:
        counts = counts / (np.sum(counts)/norm_sum)

    # convert bin_edges to x,y,z coordinates #
    xyz = []
    for i in range(len(bins)):
        b1 = np.delete(bins[i], -1)
        b2 = np.delete(bins[i],  0)
        xyz.append(np.mean(np.array([b1,b2]), axis=0))

    # smooth count histogram #
    smoothed_values_cube = gaussian_filter(counts, sigma=sigma)
    nv_oned = smoothed_values_cube.ravel()

    return (xyz[0],xyz[1],xyz[2],nv_oned)
# done get_smoothed_data

# returns euclidean distance between 2 mdapi distributions in 3D,
# assuming the distributions are on the same xyz grid
def mdapi_distance(file1, file2,
                   grid_size=DEF_GRID_SIZE,
                   sigma=DEF_SMOOTH_SIGMA,
                   norm_sum=DEF_NORMSUM,
                   write_dx_file=DEF_WRITE_DX):

    xyz1 = get_data_from_file(file1)
    xyz2 = get_data_from_file(file2)

    min_x = (np.amin(xyz1[:,0]) + np.amin(xyz2[:,0]))/2.0
    max_x = (np.amax(xyz1[:,0]) + np.amax(xyz2[:,0]))/2.0
    min_y = (np.amin(xyz1[:,1]) + np.amin(xyz2[:,1]))/2.0
    max_y = (np.amax(xyz1[:,1]) + np.amax(xyz2[:,1]))/2.0
    min_z = (np.amin(xyz1[:,2]) + np.amin(xyz2[:,2]))/2.0
    max_z = (np.amax(xyz1[:,2]) + np.amax(xyz2[:,2]))/2.0

    x1,y1,z1,c1 = get_smoothed_data(xyz1, min_x,min_y,min_z, max_x,max_y,max_z,
                                    grid_size, sigma, norm_sum)

    x2,y2,z2,c2 = get_smoothed_data(xyz2, min_x,min_y,min_z, max_x,max_y,max_z,
                                    grid_size, sigma, norm_sum)

    if write_dx_file:
        write_dx(file1+'.dx', x1, y1, z1, c1)
        write_dx(file2+'.dx', x2, y2, z2, c2)

    return np.linalg.norm(c1 - c2)
# done mdapi_distance


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

    optparser.add_option("--grid-size", action="store", type="float", dest="grid_size",
                         help="set grid size to NUM angstroms [default: %default]",
                         metavar="NUM")

    optparser.add_option("--sigma", action="store", type="float", dest="smooth_sigma",
                         help="set gaussian smoothing sigma to NUM [default: %default]",
                         metavar="NUM")

    optparser.add_option("--sum", action="store", type="float", dest="norm_sum",
                         help="normalize counts to sum to NUM [default: %default] Setting to zero turns off normalization",
                         metavar="NUM")

    optparser.add_option("-w", "--write", action="store_true", dest="write_dx_file",
                         help="write .dx files for chimera visualization [default: %default]")
    optparser.add_option("-n", "--nowrite", action="store_false", dest="write_dx_file",
                         help="don't write .dx files")


    optparser.set_defaults(verbose=True,
                           grid_size=DEF_GRID_SIZE,
                           smooth_sigma=DEF_SMOOTH_SIGMA,
                           norm_sum=DEF_NORMSUM,
                           write_dx_file=DEF_WRITE_DX)

    (options, args) = optparser.parse_args()

    if len(args) < 1:
        optparser.error("incorrect number of arguments: at least one mdapi .csv file is required")

    if len(args) > 2:
        optparser.error("incorrect number of arguments: at most two mdapi .csv files are supported")

    # 2 .csv files - calculate distance #
    if len(args) == 2:
        fn1 = args[0]
        fn2 = args[1]
        distance = mdapi_distance(fn1, fn2,
                                  grid_size=options.grid_size,
                                  sigma=options.smooth_sigma,
                                  norm_sum=options.norm_sum,
                                  write_dx_file=options.write_dx_file)

        sys.stdout.write('%s,%s,%.4f\n' % (fn1,fn2,distance))

    # 1 .csv file - write .dx for chimera viewing #
    else:
        fname = args[1]
        xyzdata = get_data_from_file(fname)

        min_x = np.amin(xyzdata[:,0])
        max_x = np.amax(xyzdata[:,0])
        min_y = np.amin(xyzdata[:,1])
        max_y = np.amax(xyzdata[:,1])
        min_z = np.amin(xyzdata[:,2])
        max_z = np.amax(xyzdata[:,2])

        x,y,z,counts = get_smoothed_data(xyzdata, min_x,min_y,min_z, max_x,max_y,max_z,
                                         grid_size=options.grid_size, sigma=options.smooth_sigma,
                                         norm_sum=options.norm_sum)

        write_dx(fname+'.dx', x, y, z, counts)

    sys.exit(0)
