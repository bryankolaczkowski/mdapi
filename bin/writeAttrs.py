#!/usr/bin/env python3

import sys
from optparse import OptionParser, OptionGroup

### set up command-line argument parsing ###
optparser = OptionParser(usage="usage: %prog [options] FILE1.resinteractions.csv [FILE2.resinteractions.csv ... ]",
                         version="%prog v0.1",
                         description="process mdapi residue interactions file(s) for attribute definition in chimera")

optparser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                     help="print some runtime info to the screen [default: %default]")
optparser.add_option("-q", "--quiet", action="store_false", dest="verbose",
                     help="no news is (probably) good news")

optparser.set_defaults(verbose=True)

(options, args) = optparser.parse_args()

if len(args) < 1:
    optparser.error("incorrect number of arguments: at least one FILE.resinteractions.csv is required")


### process input files ###
data_dictionary = {}

for infname in args:

    handle = open(infname, "r")
    datatypes = handle.readline().strip().split(",")[1:]
    for dtype in datatypes:
        if dtype not in data_dictionary.keys():
            data_dictionary[dtype] = {}

    for line in handle:
        linearr = line.strip().split(",")
        resid  = linearr[0]
        resnum = int(resid[3:])
        values = [float(x) for x in linearr[1:]]
        for i in range(len(values)):
            dtype = datatypes[i]
            value = values[i]
            if resnum not in data_dictionary[dtype].keys():
                data_dictionary[dtype][resnum] = [value]
            else:
                data_dictionary[dtype][resnum].append(value)

    handle.close()

### generate output files ###
for datatype in data_dictionary.keys():
    outfname = datatype + ".attr"
    outf = open(outfname, "w")
    outf.write("attribute: %s\n" % datatype.lower())
    outf.write("match mode: 1-to-1\n")
    outf.write("recipient: residues\n")

    reskeys = list(data_dictionary[datatype].keys())
    reskeys.sort()
    for resid in reskeys:
        val = sum(data_dictionary[datatype][resid]) / len(data_dictionary[datatype][resid])
        outf.write("\t:%s\t%.4f\n" % (resid,val))
    outf.close()
