#!/usr/bin/env python3

import sys
import os
import shutil
import glob
import time
import datetime
import math
from optparse import OptionParser, OptionGroup
from _version import __version__
from pqr2pdb  import pqr2pdb

starttime = time.time()

### read relevent environment variables ###

# maybe they set the environment variable #
MDAPIDIR = os.environ.get("MDAPIDIR")

# check the user's path for executables #
GMX      = shutil.which("gmx")
PDB2PQR  = shutil.which("pdb2pqr")
CHARMM   = shutil.which("charmm")
if not MDAPIDIR:
    MDAPIDIR = shutil.which("mdapi.py")
    if MDAPIDIR:
        MDAPIDIR = MDAPIDIR.split("/bin")[0]

# try grepping the program's location from command-line #
if not MDAPIDIR:
    MDAPIDIR = os.path.realpath(os.path.dirname(sys.argv[0]))
    if MDAPIDIR:
        MDAPIDIR = MDAPIDIR.split("/bin")[0]

# empty-out if null #
if not MDAPIDIR:
    MDAPIDIR = ""
if not GMX:
    GMX = ""
if not PDB2PQR:
    PDB2PQR = ""
if not CHARMM:
    CHARMM = ""

# remove trailing "/" if needed #
if len(MDAPIDIR) > 0 and MDAPIDIR[-1] == "/":
    MDAPIDIR = MDAPIDIR[:-1]
if len(GMX) > 0 and GMX[-1] == "/":
    GMX = GMX[:-1]
if len(PDB2PQR) > 0 and PDB2PQR[-1] == "/":
    PDB2PQR = PDB2PQR[:-1]
if len(CHARMM) > 0 and CHARMM[-1] == "/":
    CHARMM = CHARMM[:-1]


### set up command-line argument parsing ###
optparser = OptionParser(usage="usage: %prog [options] PROTFILE.pdb",
                         version=__version__,
                         description="Quick MAPping of protein molecular surface interactions")

optparser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                     help="print some runtime info to the screen [default: %default]")
optparser.add_option("-q", "--quiet", action="store_false", dest="verbose",
                     help="no news is (probably) good news")
optparser.add_option("--logs", action="store_true", dest="logs",
                     help="keep log and intermediate files after program execution [default: %default]")
optparser.add_option("--nologs", action="store_false", dest="logs",
                     help="remove log and intermediate files after program execution")
optparser.add_option("--clean", action="store_true", dest="clean",
                     help="remove unneeded files from the current directory, keeping only the input and results files")

group = OptionGroup(optparser, "Data and Helper Options", "where to find data files and helper programs")
group.add_option("--mdapidir", action="store", type="string", dest="mdapidir",
                 help="set directory containing the MDAPI installation to DIR [default: %default] Note that if you set the MDAPI installation directory this way, you must specify all the force field, topology and radius file options explicitly; you cannot use the default options unless you set the MDAPIDIR environment variable",
                 metavar="DIR")
group.add_option("--pdb2pqr", action="store", type="string", dest="pdb2pqr",
                 help="set name of pdb2pqr execuable program to FILE [default: %default]",
                 metavar="FILE")
group.add_option("--gmx", action="store", type="string", dest="gromacs",
                 help="set name of gromacs gmx executable program to FILE [default: %default]",
                 metavar="FILE")
group.add_option("--gmxforcefield", action="store", type="string", dest="gmxforcefield",
                 help="set directory containing gromacs force field data to DIR [default: %default]",
                 metavar="DIR")

group.add_option("--charmm", action="store", type="string", dest="charmm",
                 help="set name of charmm executable program to FILE [default: %default]",
                 metavar="FILE")
group.add_option("--protrtf", action="store", type="string", dest="charmmprotrtf",
                 help="set name of charmm protein topology to FILE [default: %default]",
                 metavar="FILE")
group.add_option("--protprm", action="store", type="string", dest="charmmprotprm",
                 help="set name of charmm protein force-field parameters to FILE [default: %default]",
                 metavar="FILE")
group.add_option("--probertf", action="store", type="string", dest="charmmprobertf",
                 help="set name of charmm probe topology to FILE [default: %default]",
                 metavar="FILE")
group.add_option("--probeprm", action="store", type="string", dest="charmmprobeprm",
                 help="set name of charmm probe force-field parameters to FILE [default: %default]",
                 metavar="FILE")
group.add_option("--radstr", action="store", type="string", dest="charmmradiusstr",
                 help="set name of charmm implicit solvent radius definitions to FILE [default: %default]",
                 metavar="FILE")

optparser.add_option_group(group)

# options for controlling how grid calculation is set up #
group = OptionGroup(optparser, "Grid Calculation Options", "control the grid calculations")
group.add_option("--translate", action="store", type="string", dest="center",
                 help="translate system's geometric center to X,Y,Z [default: %default]",
                 metavar="X,Y,Z")
group.add_option("--notranslate", action="store_false", dest="notranslate",
                 help="don't translate the system's geometric center [default: %default]")
group.add_option("--pH", action="store", type="float", dest="pH",
                 help="set the pH for protonation calculations to PH [default: %default]",
                 metavar="PH")
group.add_option("--pqroptimize", action="store_true", dest="pqroptimize",
                 help="use pdb2pqr to optimize the protein's hydrogen bond network [default: %default]")
group.add_option("--nopqroptimize", action="store_false", dest="pqroptimize",
                 help="don't allow pdb2pqr to optimize the hydrogen bond network; it will still perform debumping, though")
group.add_option("--margin", action="store", type="float", dest="margin",
                 help="set the grid calculation box margin to MAR angstroms away from the protein's rectangular extents [default: %default]",
                 metavar="MAR")
group.add_option("--adaptmar", action="store", type="float", dest="adaptmar",
                 help="set the grid calculation adaptive margin to MAR angstroms [default: %default]; outside this margin, the gridsize will be doubled",
                 metavar="MAR")
group.add_option("--tightmar", action="store", type="float", dest="tightmar",
                 help="set the grid calculation tight margin to MAR angstroms [default: %default]; inside this margin, the gridsize will be halved",
                 metavar="MAR")
group.add_option("--gridsize", action="store", type="float", dest="gridsize",
                 help="set the grid stepsize to NUM angstroms [default: %default]",
                 metavar="NUM")
group.add_option("--setgrid", action="store", type="string", dest="setgrid",
                 help="manually set the rectangular grid extents to VALS, where VALS=minX,minY,minZ,maxX,maxY,maxZ, a comma-delimited series of minimum and maximum x,y,z coordinates. Note that the margin will be added to these values; you would typically set them to the minimum and maximum x,y,z atom positions over a set of protein structures, so all the grids will line up [default: calculated from the input protein structure]",
                 metavar="VALS")
group.add_option("--solvent", action="store_true", dest="usesolvent",
                 help="use an implicit solvdent model when calculating system energies [default: %default]")
group.add_option("--vacuum", action="store_false", dest="usesolvent",
                 help="don't use the implicit solvent model; calculate system energies in a vacuum")
group.add_option("--soft", action="store_true", dest="softcore",
                 help="use soft-core potential calculations [default: %default]")
group.add_option("--hard", action="store_false", dest="softcore",
                 help="use hard-core potential calculations")
optparser.add_option_group(group)

# set option defaults #
gmxfff   = "charmm36.ff"
protrtf  = "top_all36_prot.rtf"
protprm  = "par_all36_prot.prm"
probertf = "probes.rtf"
probeprm = "probes.prm"
radiusf  = "radius_gbsw.str"

if MDAPIDIR:
    gmxfff   = "%s/data/%s" % (MDAPIDIR,gmxfff)
    protrtf  = "%s/data/%s" % (MDAPIDIR,protrtf)
    protprm  = "%s/data/%s" % (MDAPIDIR,protprm)
    probertf = "%s/data/%s" % (MDAPIDIR,probertf)
    probeprm = "%s/data/%s" % (MDAPIDIR,probeprm)
    radiusf  = "%s/data/%s" % (MDAPIDIR,radiusf)


optparser.set_defaults(verbose=False,
                       logs=False,
                       clean=False,
                       mdapidir=MDAPIDIR,
                       pdb2pqr=PDB2PQR,
                       gromacs=GMX,
                       gmxforcefield=gmxfff,
                       charmm=CHARMM,
                       charmmprotrtf=protrtf,
                       charmmprotprm=protprm,
                       charmmprobertf=probertf,
                       charmmprobeprm=probeprm,
                       charmmradiusstr=radiusf,
                       center="0,0,0",
                       notranslate=False,
                       pH=7.0,
                       pqroptimize=True,
                       margin=8,        # angstrom margin for grid calculation box
                       adaptmar=6,      # outside adaptmar from protein, grid size doubles
                       tightmar=0,      # inside tightmar from protein, grid size halves
                       gridsize=1.0,    # default grid step size
                       setgrid="",
                       usesolvent=True,
                       softcore=True)

(options, args) = optparser.parse_args()

if options.mdapidir != MDAPIDIR:
    MDAPIDIR = options.mdapidir
    if len(MDAPIDIR) > 0 and MDAPIDIR[-1] == "/":
        MDAPIDIR = MDAPIDIR[:-1]

if len(args) < 1:
    optparser.error("incorrect number of arguments: PROTFILE.pdb filename required")
protfilename = args[-1]


### sort out input pdb filename ###
if protfilename.find("./") == 0:
    protfilename = protfilename[2:]
protbasename = protfilename.split("/")[-1].split(".pdb")[0]

# check if input file exists (doesn't check if it's a pdb file) #
if not os.path.isfile(protfilename):
    sys.stderr.write("ERROR: input file %s does not exist.\n" % protfilename)
    sys.exit(1)

# copy protfilename into current directory, if needed #
if protfilename.find("/") != -1:
    shutil.copy(protfilename, "./")
    protfilename = protfilename.split("/")[-1]


if options.gmxforcefield[-1] == "/":
    options.gmxforcefield = options.gmxforcefield[:-1]


### set the files we need in the end ###
gmxffname  = options.gmxforcefield.split("/")[-1]
prodbase   = "%s_qmap" % protbasename
cenfname   = prodbase + ".cen.pdb"
pqrfname   = prodbase + ".pqr"
prtnfname  = prodbase + ".h.charmm.pdb"
gmxpdbfnm  = prodbase + ".gmx.pdb"
gmxtopfnm  = prodbase + ".gmx.top"
chrminfnm  = prodbase + ".pdb"
chrminpfnm = prodbase + ".charmm.inp"

# keep the original protein file    #
# and the charmm input protein file #
my_keep_files = [protfilename, chrminfnm]
# save anything that looks like a results file #
dx_files      = glob.glob("*.dx")
csv_files     = glob.glob("*.csv")
qmapres_files = glob.glob("*.qmapout")
results_files = dx_files + csv_files + qmapres_files
# we'll try to save scripts, as well #
py_files      = glob.glob("*.py")
ba_files      = glob.glob("*.bash")
sh_files      = glob.glob("*.sh")
program_files = py_files + ba_files + sh_files

keep_files    = my_keep_files + results_files + program_files

### check for --clean option ###
if options.clean:
    for f in glob.glob("*"):
        if f not in keep_files:
            if os.path.isfile(f) or os.path.islink(f):
                os.remove(f)
    sys.exit(0)


### check for missing environment variables after option parsing ###
MDAPIDIR = options.mdapidir
if not MDAPIDIR:
    sys.stderr.write("ERROR: MDAPIDIR is not defined\n")
    sys.stderr.write("  please either export MDAPIDIR=\"DIR\" to set the MDAPIDIR environment\n")
    sys.stderr.write("  variable to your MDAPI installation directory, or provide the installation\n")
    sys.stderr.write("  directory using the --mdapidir=DIR option\n")

GMX = options.gromacs
if not GMX:
    sys.stderr.write("ERROR: GROMACS gmx executable is not in your PATH\n")
    sys.stderr.write("  please either upate your PATH environment variable to include your gromacs\n")
    sys.stderr.write("  installation, or provide the loaction of gmx using the --gromacs=GMX option\n")

PDB2PQR = options.pdb2pqr
if not PDB2PQR:
    sys.stderr.write("ERROR: pdb2pqr executable is not in your PATH\n")
    sys.stderr.write("  please either update your PATH environment variable to include your pdb2pqr\n")
    sys.stderr.write("  installation, or provide the location of pdb2pqr using the --pdb2pqr=PDB2PQR option\n")

CHARMM = options.charmm
if not CHARMM:
    sys.stderr.write("ERROR: charmm executable is not in your PATH\n")
    sys.stderr.write("  please either update your PATH environment variable to include your charmm\n")
    sys.stderr.write("  installation, or provide the location of charmm using the --charmm=CHARMM option\n")

if not MDAPIDIR or not GMX or not PDB2PQR or not CHARMM:
    sys.exit(1)


################################################################################
# YAY! WE GET TO DO SOMETHING! - GENERATE CHARMM INPUT                         #

### link force field files to local directory for gromacs ###
if options.verbose:
    sys.stdout.write("linking force field files to local directory for gromacs... ")
    sys.stdout.flush()

if os.path.islink(gmxffname):
    os.remove(gmxffname)
elif os.path.isdir(gmxffname):
    shutil.rmtree(gmxffname)
if not os.path.isdir(options.gmxforcefield):
    sys.stdout.write("ERROR: force field %s not found\n" % (options.gmxforcefield))
    sys.exit(1)
os.symlink(options.gmxforcefield, gmxffname)

if options.verbose:
    sys.stdout.write("done.\n")


## translate system to specified center ##
if options.notranslate:
    shutil.copy(protfilename, cenfpre)
else:
    if options.verbose:
        sys.stdout.write("translating system center to %s... " % options.center)
        sys.stdout.flush()

    (x,y,z) = options.center.split(",")
    cmd = "%s editconf -f %s -o %s -center %s %s %s > gmx.center.log 2>&1" % (GMX, protfilename, cenfname, x,y,z)
    os.system(cmd)

    if options.verbose:
        sys.stdout.write("done.\n")


### calculate protonation states ###
if options.verbose:
    sys.stdout.write("calculating protonation states... ")
    sys.stdout.flush()

noopt = "--noopt"
if options.pqroptimize:
    noopt = ""
cmd = "%s %s --ff=charmm --ph-calc-method=propka --with-ph=%f --ffout=charmm %s %s > pdb2pqr.log 2>&1" % (PDB2PQR, noopt, options.pH, cenfname, pqrfname)
os.system(cmd)

# fix pdb2pqr output #
outf = open(prtnfname, "w")
outf.write(pqr2pdb(pqrfname,True))
outf.close()

if options.verbose:
    sys.stdout.write("done.\n")


### build initial gromacs topology ###
if options.verbose:
    sys.stdout.write("generating initial topology... ")
    sys.stdout.flush()

# run pdb2gmx #
gmxffdesc = gmxffname.split(".ff")[0]
cmd = "%s pdb2gmx -f %s -o %s -p %s -ff %s -water tip3p > gmx.log 2>&1" % (GMX, prtnfname, gmxpdbfnm, gmxtopfnm, gmxffdesc)
os.system(cmd)

# check pdb2gmx output #
if not os.path.exists(gmxpdbfnm):
    sys.stderr.write("ERROR: pdb2gmx failed to generate initial topology\n")
    sys.exit(1)

nanlines = []
handle = open(gmxpdbfnm, "r")
for line in handle:
    if line.find("nan") > -1:
        nanlines.append(line)
handle.close()

if len(nanlines) > 0:
    sys.stderr.write("ERROR: pdb2gmx encountered missing heavy atoms in the following residues:\n")
    for nanl in nanlines:
        sys.stderr.write("  %s" % nanl)
    sys.exit(1)

# convert gromacs topology to charmm-friendly format #
outf = open(chrminfnm, "w")
outf.write(pqr2pdb(gmxpdbfnm,False))
outf.close()

if options.verbose:
    sys.stdout.write("done.\n")


################################################################################
# CALCULATE ENERGY GRIDS IN CHARMM                                             #

grids = ["uncharged", "charged"]

# nonbonded term commands #
# nbond_cmds[0] -> for calculations in a vacuum
# nbond_cmds[1] -> for calculations using the implicit solvent model
softstr = ""
if options.softcore:
    softstr = "soft"
nbond_cmds = ["nbonds nbxmod 5 atom rdiel fswitch vatom vdistance vfswitch cutnb 14.0 ctofnb 12.0 ctonnb 10.0 %s eps 4.0 e14fac 1.0 wmin 1.5" % softstr,
              "nbonds nbxmod 5 atom cdiel switch  vatom vdistance vswitch  cutnb 14.0 ctofnb 12.0 ctonnb 10.0 %s eps 1.0 e14fac 1.0 wmin 1.5" % softstr]

charmmstr = """* GRID energy calculation using CHARMM
*

! the log file gets ginormous with default print level
! comment this if you need to debug the charmm run
! setting print level < 2 will likely break things
prnlev 2

!--- read topology and force field files --------------------------------------!

open unit 20 read card name "%s"
read rtf card unit 20
close unit 20
open unit 20 read card name "%s"
read rtf card unit 20 append
close unit 20

open unit 20 read card name "%s"
read param card unit 20
close unit 20
open unit 20 read card name "%s"
read param card unit 20 append
close unit 20

!--- done read topology and force field files ---------------------------------!


!--- generate protein system --------------------------------------------------!

! read sequence
open unit 21 read form name "%s"
read sequ pdb unit 21

generate prot setup

! read coordinates
rewind unit 21
read coord pdb unit 21
close unit 21

! repair topology as best we can
ic fill preserve
ic parameter
ic build
hbuild select hydrogen end

! set up the grid coordinates
set extraspace = %f   ! extra spacing to add around the protein in angstroms
%s                    ! (possibly) calculate grid coordinates from the protein
set nxmin = %s                      ! minimum X atom
set nymin = %s                      ! minimum Y atom
set nzmin = %s                      ! minimum Z atom
set nxmax = %s                      ! maximum X atom
set nymax = %s                      ! maximum Y atom
set nzmax = %s                      ! maximum Z atom
calc maxX = @nxmax + @extraspace    ! maximum X coordinate
calc maxY = @nymax + @extraspace    ! maximum Y coordinate
calc maxZ = @nzmax + @extraspace    ! maximum Z coordinate
calc minX = @nxmin - @extraspace    ! minimum X coordinate
calc minY = @nymin - @extraspace    ! minimum Y coordinate
calc minZ = @nzmin - @extraspace    ! minimum Z coordinate
calc cenX = (@minX + @maxX) / 2.0   ! center  X coordinate
calc cenY = (@minY + @maxY) / 2.0   ! center  Y coordinate
calc cenZ = (@minZ + @maxZ) / 2.0   ! center  Z coordinate

!--- done generate protein system ---------------------------------------------!


!--- calculate energy grid ----------------------------------------------------!

! set implicit solvent parameters or vacuum parameters
%s

! set nonbonded calculation parameters
%s

set theprobe = uncharged
label probeloop

if @theprobe eq uncharged then
  set probeid = XRC
endif

if @theprobe eq charged then
  set probeid = XRN
  delete atom select segid prob end
endif

! generate probe atom sequence
read sequ card
* probe
1 @probeid
generate prob setup

! set position of probe atom
scalar x set @minX select segid prob end
scalar y set @minY select segid prob end
scalar z set @minZ select segid prob end

! set probe radius for implicit solvent model
scalar wmain set 2.3 sele resn XRC end
scalar wmain set 2.3 sele resn XRN end
scalar wmain set 1.8 sele resn XRO end

! get this probe's radius
scalar wmain stat select segid prob end
set myRadius = ?save

! initialize energies
energy
set baseVDW = ?vdW
set baseELE = ?elec
set baseENE = ?ener

open unit 12 write card name %s.@theprobe.qmapout
write title unit 12
*# qmap grid energy difference calculations
*#
*# min x,y,z @minX @minY @minZ
*# max x,y,z @maxX @maxY @maxZ
*# cen x,y,z @cenX @cenY @cenZ
*# grid spacing %f
*# initial energies @baseVDW @baseELE @baseENE
*#
*# x,y,z = grid position
*# dist  = minimum distance to protein atoms in angstroms
*#
*# potential energy differences (delta-G)
*#   vdw    = van der waals
*#   elec   = electrostatic
*#   totpot = total potential energy
*# x y z dist vdw elec totpot
echu 12

set x = @minX
label ix
  set y = @minY
  label iy
    set z = @minZ
    label iz

    ! set new probe position
    scalar x set @x select segid prob end
    scalar y set @y select segid prob end
    scalar z set @z select segid prob end

    ! get distance of probe from protein
    coor mindist sele segid prot end sele segid prob end

    ! only calculate system energy if probe is outside of protein
    if ?mind gt @myRadius then
        ! calculate new system energy
        energy
        calc diffVDW = ?vdW  - @baseVDW
        calc diffELE = ?elec - @baseELE
        calc diffENE = ?ener - @baseENE

        ! print energy results to file
        echo @x @y @z ?mind @diffVDW @diffELE @diffENE
    endif

    ! set adaptive grid spacing based on probe distance from protein
    set ginc = %f
    if ?mind gt %f set ginc = %f
    if ?mind le %f set ginc = %f

    incr z by @ginc
    if z le @maxZ goto iz
  incr y by @ginc
  if y le @maxY goto iy
incr x by @ginc
if x le @maxX goto ix

if @theprobe eq uncharged then
  set theprobe = charged
  goto probeloop
endif

close unit 12
stop

"""

if options.verbose:
    sys.stdout.write("calculating interaction grids... ")
    sys.stdout.flush()

# set the grid coordinates strings #
gextarr = ["coord stat", "?xmin", "?ymin", "?zmin", "?xmax", "?ymax", "?zmax"]
if options.setgrid:
    gextarr = [""] + options.setgrid.split(",")


# set the nonbond parameter string #
nbparamstr = nbond_cmds[0]
impsolparm = ""
if options.usesolvent:
    nbparamstr = nbond_cmds[1]
    impsolparm = "stream \"%s\"\ngbsw nang 50" % (options.charmmradiusstr)
    # charmm documentation suggests 'gbsw sgamma 0.005 nang 50' for the #
    # self-consistent implicit-solvent model, but this can cause a      #
    # segfault. Setting sgamma 0.0 (the default) seems to work          #

mycharmmstr = charmmstr % (options.charmmprobertf,options.charmmprotrtf, options.charmmprobeprm,options.charmmprotprm, chrminfnm, options.margin, gextarr[0], gextarr[1], gextarr[2], gextarr[3], gextarr[4], gextarr[5], gextarr[6], impsolparm, nbparamstr, protbasename, options.gridsize, options.gridsize, options.adaptmar, options.gridsize*2.0, options.tightmar, options.gridsize/2.0)

outf = open(chrminpfnm, "w")
outf.write(mycharmmstr)
outf.close()

cmd = "%s < %s > charmm.log 2>&1" % (CHARMM, chrminpfnm)
os.system(cmd)

if options.verbose:
    sys.stdout.write("done.\n")

# convert charmm output file name,       #
# which currently imposes all-lower-case #
# to protbasename                        #
for qmapname in grids:
    os.system("mv %s.%s.qmapout %s.%s.qmapout" % (protbasename.lower(), qmapname, protbasename, qmapname))

# quietly remove gmx backup files #
os.system("rm -f \#*\#")

# remove log and intermediate files #
if not options.logs:
    os.system("rm -f %s" % gmxffname)
    os.system("rm -f pdb2pqr.log")
    os.system("rm -f %s" % pqrfname)
    os.system("rm -f %s.propka" % prodbase)
    os.system("rm -f gmx.log")
    os.system("rm -f gmx.center.log")
    os.system("rm -f %s" % cenfname)
    os.system("rm -f posre.itp")
    os.system("rm -f %s" % prtnfname)
    os.system("rm -f %s" % gmxpdbfnm)
    os.system("rm -f %s" % gmxtopfnm)
    os.system("rm -f charmm.log")
    os.system("rm -f %s" % chrminpfnm)

runtime = time.time() - starttime

if options.verbose:
    sys.stdout.write("finished in %s.\n" % datetime.timedelta(seconds=runtime))
