#!/usr/bin/env python3

import sys
import os
import shutil
import glob
import math
from optparse import OptionParser, OptionGroup
from _version import __version__
from pqr2pdb  import pqr2pdb


### HELPER FUNCTIONS ------------------------------------------------------- ###

# helper function to calculate the volume #
# of a typical triclinic box              #
def triclinic_vol(x,y,z,a,b,g):
    bv    = x*y*z
    cosa  = math.cos(math.radians(a))
    cosb  = math.cos(math.radians(b))
    cosg  = math.cos(math.radians(g))
    cosa2 = cosa*cosa
    cosb2 = cosb*cosb
    cosg2 = cosg*cosg
    t2    = 2.0 * cosa * cosb * cosg
    srprt = math.sqrt(1.0-cosa2-cosb2-cosg2+t2)
    return bv * srprt


# helper function to convert .xvg time,x,y,z location format #
# to .csv time,x,y,z location format                         #
def xvgxyz_to_csvxyz(inxvgfn, outcsvfn):
    ## convert .xvg to time,x,y,z format in angstroms ##
    #
    # XVG format will be like:
    #
    # time1 X1 Y1 Z1 X2 Y2 Z2 X3 Y3 Z3 ...
    # time2 X1 Y1 Z1 X2 Y2 Z2 X3 Y3 Z3 ...
    # ...
    #
    # NOTE: these will be in nm; we should probably multiply by 10 to get angstroms.
    #
    handle = open(inxvgfn, "r")
    outf   = open(outcsvfn, "w")
    outf.write("ps,xAng,yAng,zAng\n")
    for line in handle:
        if len(line) < 1 or line[0] == "#" or line[0] == "@":
            continue
        linearr = line.split()
        time = float(linearr[0])
        xyz  = [float(x)*10.0 for x in linearr[1:]]
        for i in range(0,len(xyz),3):
            x = xyz[i]
            y = xyz[i+1]
            z = xyz[i+2]
            outf.write("%.2f,%.2f,%.2f,%.2f\n" % (time,x,y,z))
    handle.close()
    outf.close()
    return


# helper function to create selection format for #
# a given residue name and list of atoms         #
def get_sel_format(res,atoms):
    atomstr = "( " + " or ".join([ 'name "{}"'.format(x) for x in atoms]) + " )"
    return "( resname \"%s\" and %s )" % (res, atomstr)


# contact types for which we can identify #
# hydrogen bonds or apolar contacts       #
INTERACTION_CONTACT_TYPES = ("HBACC", "HBDON", "APOLAR")

# helper function to get/check interaction groups for     #
# hydrogen bonds or apolar contacts                       #
# calculates appropriate groups for h-bond donor/acceptor #
# interactions or apolar contact interactions             #
# returns group and parameters needed for gmx hbond       #
def get_interaction_group(inttype, fragmapids, cmdline_options):
    others = ""
    hbondoptions = "-r %f -a %f" % (cmdline_options.hbond_rad/10.0, cmdline_options.hbond_angle)

    if inttype == INTERACTION_CONTACT_TYPES[0]:
        if "HBDON" in fragmapids:
            others = "HBDON"
        elif "PHBDON" in fragmapids:
            others = "PHBDON"
        else:
            sys.stderr.write("WARNING: can't calculate interactions for type %s\n" % inttype)

    elif inttype == INTERACTION_CONTACT_TYPES[1]:
        if "HBACC" in fragmapids:
            others = "HBACC"
        elif "NHBACC" in fragmapids:
            others = "NHBACC"
        else:
            sys.stderr.write("WARNING: can't calculate interactions for type %s\n" % inttype)

    elif inttype == INTERACTION_CONTACT_TYPES[2]:
        if "APOLAR" in fragmapids:
            others = "APOLAR"
            hbondoptions = "-r %f -contact" % (cmdline_options.contact_rad/10.0)
        else:
            sys.stderr.write("WARNING: can't calculate interactions for type %s\n" % inttype)

    else:
        sys.stderr.write("WARNING: can't calculate interactions for type %s\n" % inttype)

    return (others, hbondoptions)

### END HELPER FUNCTIONS --------------------------------------------------- ###


### read relevent environment variables ###

# maybe they set the environment variable #
MDAPIDIR = os.environ.get("MDAPIDIR")

# check the user's path for executables #
GMX      = shutil.which("gmx")
PDB2PQR  = shutil.which("pdb2pqr")
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

# remove trailing "/" if needed #
if len(MDAPIDIR) > 0 and MDAPIDIR[-1] == "/":
    MDAPIDIR = MDAPIDIR[:-1]
if len(GMX) > 0 and GMX[-1] == "/":
    GMX = GMX[:-1]
if len(PDB2PQR) > 0 and PDB2PQR[-1] == "/":
    PDB2PQR = PDB2PQR[:-1]


### set up command-line argument parsing ###
optparser = OptionParser(usage="usage: %prog [options] setup|run|analyze PROTFILE.pdb",
                         version=__version__,
                         description="Molecular Dynamics Analysis of Protein Interactions")

optparser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                     help="print some runtime info to the screen [default: %default]")
optparser.add_option("-q", "--quiet", action="store_false", dest="verbose",
                     help="no news is (probably) good news")
optparser.add_option("--clean", action="store_true", dest="clean",
                     help="remove unneeded files from the current directory, keeping only the files needed to run subsequent analysis steps")
optparser.add_option("--wipe", action="store_true", dest="wipe",
                     help="completely wipe out the current directory, rather than running analysis.\nWARNING: this will delete all files in the directory, except the PROTFILE.pdb file!")

group = OptionGroup(optparser, "Data and Helper Options", "where to find data files and helper programs")
group.add_option("--forcefield", action="store", type="string", dest="forcefield",
                 help="set directory containing force field data to DIR [default: %default]",
                 metavar="DIR")
group.add_option("--mdapidir", action="store", type="string", dest="mdapidir",
                 help="set directory containing the MDAPI installation to DIR [default: %default] Note that if you set the MDAPI installation directory this way, you must specify the force field, fragmap and protein options explicitly using the --forcefield=DIR, --fragmap=FILE and --protein=FILE options; you cannot use the default options unless you set the MDAPIDIR environment variable",
                 metavar="DIR")
group.add_option("--pdb2pqr", action="store", type="string", dest="pdb2pqr",
                 help="set name of pdb2pqr execuable program to FILE [default: %default]",
                 metavar="FILE")
group.add_option("--gmx", action="store", type="string", dest="gromacs",
                 help="set name of gromacs gmx executable program to FILE [default: %default]",
                 metavar="FILE")
optparser.add_option_group(group)

# options for controlling how MD is set up #
group = OptionGroup(optparser, "Molecular Dynamics Options", "control the MD run")
group.add_option("--pH", action="store", type="float", dest="pH",
                 help="set the pH for protonation calculations to PH [default: %default]",
                 metavar="PH")
group.add_option("--pqroptimize", action="store_true", dest="pqroptimize",
                 help="use pdb2pqr to optimize the protein's hydrogen bond network [default: %default]")
group.add_option("--nopqroptimize", action="store_false", dest="pqroptimize",
                 help="don't allow pdb2pqr to optimize the hydrogen bond network; it will still perform debumping, though")
group.add_option("--neutral", action="store_true", dest="neutralize",
                 help="neutralize the system charge prior to MD simulation [default: %default]")
group.add_option("--noneutral", action="store_false", dest="neutralize",
                 help="do not neutralize the system charge before MD simulation")
group.add_option("--margin", action="store", type="float", dest="margin",
                 help="set the simulation box margin to MAR angstroms [default: %default]",
                 metavar="MAR")
group.add_option("--boxtype", action="store", type="choice", choices=("triclinic", "cubic", "dodecahedron", "octahedron"),
                 dest="boxtype", help="set simulation box to TYPE [default: %default]",
                 metavar="TYPE")
group.add_option("--threads", action="store", type="int", dest="mdthreads",
                 help="set number of MD threads to THRDS [default: %default]", metavar="THRDS")
group.add_option("--rseed", action="store", type="int", dest="randseed",
                 help="set the random number seed to SEED [default: generated randomly]",
                 metavar="SEED")
group.add_option("--continue", action="store_true", dest="continuemd",
                 help="continue previous production MD run [default: %default]")
optparser.add_option_group(group)

# options for fragment types and densities #
group = OptionGroup(optparser, "Fragment Options", "how much and what types of fragments to include")
group.add_option("--fragdensity", action="store", type="float", dest="fragdensity",
                 help="set fragment density to DEN [default: %default]", metavar="DEN")
group.add_option("--fragmap", action="store", type="string", dest="fragmap",
                 help="sets fragmap definition to FILE [default: %default] Note that all fragments must be available in the forcefield/mol/ directory as both .pdb and .itp files",
                 metavar="FILE")
optparser.add_option_group(group)

# options for results analysis #
group = OptionGroup(optparser, "Analysis Options", "what types of results analyses to do")
group.add_option("--center", action="store_true", dest="center",
                 help="translate protein to x=0,y=0,z=0 [default: %default]")
group.add_option("--nocenter", action="store_false", dest="center",
                 help="don't translate protein to x=0,y=0,z=0; leave it wherever it is")
group.add_option("--protein", action="store", type="string", dest="proteindef",
                 help="sets protein definition to FILE [default: %default]", metavar="FILE")
group.add_option("--burnin", action="store", type="int", dest="burnin",
                 help="remove the first NS nanoseconds from the production MD run [default: %default]",
                 metavar="NS")
group.add_option("--fragmaps", action="store_true", dest="calcfragmaps",
                 help="calculate fragmap and protein exclusion map coordinates [default: %default]")
group.add_option("--nofragmaps", action="store_false", dest="calcfragmaps",
                 help="don't calculate any fragmap or protein exclusion coordinates")
group.add_option("--proteinmaps", action="store_true", dest="calcproteinmaps",
                 help="calculate protein interaction map coordinates [default: %default]")
group.add_option("--noproteinmaps", action="store_false", dest="calcproteinmaps",
                 help="don't calculate any protein interaction map coordinates")
group.add_option("--residues", action="store_true", dest="calcresidues",
                 help="calculate per-residue scores for hydrogen bonding and hydrophobic interactions [default: %default]")
group.add_option("--noresidues", action="store_false", dest="calcresidues",
                 help="don't calculate any per-residue interaction scores")
group.add_option("--hbond_rad", action="store", type="float", dest="hbond_rad",
                 help="set hydrogen bond radius to NUM angstroms [default: %default]", metavar="NUM")
group.add_option("--hbond_angle", action="store", type="float", dest="hbond_angle",
                 help="set hydrogen bond angle to NUM [default: %default]", metavar="NUM")
group.add_option("--contact_rad", action="store", type="float", dest="contact_rad",
                 help="set hydrophobic contact radius to NUM angstroms [default: %default]", metavar="NUM")
optparser.add_option_group(group)

# set option defaults #
defff = "charmm36.ff"
deffm = "standard.fragmapdef"
defpr = "standard.proteindef"
if MDAPIDIR:
    defff = "%s/data/%s" % (MDAPIDIR,defff)
    deffm = "%s/data/%s" % (MDAPIDIR,deffm)
    defpr = "%s/data/%s" % (MDAPIDIR,defpr)

optparser.set_defaults(verbose=True,
                       clean=False,
                       wipe=False,
                       pH=7.0,
                       pqroptimize=True,
                       neutralize=True,
                       margin=15,               # 15 angstrom margin for sim box
                       boxtype="dodecahedron",  # rhombic dodecahedron sim box
                       mdthreads=8,
                       randseed=int(os.urandom(2).hex(), 16),
                       continuemd=False,
                       forcefield=defff,
                       mdapidir=MDAPIDIR,
                       pdb2pqr=PDB2PQR,
                       gromacs=GMX,
                       fragdensity=0.001138,    # ~1.9M; >= ~0.85M for each fragment
                       fragmap=deffm,
                       center=True,
                       proteindef=defpr,
                       burnin=10,
                       calcfragmaps=True,
                       calcproteinmaps=True,
                       calcresidues=True,
                       hbond_rad=3.5,
                       hbond_angle=30.0,
                       contact_rad=3.8)


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


if options.forcefield[-1] == "/":
    options.forcefield = options.forcefield[:-1]


### set the files we need in the end ###
gmxffname    = options.forcefield.split("/")[-1]
silcstopf    = "%s_mdapi.top" % protbasename
posrecaf     = "posre_protein_ca.itp"
equilbase    = "%s_mdapi.equil" % protbasename
equilndxf    = "%s.ndx" % equilbase
equilrecpdbf = "%s_mdapi.equilcen.pdb" % protbasename

prodbase     = "%s_mdapi.prodmd" % protbasename
prodtprf     = "%s.tpr" % prodbase
prodxtcf     = "%s.xtc" % prodbase
prodedrf     = "%s.edr" % prodbase
prodcptf     = "%s.cpt" % prodbase

pbcfitxtcf   = "%s.pbc.fit.xtc" % prodbase
fragmapndxf  = "%s.ndx" % prodbase
resintf      = "%s.resinteractions.csv" % prodbase
vispdbf      = "%s.central.pdb" % prodbase

setup_keep_files   = [protfilename, gmxffname, silcstopf, equilndxf, equilrecpdbf, posrecaf, "posre.itp"]
run_keep_files     = [prodtprf, prodxtcf, prodedrf, prodcptf]
analyze_keep_files = [pbcfitxtcf, fragmapndxf, resintf, vispdbf]
results_files_csv  = glob.glob("*.csv")   # dynamically identify .csv  files #
results_files_dx   = glob.glob("*.dx")    # dynamically identify .dx   files #
results_files_attr = glob.glob("*.attr")  # dynamically identify .attr files #
qmap_keep_files    = glob.glob("*_qmap.pdb") + glob.glob("*.qmapout")
keep_files         = setup_keep_files + run_keep_files + analyze_keep_files + results_files_csv + results_files_dx + results_files_attr + qmap_keep_files

# we'll try to save scripts, as well #
py_files      = glob.glob("*.py")
ba_files      = glob.glob("*.bash")
sh_files      = glob.glob("*.sh")
program_files = py_files + ba_files + sh_files

keep_files += program_files


### set other files we'll use ###
pqrf         = "%s.pqr" % protbasename
pdb2gmxf     = "%s_4pdb2gmx_ph.pdb" % protbasename
gmxpdbf      = "%s_gmx.pdb" % protbasename
gmxtopf      = "%s_gmx.top" % protbasename
gmxtprf      = "%s_gmx.tpr" % protbasename
gmxboxpdbf   = "%s_gmx_box.pdb" % protbasename
silcspdbf    = "%s_mdapi.pdb" % protbasename
silcstprf    = "%s_mdapi.tpr" % protbasename
eminbase     = "%s_mdapi.emin" % protbasename
emintprf     = "%s.tpr" % eminbase
eminpdbf     = "%s.pdb" % eminbase
equiltprf    = "%s.tpr" % equilbase
equilpdbf    = "%s.pdb" % equilbase


### check for --wipe option ###
if options.wipe:
    for f in glob.glob("*"):
        if f not in program_files:
            if os.path.isfile(f) and f != protfilename:
                os.remove(f)
            elif os.path.islink(f):
                os.remove(f)
            elif os.path.isdir(f):
                shutil.rmtree(f)
    sys.exit(0)


### check for --clean option ###
if options.clean:
    for f in glob.glob("*"):
        if f not in keep_files:
            if os.path.isfile(f) or os.path.islink(f):
                os.remove(f)
            elif os.path.isdir(f):
                shutil.rmtree(f)
    sys.exit(0)


### done --wipe or --clean; we'll need a setup|run|analyze option now ###
if len(args) < 2:
    optparser.error("incorrect number of arguments: setup|run|analyze and PROTFILE.pdb filename required")
runtype = args[0].lower()

valid_runtypes = ["setup", "run", "analyze"]

if runtype not in valid_runtypes:
    errstr = "incorrect runtype: %s must be" % runtype
    for vrt in valid_runtypes[:-1]:
        errstr += " %s" % vrt
    errstr += " or %s" % valid_runtypes[-1]
    optparser.error(errstr)

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
    sys.stderr.write("  installation, or provide the loaction of gmx using the --gmx=GMX option\n")

PDB2PQR = options.pdb2pqr
if not PDB2PQR and runtype == valid_runtypes[0]:
    sys.stderr.write("ERROR: pdb2pqr executable is not in your PATH\n")
    sys.stderr.write("  please either update your PATH environment variable to include your pdb2pqr\n")
    sys.stderr.write("  installation, or provide the location of pdb2pqr using the --pdb2pqr=PDB2PQR option\n")

if not MDAPIDIR or not GMX or (not PDB2PQR and runtype == valid_runtypes[0]):
    sys.exit(1)

if options.verbose:
    sys.stdout.write("executing %s!\n" % runtype)


### parse fragmap definition for setup or analysis ###
fragmap_molecule_name_to_label = {}
fragmap_selection_name_to_label_atoms = {}
fragmap_fragmap_name_to_selection_names = {}

if runtype == valid_runtypes[0] or runtype == valid_runtypes[2]:
    if not os.path.isfile(options.fragmap):
        sys.stderr.write("ERROR: fragmap definition file %s not found\n" % options.fragmap)
        sys.exit(1)

    handle = open(options.fragmap, "r")
    for line in handle:
        if len(line) == 0 or line[0] == "#":
            continue
        linearr = line.split()
        if len(linearr) < 3:
            continue
        deftype = linearr[1]
        if deftype == "MOLECULE":
            mname  = linearr[2]
            mlabel = linearr[3]
            fragmap_molecule_name_to_label[mname] = mlabel
        elif deftype == "SELECTION":
            sname  = linearr[2]
            mlabel = linearr[3]
            atoms  = linearr[4:]
            fragmap_selection_name_to_label_atoms[sname] = (mlabel,atoms)
        elif deftype == "FRAGMAP":
            fname  = linearr[2]
            snames = linearr[3:]
            fragmap_fragmap_name_to_selection_names[fname] = snames
        else:
            sys.stderr.write("ERROR: fragmap definition file %s is invalid\n" % options.fragmap)
            sys.stderr.write("        found definition for %s\n" % deftype)
            sys.exit(1)
    handle.close()


################################################################################
# SETUP                                                                        #
if runtype == valid_runtypes[0]:

    ### link force field files to local directory for gromacs ###
    if options.verbose:
        sys.stdout.write("linking force field files to local directory for gromacs... ")
        sys.stdout.flush()
    if os.path.islink(gmxffname):
        os.remove(gmxffname)
    elif os.path.isdir(gmxffname):
        shutil.rmtree(gmxffname)
    if not os.path.isdir(options.forcefield):
        sys.stdout.write("ERROR: force field %s not found\n" % (options.forcefield))
        sys.exit(1)
    os.symlink(options.forcefield, gmxffname)
    if options.verbose:
        sys.stdout.write("done.\n")


    ### calculate protonation states ###
    if options.verbose:
        sys.stdout.write("calculating protonation states... ")
        sys.stdout.flush()
    noopt = "--noopt"
    if options.pqroptimize:
        noopt = ""
    cmd = "%s %s --ff=charmm --ph-calc-method=propka --with-ph=%f --ffout=charmm %s %s > pdb2pqr.log 2>&1" % (PDB2PQR, noopt, options.pH, protfilename, pqrf)
    os.system(cmd)
    # fix pdb2pqr output #
    outf = open(pdb2gmxf, "w")
    outf.write(pqr2pdb(pqrf,True))
    outf.close()
    if options.verbose:
        sys.stdout.write("done.\n")


    ### build initial gromacs topology ###
    if options.verbose:
        sys.stdout.write("generating initial topology... ")
        sys.stdout.flush()

    # run pdb2gmx #
    gmxffdesc = gmxffname.split(".ff")[0]
    cmd = "%s pdb2gmx -f %s -o %s -p %s -ff %s -water tip3p -merge all > gmx.pdb2gmx.log 2>&1" % (GMX, pdb2gmxf, gmxpdbf, gmxtopf, gmxffdesc)
    os.system(cmd)

    # check pdb2gmx output #
    if not os.path.exists(gmxpdbf):
        sys.stderr.write("ERROR: pdb2gmx failed to generate initial topology\n")
        sys.exit(1)

    nanlines = []
    handle = open(gmxpdbf, "r")
    for line in handle:
        if line.find("nan") > -1:
            nanlines.append(line)
    handle.close()

    if len(nanlines) > 0:
        sys.stderr.write("ERROR: pdb2gmx encountered missing heavy atoms in the following residues:\n")
        for nanl in nanlines:
            sys.stderr.write("  %s" % nanl)
        sys.exit(1)

    if options.verbose:
        sys.stdout.write("done.\n")


    ### generate initial empty box ###
    if options.verbose:
        sys.stdout.write("generating empty system box... ")
        sys.stdout.flush()

    # make empty box #
    realmargin = float(options.margin) / 10.0  # convert angstroms to nm #
    cmd = "%s editconf -f %s -o %s -d %f -bt %s -c > gmx.editconf.log 2>&1" % (GMX, gmxpdbf, gmxboxpdbf, realmargin, options.boxtype)
    os.system(cmd)

    if options.verbose:
        sys.stdout.write("done.\n")


    ### generate positional restraints on all backbone alpha-Carbon atoms ###
    if options.verbose:
        sys.stdout.write("generating positional restraints for protein C-alpha atoms... ")
        sys.stdout.flush()

    # calculate TPR file for positional restraints #
    cmd = "%s grompp -f %s/data/simple.mdp -c %s -p %s -o %s > gmx.grompp.log 2>&1" % (GMX, MDAPIDIR, gmxboxpdbf, gmxtopf, gmxtprf)
    os.system(cmd)

    if not os.path.isfile(gmxtprf):
        sys.stderr.write("ERROR grompp failed to generate tpr file for positional restraints\n")
        sys.exit(1)

    # select alpha-Carbon atoms #
    outf = open("select.selection", "w")
    outf.write("(atomname CA or atomname \"C1'\")")
    outf.close()
    cmd = "%s select -f %s -s %s -on selection.ndx -sf select.selection > gmx.select.log 2>&1" % (GMX, gmxboxpdbf, gmxtprf)
    os.system(cmd)

    # generate positional restraints #
    cmd = "%s genrestr -f %s -fc 50.208 50.208 50.208 -o %s -n selection.ndx > gmx.genrestr.log 2>&1" % (GMX, gmxboxpdbf, posrecaf)
    os.system(cmd)

    # insert positional restraints into gromacs topology file #
    outf = open("tmp.top", "w")
    handle = open(gmxtopf, "r")
    for line in handle:
        if line == "; Include water topology\n":
            outf.write("; Include C-alpha positional restraints file\n")
            outf.write("#ifdef POSRE_CA\n")
            outf.write("#include \"%s\"\n" % posrecaf)
            outf.write("#endif\n\n")
        outf.write(line)
    handle.close()
    outf.close()
    shutil.move("tmp.top", gmxtopf)

    if options.verbose:
        sys.stdout.write("done.\n")


    ### insert SILCS fragment itps into the gromacs topology file ###
    if options.verbose:
        sys.stdout.write("adding fragment definitions to the topology... ")
        sys.stdout.flush()

    # set fragments #
    fragments = fragmap_molecule_name_to_label.keys()

    if len(fragments) == 0:
        sys.stderr.write("ERROR: fragmap definition from %s empty\n" % options.fragmap)
        sys.exit(1)

    fragerr = False
    for frag in fragments:
        if not os.path.isfile("%s/mol/%s.itp" % (gmxffname, frag)) or not os.path.isfile("%s/mol/%s.pdb" % (gmxffname, frag)):
            sys.stderr.write("ERROR: no definition for fragment %s in force field %s\n" % (frag, gmxffname))
            fragerr = True
    if fragerr:
        sys.exit(1)

    # insert fragment files into gromacs topology file #
    outf = open("tmp.top", "w")
    handle = open(gmxtopf, "r")
    for line in handle:
        if line == "; Include water topology\n":
            outf.write("; Include MDAPI fragments\n")
            for frag in fragments:
                outf.write("#include \"./%s/mol/%s.itp\"\n" % (gmxffname, frag))
            outf.write("\n")
        outf.write(line)
    handle.close()
    outf.close()
    shutil.move("tmp.top", gmxtopf)

    if options.verbose:
        sys.stdout.write("done.\n")


    ### solvate system box with fragments and waters ###
    if options.verbose:
        sys.stdout.write("generating solvated system... ")
        sys.stdout.flush()

    # calculate fragment information #
    nfragtypes = len(fragments)
    xtala = -1.0
    xtalb = -1.0
    xtalc = -1.0
    angla = -1.0
    anglb = -1.0
    anglg = -1.0

    handle = open(gmxboxpdbf, "r")
    for line in handle:
        if line.find("CRYST1") != -1:
            linearr = line.split()
            xtala = float(linearr[1])
            xtalb = float(linearr[2])
            xtalc = float(linearr[3])
            angla = float(linearr[4])
            anglb = float(linearr[5])
            anglg = float(linearr[6])
            break
    handle.close()

    nfrags_per_box = round((triclinic_vol(xtala,xtalb,xtalc,angla,anglb,anglg) * options.fragdensity) / float(nfragtypes))

    shutil.copy(gmxboxpdbf, silcspdbf)
    shutil.copy(gmxtopf, silcstopf)

    frags_done=1
    for frag in fragments:
        # insert fragment into pdb box #
        cmd = "%s insert-molecules -seed %d -f %s -ci %s/mol/%s.pdb -o temp.pdb -nmol %d > gmx.insert_%s.log 2>&1" % (GMX, options.randseed, silcspdbf, gmxffname, frag, nfrags_per_box, frag)
        os.system(cmd)

        # get fragment label #
        fraglabel = fragmap_molecule_name_to_label[frag]

        if not fraglabel:
            sys.stderr.write("ERROR: can't find label for fragment %s\n" % frag)
            sys.exit(1)

        # calculate real number of fragments inserted #
        real_frag_count = 0
        last_num = -1
        handle = open("temp.pdb", "r")
        for line in handle:
            linearr = line.split()
            if len(linearr) > 4 and linearr[0] == "ATOM" and linearr[3] == fraglabel:
                curr_num = int(linearr[4])
                if curr_num != last_num:
                    real_frag_count += 1
                    last_num = curr_num
        handle.close()

        # write fragment entry into the topology file #
        outf = open(silcstopf, "a")
        outf.write("%s               %d\n" % (fraglabel, real_frag_count))
        outf.close()

        # copy pdb file with fragments into silcs pdb file #
        shutil.move("temp.pdb", silcspdbf)

        if options.verbose:
            sys.stdout.write("\n  finished %s (%d/%d) %d fragments placed (%d attempted)" % (frag, frags_done, nfragtypes, real_frag_count, nfrags_per_box))
        frags_done += 1

    # insert waters #
    cmd = "%s solvate -cp %s -o temp.pdb -p %s > gmx.solvate.log 2>&1" % (GMX, silcspdbf, silcstopf)
    os.system(cmd)
    shutil.move("temp.pdb", silcspdbf)

    if options.verbose:
        sys.stdout.write("\n  finished SOL\n")

    if options.verbose:
        sys.stdout.write("done.\n")


    ### neutralize system charge ###
    if options.neutralize:
        if options.verbose:
            sys.stdout.write("neutralizing system charge... ")
            sys.stdout.flush()

        # calculate TPR file for genion #
        cmd = "%s grompp -f %s/data/simple.mdp -c %s -p %s -o %s > gmx.grompp.genion.log 2>&1" % (GMX, MDAPIDIR, silcspdbf, silcstopf, silcstprf)
        os.system(cmd)

        # replace water with NA/CL ions to neutralize charge #
        cmd = "echo SOL | %s genion -s %s -p %s -o tmp.pdb -pname NA -nname CL -neutral > gmx.genion.log 2>&1" % (GMX, silcstprf, silcstopf)
        os.system(cmd)

        # gmx deletes the CRYST1 box, so we'll need to put it back in! #
        cryst1 = ""
        handle = open(silcspdbf, "r")
        for line in handle:
            linearr = line.split()
            if len(linearr) > 0 and linearr[0] == "CRYST1":
                cryst1 = line.strip()
                break
        handle.close()

        handle  = open("tmp.pdb", "r")
        outfile = open(silcspdbf, "w")
        for line in handle:
            linearr = line.split()
            if len(linearr) > 0 and linearr[0] == "MODEL":
                outfile.write("%s\n" % cryst1)
            outfile.write(line)
        handle.close()
        outfile.close()
        os.remove("tmp.pdb")

        if options.verbose:
            sys.stdout.write("done.\n")


    ### energy-minimize SILCS box ###
    if options.verbose:
        sys.stdout.write("emergy-minimizing initial conformation... ")
        sys.stdout.flush()
    cmd = "%s grompp -f %s/data/emin.mdp -c %s -r %s -p %s -o %s -maxwarn 1 > gmx.grompp.emin.log 2>&1" % (GMX, MDAPIDIR, silcspdbf, silcspdbf, silcstopf, emintprf)
    os.system(cmd)

    if not os.path.isfile(emintprf):
        sys.stderr.write("ERROR: grompp failed to generate tpr file for energy minimization\n")
        sys.exit(1)

    cmd = "%s mdrun -nt %d -s %s -deffnm %s -c %s > gmx.mdrun.emin.log 2>&1" % (GMX, options.mdthreads, emintprf, eminbase, eminpdbf)
    os.system(cmd)
    if not os.path.isfile(eminpdbf):
        sys.stderr.write("ERROR: mdrun energy minimization failed\n")
        sys.exit(1)
    if options.verbose:
        sys.stdout.write("done.\n")


    ### equilibrate SILCS box ###
    if options.verbose:
        sys.stdout.write("equilibrating initial conformation... ")
        sys.stdout.flush()

    # create index file for solvent and solute #
    outf = open("mdapi_select.txt", "w")
    for frag in fragments:
        outf.write("r %s | " % frag)
    outf.write("r SOL | r NA | r CL\nq\n")
    outf.close()

    cmd = "%s make_ndx -f %s -o %s > gmx.make_ndx.equil.log < mdapi_select.txt 2>&1" % (GMX, eminpdbf, equilndxf)
    os.system(cmd)

    mdapi_index = sum([line.find("[") != -1 for line in open(equilndxf).readlines()])
    prot_lipid_index = mdapi_index;
    mdapi_index = mdapi_index-1;
    outf = open("mdapi_select_2.txt", "w")
    outf.write("name %d SOLVENT\n" % mdapi_index)
    outf.write("! %d\n" % mdapi_index)
    outf.write("name %d SOLUTE\n" % prot_lipid_index)
    outf.write("q\n")
    outf.close()

    cmd = "%s make_ndx -f %s -n %s -o %s > gmx.make_ndx.equil.log < mdapi_select_2.txt 2>&1" % (GMX, eminpdbf, equilndxf, equilndxf)
    os.system(cmd)

    if not os.path.isfile(equilndxf):
        sys.stderr.write("ERROR: make_ndx failed to generate index file for equilibration\n")
        sys.exit(1)

    # generate equilibration tpr file #
    cmd = "%s grompp -f %s/data/equil.mdp -n %s -c %s -r %s -p %s -o %s -maxwarn 1 > gmx.grompp.equil.log 2>&1" % (GMX, MDAPIDIR, equilndxf, eminpdbf, eminpdbf, silcstopf, equiltprf)
    os.system(cmd)

    if not os.path.isfile(equiltprf):
        sys.stderr.write("ERROR: grompp failed to generate tpr file for equilibration\n")
        sys.exit(1)

    # equilibrate system #
    cmd = "%s mdrun -nt %d -s %s -deffnm %s -c %s > gmx.mdrun.equil.log 2>&1" % (GMX, options.mdthreads, equiltprf, equilbase, equilpdbf)
    os.system(cmd)

    if not os.path.isfile(equilpdbf):
        sys.stderr.write("ERROR: mdrun failed to equilibrate system\n")
        sys.exit(1)

    if options.verbose:
        sys.stdout.write("done.\n")


    ### re-center equilibration pdb ###
    if options.verbose:
        sys.stdout.write("compacting equilibrated system... ")
        sys.stdout.flush()

    cmd = "echo System | %s trjconv -f %s -s %s -o %s -pbc mol -ur compact > gmx.trjconv.equilcen.log 2>&1" % (GMX, equilpdbf, equiltprf, equilrecpdbf)
    os.system(cmd)

    if not os.path.isfile(equilrecpdbf):
        sys.stderr.write("ERROR: trjconv failed to compact equilibrated system\n")
        sys.exit(1)

    if options.verbose:
        sys.stdout.write("done.\n")
        sys.stdout.write("finished %s; final system written to file %s for visualization.\n" % (valid_runtypes[0], equilrecpdbf))
# SETUP                                                                        #
################################################################################


################################################################################
# RUN                                                                          #
elif runtype == valid_runtypes[1]:

    # checkpoint interval (minutes) #
    cpinterval = 60

    ### continue existing production MD run ###
    if options.continuemd:
        if options.verbose:
            sys.stdout.write("continuing existing production MD run... ")
            sys.stdout.flush()

        # check that existing MD run files exist #
        missing_files = []
        for reqfile in run_keep_files:
            if not os.path.exists(reqfile):
                missing_files.append(reqfile)
        if len(missing_files) > 0:
            sys.stderr.write("ERROR: required run files not found:\n")
            for mf in missing_files:
                sys.stdout.write("  %s\n" % mf)
            sys.stdout.write("you can't continue a run that hasn't been started yet\n")
            sys.exit(1)

        # continue existing production MD run #
        cpifname = "%s.cpt" % prodbase

        warnings = 0
        if not os.path.isfile(cpifname):
            warnings = 1

        else:
            cmd = "%s check -f %s > gmxcheck.out 2>&1" % (GMX, cpifname)
            os.system(cmd)

            warnings = 0
            handle = open("gmxcheck.out", "r")
            for line in handle:
                if line.find("Fatal error") == 0:
                    warnings += 1
            handle.close()
            os.remove("gmxcheck.out")

        if warnings > 0:
            prevcpifname = "%s_prev.cpt" % prodbase

            if not os.path.isfile(prevcpifname):
                sys.stderr.write("ERROR: no valid checkpoint file; continuation of production MD run not available\n")
                sys.exit(1)

            cmd = "%s check -f %s > gmxcheck.out 2>&1" % (GMX, prevcpifname)
            os.system(cmd)

            warnings = 0
            handle = open("gmxcheck.out", "r")
            for line in handle:
                if line.find("Fatal error") == 0:
                    warnings += 1
            handle.close()
            os.remove("gmxcheck.out")

            if warnings > 0:
                sys.stderr.write("ERROR: no valid checkpoint file; continuation of production MD run not available\n")
                sys.exit(1)

            shutil.move(prevcpifname, cpifname)

        cmd = "%s mdrun -nt %d -s %s -cpi %s -deffnm %s -cpt %d >> gmx.mdrun.prodmd.log 2>&1" % (GMX, options.mdthreads, prodtprf, cpifname, prodbase, cpinterval)
        os.system(cmd)

    ### execute new production MD run ###
    else:
        if options.verbose:
            sys.stdout.write("executing new production MD run... ")
            sys.stdout.flush()

        # check that required files from setup exist #
        missing_files = []
        for reqfile in setup_keep_files:
            if not os.path.exists(reqfile):
                missing_files.append(reqfile)
        if len(missing_files) > 0:
            sys.stderr.write("ERROR: required setup files not found:\n")
            for mf in missing_files:
                sys.stderr.write("  %s\n" % mf)
            sys.stderr.write("make sure you have successfully setup the system\n")
            sys.stderr.write("before executing production MD\n")
            sys.exit(1)

        # start production MD run #
        cmd = "%s grompp -f %s/data/prod.mdp -n %s -c %s -r %s -p %s -o %s -maxwarn 1 > gmx.grompp.prodmd.log 2>&1" % (GMX, MDAPIDIR, equilndxf, equilrecpdbf, equilrecpdbf, silcstopf, prodtprf)
        os.system(cmd)

        if not os.path.isfile(prodtprf):
            sys.stderr.write("ERROR: grompp failed to generate production MD topology\n")
            sys.exit(1)

        cmd = "%s mdrun -nt %d -s %s -deffnm %s -cpt %d > gmx.mdrun.prodmd.log 2>&1" % (GMX, options.mdthreads, prodtprf, prodbase, cpinterval)
        os.system(cmd)

    if options.verbose:
        sys.stdout.write("done.\n")
        sys.stdout.write("finished %s.\n" % valid_runtypes[1])
# RUN                                                                          #
################################################################################


################################################################################
# ANALYZE                                                                      #
elif runtype == valid_runtypes[2]:

    # check that required files from run exist #
    missing_files = []
    for reqfile in run_keep_files:
        if not os.path.exists(reqfile):
            missing_files.append(reqfile)
    if len(missing_files) > 0:
        sys.stderr.write("ERROR: required production MD output files not found:\n")
        for mf in missing_files:
            sys.stderr.write("  %s\n" % mf)
        sys.stderr.write("make sure you have successfully completed production MD\n")
        sys.stderr.write("before analyzing the results\n")
        sys.exit(1)

    ### center MD trajectory and remove burnin ###
    if options.verbose:
        sys.stdout.write("processing MD trajectory... ")
        sys.stdout.flush()
    boxcencmd = ""
    if options.center:
        boxcencmd = "-boxcenter zero"
    cmd = "echo Protein System | %s trjconv -f %s -s %s -o %s -tu ns -b %d -t0 0 -pbc mol -ur compact -center %s > gmx.trjconv.pbc.log 2>&1" % (GMX, prodxtcf, prodtprf, pbcfitxtcf, options.burnin, boxcencmd)
    os.system(cmd)
    if not os.path.isfile(pbcfitxtcf):
        sys.stderr.write("ERROR: problem processing MD trajectory\n")
        sys.exit(1)
    if options.verbose:
        sys.stdout.write("done.\n")


    ### write 'central' protein from fit trajectory ###
    ### as pdb file for visualization               ###
    cmd = "echo Protein Protein | %s cluster -f %s -s %s -cl %s -skip 100 -cutoff 100.0 -nofit > gmx.cluster.log 2>&1" % (GMX, pbcfitxtcf, prodtprf, vispdbf)
    os.system(cmd)
    if options.verbose:
        sys.stdout.write("wrote central protein to file %s for visualization.\n" % vispdbf)

    ### create index file that has fragmap and residue indices ###
    if options.verbose:
        sys.stdout.write("generating fragmap, protein and residue indices... ")
        sys.stdout.flush()

    ## parse protein definition file ##
    protein_defn = {}
    handle = open(options.proteindef, "r")
    for line in handle:
        if len(line) > 0 and line[0] != "#":
            linearr = line.split()
            if len(linearr) > 2:
                type  = linearr[0].upper()
                res   = linearr[1].upper()
                atoms = linearr[2:]
                if type not in protein_defn.keys():
                    protein_defn[type] = {}
                protein_defn[type][res] = atoms
    handle.close()

    ## write fragmap, protein and residue selecion definitions ##
    fragmap_ids = []   # save fragmap IDs for later
    protein_ids = []   # PROTEIN_APOLAR, PROTEIN_HBACC and PROEIN_HBDON
    residue_ids = []   # save residue IDs for later

    # write fragmap selections #
    fragmapself = "fragmaps_residues.selectiontxt"
    outf = open(fragmapself, "w")
    for (fragmapname,selectionnames) in fragmap_fragmap_name_to_selection_names.items():
        var = fragmapname.upper()
        fragmap_ids.append(var)
        outf.write("%s = ( " % var)
        allres = []
        for selname in selectionnames:
            (res,atoms) = fragmap_selection_name_to_label_atoms[selname]
            allres.append(get_sel_format(res,atoms))
        outf.write(" or ".join(allres))
        outf.write(" ) ;\n")

    # write protein selections #
    for type in protein_defn.keys():
        protein_type_id = "PROTEIN_" + type
        protein_ids.append(protein_type_id)
        outf.write("%s = ( " % protein_type_id)
        allres = []
        for (res,atoms) in protein_defn[type].items():
            allres.append(get_sel_format(res,atoms))
        outf.write(" or ".join(allres))
        outf.write(" ) ;\n")

    # write residue selections #
    handle = open(silcstopf, "r")
    for line in handle:
        linearr = line.split()
        if len(linearr) > 4 and linearr[0] == ";" and linearr[1] == "residue":
            resnr = int(linearr[2])
            resnm = linearr[3].upper()
            id = "%s%d" % (resnm, resnr)
            for type in protein_defn.keys():
                if resnm in protein_defn[type].keys():
                    myid = "%s_%s" % (id, type)
                    residue_ids.append(myid)
                    outf.write("%s = ( resnr %d and " % (myid, resnr))
                    atomstr = "( " + " or ".join([ 'name "{}"'.format(x) for x in protein_defn[type][resnm]]) + " )"
                    outf.write(atomstr)
                    outf.write(" ) ;\n")
    handle.close()

    # write fragmaps #
    for fragmapname in fragmap_ids:
        outf.write("%s ;\n" % fragmapname)

    # write protein terms #
    for proteintypename in protein_ids:
        outf.write("%s ;\n" % proteintypename)

    # write residues #
    for residuename in residue_ids:
        outf.write("%s ;\n" % residuename)

    outf.close()

    ## create fragmap and residues index ##
    fmapndx = "fragmaps_residues"
    cmd = "%s select -f %s -s %s -on %s -sf %s > gmx.select.fmapndx.log 2>&1" % (GMX, pbcfitxtcf, prodtprf, fmapndx, fragmapself)
    os.system(cmd)
    fmapndx += ".ndx"

    ## combine fragmap and original index files ##
    cmd = "cat %s %s > %s" % (equilndxf, fmapndx, fragmapndxf)
    os.system(cmd)

    if options.verbose:
        sys.stdout.write("done.\n")


    ### calculate fragmaps ###
    if options.calcfragmaps:
        if options.verbose:
            sys.stdout.write("calculating fragmap coordinates...\n")

        donefm = 0
        totfm  = len(fragmap_ids) + 1
        # add default 'protein exclusion' map #
        protselname = "PROTEIN"
        for selname in fragmap_ids + [protselname]:
            fragmapname = selname.lower()

            ## get time,x,y,z loactions in .xvg format using gromacs ##
            selrpos = "dyn_mol_com"
            seltype = "dyn_mol_com"
            # protein needs to select the atoms, not the (sub)molecule #
            if selname == protselname:
                selrpos = "atom"
                seltype = "atom"
            fragmapxvgf = "%s_%s" % (protbasename, fragmapname)
            cmd = "%s trajectory -f %s -s %s -n %s -select %s -ox %s -selrpos %s -seltype %s -x -y -z > gmx.trajectory.%s.log 2>&1" % (GMX, pbcfitxtcf, prodtprf, fragmapndxf, selname, fragmapxvgf, selrpos, seltype, fragmapname)
            os.system(cmd)
            fragmapxvgf += ".xvg"

            outfn  = "%s.%s.csv" % (prodbase, fragmapname)
            xvgxyz_to_csvxyz(fragmapxvgf, outfn)

            donefm += 1
            if options.verbose:
                sys.stdout.write("  done %s (%d/%d)\n" % (fragmapname, donefm, totfm))

        if options.verbose:
            sys.stdout.write("done.\n")
    ### done calculating fragmaps ###


    ### calculate protein interaction maps ###
    if options.calcproteinmaps:
        if options.verbose:
            sys.stdout.write("calculating protein interaction group coordinates...\n")

        doneprot = 0
        totprot  = len(protein_ids)
        for selname in protein_ids:
            protmapname = selname.lower()

            # sort out protein interaction type #
            prot_int_type = selname.split("_")[-1]
            if prot_int_type not in INTERACTION_CONTACT_TYPES:
                sys.stderr.write("WARNING: can't calculate protein-fragment interactions for protein type %s\n" % selname)
                continue

            # sort out fragment type #
            (others, hbondoptions) = get_interaction_group(prot_int_type, fragmap_ids, options)
            if others == "":
                continue

            # calculate atoms forming hydrogen bonds #
            # or apolar contacts                     #
            pintndxf = protmapname + ".ndx"
            cmd = "echo %s %s | %s hbond -f %s -s %s -n %s -hbn %s %s > gmx.hbond.%s.log 2>&1" % (selname, others, GMX, pbcfitxtcf, prodtprf, fragmapndxf, pintndxf, hbondoptions, protmapname)
            os.system(cmd)

            # find protein interacting atoms from the hbond index file #
            if prot_int_type == INTERACTION_CONTACT_TYPES[0]:
                ndxlabel = "hbonds_%s-%s" % (selname,others)
                ndxindex = 2
            elif prot_int_type == INTERACTION_CONTACT_TYPES[1]:
                ndxlabel = "hbonds_%s-%s" % (selname,others)
                ndxindex = 0
            elif prot_int_type == INTERACTION_CONTACT_TYPES[2]:
                ndxlabel = "contacts_%s-%s" % (selname,others)
                ndxindex = 1

            protatoms = set([])
            handle = open(pintndxf, "r")
            foundlabel = False
            line = handle.readline()
            while line and not foundlabel:
                if len(line) > 0 and line[0] == "[":
                    label = line.split()[1]
                    if label == ndxlabel:
                        foundlabel = True
                line = handle.readline()
            while line and len(line) > 0 and line[0] != "[":
                linearr = line.strip().split()
                protatoms.add(linearr[ndxindex])
                line = handle.readline()
            handle.close()

            protatomslist = list(protatoms)
            protatomslist.sort()

            # add protein interacting atoms to the index file #
            pint_ndx_label = selname + "_FRAGINTERACTING"
            handle = open(fragmapndxf, "a")
            handle.write("[ %s ]\n" % pint_ndx_label)
            for patom in protatomslist:
                handle.write("\t%s" % patom)
            handle.write("\n")
            handle.close()

            ## get time,x,y,z loactions in .xvg format using gromacs ##
            selrpos = "atom"
            seltype = "atom"
            pmapxvgf = "%s_%s" % (protbasename, protmapname)
            cmd = "%s trajectory -f %s -s %s -n %s -select %s -ox %s -selrpos %s -seltype %s -x -y -z > gmx.trajectory.%s.log 2>&1" % (GMX, pbcfitxtcf, prodtprf, fragmapndxf, pint_ndx_label, pmapxvgf, selrpos, seltype, protmapname)
            os.system(cmd)
            pmapxvgf += ".xvg"

            outfn  = "%s.%s.csv" % (prodbase, protmapname)
            xvgxyz_to_csvxyz(pmapxvgf, outfn)

            doneprot += 1
            if options.verbose:
                sys.stdout.write("  done %s (%d/%d)\n" % (protmapname,doneprot,totprot))

        if options.verbose:
            sys.stdout.write("done.\n")
    ### done calculating protein interaction maps ###


    ### calculate residue interactions ###
    if options.calcresidues:
        if options.verbose:
            sys.stdout.write("calculating residue-fragment interactions...\n")

        interaction_proportions = {}
        residues = []
        for resid in residue_ids:
            res = resid.split("_")[0]
            if res not in residues:
                residues.append(res)
        for residue in residues:
            if residue not in interaction_proportions.keys():
                interaction_proportions[residue] = {}
                for inttype in INTERACTION_CONTACT_TYPES:
                    interaction_proportions[residue][inttype] = (0.0,0.0)

        doneres = 0
        totres  = len(residue_ids)
        for resid in residue_ids:
            residue = resid.split("_")[0]
            inttype = resid.split("_")[1]

            (others, hbondoptions) = get_interaction_group(inttype, fragmap_ids, options)
            if others == "":
                continue

            grpintf  = "interactions.%s" % resid
            cmd = "echo %s %s | %s hbond -f %s -s %s -n %s -num %s %s > gmx.hbond.%s.log 2>&1" % (resid, others, GMX, pbcfitxtcf, prodtprf, fragmapndxf, grpintf, hbondoptions, resid)
            os.system(cmd)
            grpintf += ".xvg"

            ## parse .xvg file to get proportion of time with interactions ##
            contacts   = 0
            sum        = 0
            timeframes = 0
            handle = open(grpintf, "r")
            for line in handle:
                if len(line) > 0 and line[0] != "#" and line[0] != "@":
                    linearr = line.split()
                    if len(linearr) > 1:
                        hits = int(linearr[1])
                        sum += hits
                        if hits > 0:
                            contacts += 1
                        timeframes += 1
            handle.close()

            interaction_proportions[residue][inttype] = (float(contacts)/float(timeframes), float(sum)/float(timeframes))

            doneres += 1
            if options.verbose:
                sys.stdout.write("  done %s (%d/%d)\n" % (resid,doneres,totres))

        outf = open(resintf, "w")
        outf.write("residue")
        for inttype in INTERACTION_CONTACT_TYPES:
            outf.write(",Ptime_%s,Ave_%s" % (inttype,inttype))
        outf.write("\n")
        for residue in residues:
            outf.write(residue)
            for inttype in INTERACTION_CONTACT_TYPES:
                (c,s) = interaction_proportions[residue][inttype]
                outf.write(",%f,%f" % (c,s))
            outf.write("\n")
        outf.close()

        if options.verbose:
            sys.stdout.write("done.\n")
    ### done calculating residue interactions ###


    if options.verbose:
        sys.stdout.write("finished %s.\n" % valid_runtypes[2])
# ANALYZE                                                                      #
################################################################################


# quietly clean up gromacs backups #
for f in glob.glob("#*#"):
    os.remove(f)
