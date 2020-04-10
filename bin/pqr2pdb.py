#!/usr/bin/env python3

charmm_names_to_fix = {"SPP":"ASPP", "LUP":"GLUP"}
disulfide_atoms     = ["1CB", "1SG"]
disulfide_residues  = ["ISU"]

def parseResidue(the_residue, four_letter_codes):
    result = ""

    good_resi_name = the_residue[0][4] + " " # add spacer for 'normal' 3-letter names #
    for the_line in the_residue:
        the_name = the_line[4]
        if the_name in charmm_names_to_fix.keys():
            if four_letter_codes:
                good_resi_name = charmm_names_to_fix[the_name]
            else:
                good_resi_name = charmm_names_to_fix[the_name][0:3] + " "
            break

    for (atom,atom_num,atom_nam,element,resi_nam,resi_num,chain_id,x,y,z,occupancy,tempfactor) in the_residue:
        # sort out disulfides #
        if resi_nam.strip() in disulfide_residues and atom_nam.strip() in disulfide_atoms:
            atom_nam = "  " + atom_nam.strip()[1:] + " "
            element  = "           " + atom_nam.strip()[0]
        # write pdb line #
        result += atom
        result += atom_num
        result += atom_nam
        result += " "     # blank alternate location
        result += good_resi_name
        result += chain_id
        result += resi_num
        result += "    "  # blank insertion code, plus spacer columns
        result += x
        result += y
        result += z
        result += occupancy
        result += tempfactor
        result += element
        result += "\n"
    return result


def pqr2pdb(pqrfname, four_letter_codes):
    result = ""

    last_res_id = 0
    full_residue = []

    handle = open(pqrfname, "r")
    for line in handle:
        atom = line[0:6]

        # we'll parse only the ATOM entries #
        if atom.strip() != "ATOM":
            continue

        atom_num = line[6:11]
        atom_nam = line[11:16]
        element  = "           " + atom_nam.strip()[0]  # hack to get element symbol
        resi_nam = line[17:20]
        resi_num = line[22:26]
        chain_id = line[21]
        # hack to fix missing chain ID #
        if chain_id == " ":
            chain_id = "A"

        # coordinates (angstroms) #
        x = line[30:38]
        y = line[38:46]
        z = line[46:54]

        # 'bogus' occupancy and b-factor #
        occupancy  = "  1.00"
        tempfactor = "  0.00"

        # store residue information #
        curr_res_id = int(resi_num.strip())
        if curr_res_id != last_res_id:
            if len(full_residue) > 0:
                result += parseResidue(full_residue, four_letter_codes)
            full_residue = []
            last_res_id  = curr_res_id

        full_residue.append((atom,atom_num,atom_nam,element,resi_nam,resi_num,chain_id,x,y,z,occupancy,tempfactor))

    handle.close()
    result += "TER\nEND\n"
    return result

if __name__ == "__main__":
    import sys
    from optparse import OptionParser

    optparser = OptionParser(usage="usage: %prog [options] INFILE.pqr",
                             version="%prog v0.5",
                             description="fixes a pqr file built using pdb2pqr with --ffout=charmm option. repairs residue names ASPP and GLUP and converts to a pretty strict pdb format")

    optparser.add_option("-4", "--four-letter-names", action="store_true", dest="fourlettercodes",
                         help="use 4-letter residue names for ASPP and GLUP (for silcs, etc) [default: %default]")
    optparser.add_option("-3", "--three-letter-names", action="store_false", dest="fourlettercodes",
                         help="only use 3-letter residue names (for charmm)")

    optparser.set_defaults(fourlettercodes=True)

    (options, args) = optparser.parse_args()

    if len(args) < 1:
        optparser.error("incorrect number of arguments: INFILE.pqr filename required")
    protfilename = args[-1]

    sys.stdout.write(pqr2pdb(protfilename, options.fourlettercodes))
    sys.exit(0)
