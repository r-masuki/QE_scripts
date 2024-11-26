import numpy as np
from pymatgen.core import Structure

from spglib import get_symmetry_dataset, get_spacegroup_type
from spglib import get_magnetic_symmetry_dataset, get_magnetic_spacegroup_type  # develop branch
from spinspg import get_spin_symmetry

from ase import Atoms
from ase.io import read
import re
import sys

def write_qe_input(qe_namelists, qe_cards, filename):
    
    with open(filename, "w") as f:
        # write control namelist
        for namelist in qe_namelists.keys():
            f.write("&"+namelist+"\n")
            for key, val in qe_namelists[namelist].items():
                f.write("  "+key+" = "+str(val)+"\n")
            f.write("/\n")
        # write atomic species
        f.write("ATOMIC_SPECIES\n")
        for species in qe_cards["atomic_species"]:
            f.write("  "+species[0]+" "+species[1]+" "+species[2]+"\n")
        # write atomic positions
        f.write("ATOMIC_POSITIONS {crystal}\n")
        for atom in qe_cards["atomic_positions"]:
            f.write("  "+atom["element"]+" ")
            for coord in atom["position"]:
                f.write(str(coord)+" ")
            f.write("\n")
        # write k-points
        f.write("K_POINTS  "+qe_cards["k_points"]["type"]+"\n")
        for point in qe_cards["k_points"]["points"]:
            f.write("  ")
            for coord in point:
                f.write(str(coord)+" ")
            f.write("\n")

