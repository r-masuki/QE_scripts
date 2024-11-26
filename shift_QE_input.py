import numpy as np
from pymatgen.core import Structure

from spglib import get_symmetry_dataset, get_spacegroup_type
from spglib import get_magnetic_symmetry_dataset, get_magnetic_spacegroup_type  # develop branch
from spinspg import get_spin_symmetry

from ase import Atoms
from ase.io import read
import re
import sys

from parse_qe_input import parse_namelist, parse_qe_cards
from write_qe_input import write_qe_input
    
              
# parse input file of quantum espresso
args = sys.argv
filename = args[1]

qe_namelists = parse_namelist(filename, ["control", "system", "electrons"])
qe_cards = parse_qe_cards(filename)

# print(qe_cards)

shift_fractional = [-0.25, -0.125, -0.5]

for iat in range(len(qe_cards["atomic_positions"])):
    for ixyz in range(3):
        qe_cards["atomic_positions"][iat]["position"][ixyz] += shift_fractional[ixyz]

write_qe_input(qe_namelists, qe_cards, "qe.in")

