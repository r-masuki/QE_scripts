import numpy as np
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor

from spglib import get_symmetry_dataset, get_spacegroup_type
from spglib import get_magnetic_symmetry_dataset, get_magnetic_spacegroup_type  # develop branch
from spinspg import get_spin_symmetry

from ase import Atoms
from ase.io import read
import sys

# parse input file of quantum espresso
args = sys.argv
input_filename = args[1]
output_filename = args[2]

# read the structure from the espresso input using ase.io
atoms = read(input_filename, format="espresso-in")
structure = AseAtomsAdaptor.get_structure(atoms)
structure.to(filename=output_filename)
