import numpy as np
from pymatgen.core import Structure

from spglib import get_symmetry_dataset, get_spacegroup_type
from spglib import get_magnetic_symmetry_dataset, get_magnetic_spacegroup_type  # develop branch
from spinspg import get_spin_symmetry

from ase import Atoms
from ase.io import read
import re
import sys

from parse_qe_input import parse_namelist, parse_qe_cards, get_noncol_magmoms

# parse input file of quantum espresso
args = sys.argv
filename = args[1]

qe_namelists = parse_namelist(filename, ["control", "system", "electrons"])
qe_cards = parse_qe_cards(filename)


atoms = read(filename, format="espresso-in")

# get symmetry operations
lattice = atoms.get_cell()
positions = atoms.get_scaled_positions()
numbers = atoms.get_atomic_numbers()
magmoms = get_noncol_magmoms(qe_namelists, qe_cards)
cell = (lattice, positions, numbers, magmoms)

Amat = np.array(lattice)
Ainv = np.linalg.inv(Amat)

dataset_as_noncollinear = get_magnetic_symmetry_dataset(cell, symprec=1e-2, is_axial=True)

uni_number = dataset_as_noncollinear["uni_number"]  # Range from 1 to 1651
magnetic_spacegroup_type = get_magnetic_spacegroup_type(uni_number)
nsym = len(dataset_as_noncollinear["rotations"])

# print space group information
# print("magnetic spacegroup type : ", magnetic_spacegroup_type)

pattern_tmp = re.compile(r'[\[\]]')

print("number of symmetry operations : \n", nsym)

for isym in range(nsym):
    print("Symmetry operation : ", isym)
    print("Rotation in real space (fractional):")
    print(re.sub(pattern_tmp, " ", 
                 str(dataset_as_noncollinear["rotations"][isym])))
    
    print("Translation :")
    print(re.sub(pattern_tmp, " ", 
                 str(dataset_as_noncollinear["translations"][isym]))
    )

    print("Rotation in spin space (Cartesian) :")
    print(re.sub(pattern_tmp, " ", 
                 str(np.dot(Amat.transpose(), 
                 np.dot(dataset_as_noncollinear["rotations"][isym], 
                 Ainv.transpose())))
                 )
    )

    print("Time reversal : \n", str(dataset_as_noncollinear["time_reversals"][isym]))


# Find spin symmetry operations
sog, rotations, translations, spin_rotations = get_spin_symmetry(
    lattice, positions, numbers, magmoms
)
nsym = len(rotations)

"""print(f"Spin-only group: {sog}")

# Some operations have nontrivial spin rotations
for isym in range(nsym):
    print(f"Spin symmetry operation : {isym}")
    print("Rotation :")
    print(rotations[isym])
    print("Translation :")
    print(translations[isym])
    print("Spin rotation :")
    print(spin_rotations[isym])"""
