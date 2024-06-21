import numpy as np
from pymatgen.core import Structure

from spglib import get_symmetry_dataset, get_spacegroup_type
from spglib import get_magnetic_symmetry_dataset, get_magnetic_spacegroup_type  # develop branch

from ase import Atoms
from ase.io import read
import re


def strip_comment(line):
    return line.split("!")[0].strip()

def parse_namelist(file_path, namelist_names):
    namelists = {}
    current_namelist = None

    with open(file_path, 'r') as file:
        lines = file.readlines()

    for line_tmp in lines:
        line_tmp = strip_comment(line_tmp.strip())
        for line in [x.strip() for x in re.split(',', line_tmp) if x != ""]:
            # line = line.strip()
            
            # Check if the line starts a new namelist
            if line.startswith('&'):
                current_namelist = line[1:]
                if current_namelist in namelist_names:
                    namelists[current_namelist] = {}
                else:
                    current_namelist = None
            # Check if the line ends the current namelist
            elif line == '/' and current_namelist:
                current_namelist = None
            # Parse key-value pairs within the current namelist
            elif current_namelist:
                key_value_pair = re.split(r'\s*=\s*', line.rstrip(','))
                if len(key_value_pair) == 2:
                    key, value = key_value_pair
                    # Remove quotes from strings and convert numeric values
                    value = value.strip()
                    if value.startswith("'") and value.endswith("'"):
                        value = value[1:-1]
                    elif value.replace('.', '', 1).isdigit() or \
                        (value.startswith('-') and value[1:].replace('.', '', 1).isdigit()):
                        value = float(value) if '.' in value else int(value)
                    namelists[current_namelist][key.strip()] = value

    return namelists

def parse_qe_cards(file_path):
    atomic_species = []
    atomic_positions = []
    k_points = {}

    with open(file_path, 'r') as file:
        lines = file.readlines()

    current_card = None
    for line in lines:
        line = line.strip()

        # Detect the start of a new card
        atomic_positions_pattern = re.compile(r'^ATOMIC_POSITIONS.*')
        if line == 'ATOMIC_SPECIES':
            current_card = 'ATOMIC_SPECIES'
            continue
        elif line.startswith('ATOMIC_POSITIONS'):
            current_card = 'ATOMIC_POSITIONS'
            continue
        elif line.startswith('K_POINTS'):
            current_card = 'K_POINTS'
            k_points['type'] = line.split()[1]
            continue
        elif line.startswith('&') or line == '/':
            current_card = None

        # Parse the content of the current card
        if current_card == 'ATOMIC_SPECIES':
            parts = line.split()
            if len(parts) == 3:
                atomic_species.append(parts)
            else:
                raise ValueError(f'Invalid line in ATOMIC_SPECIES card: {line}')
        elif current_card == 'ATOMIC_POSITIONS':
            parts = line.split()
            if len(parts) >= 4:
                element = parts[0]
                position = list(map(float, parts[1:4]))
                atomic_positions.append({
                    'element': element,
                    'position': position
                })
                # We don't read magnetization from ATOMIC_POSITIONS card in the current implementation
                # atomic_positions.append({
                #     'element': element,
                #     'position': position,
                #     'magnetic_moment': list(map(float, parts[4:])) if len(parts) > 4 else None
                # })
        elif current_card == 'K_POINTS':
            if 'points' not in k_points:
                k_points['points'] = []
            k_points['points'].append(list(map(int, line.split())))

    return {
        'atomic_species': atomic_species,
        'atomic_positions': atomic_positions,
        'k_points': k_points
    }


def parse_noncolmag_table(filename, qe_ntyp):


    with open(filename, "r") as f:
        lines = f.readlines()

    magmom_table = np.zeros(qe_ntyp)
    angle1_table = np.zeros(qe_ntyp)
    angle2_table = np.zeros(qe_ntyp)

    pattern_mags = re.compile(rf'starting_magnetization\((\d+)\)\s*=\s*([+-]?[0-9]*\.?[0-9]+)')
    pattern_angle1 = re.compile(rf'angle1\((\d+)\)\s*=\s*([+-]?[0-9]*\.?[0-9]+)')
    pattern_angle2 = re.compile(rf'angle2\((\d+)\)\s*=\s*([+-]?[0-9]*\.?[0-9]+)')

    for line in lines:
        match = pattern_mags.match(line.strip())
        if match:
            index = int(match.group(1))
            value = float(match.group(2))
            magmom_table[index-1] = value
        match = pattern_angle1.match(line.strip())
        if match:
            index = int(match.group(1))
            value = float(match.group(2))
            angle1_table[index-1] = value
        match = pattern_angle2.match(line.strip())
        if match:
            index = int(match.group(1))
            value = float(match.group(2))
            angle2_table[index-1] = value
    
    return magmom_table, angle1_table, angle2_table

# parse input file of quantum espresso
filename = "pw.Mn3Sn.scf.in"
qe_namelists = parse_namelist(filename, ["control", "system", "electrons"])
qe_cards = parse_qe_cards(filename)

print(qe_namelists["system"].keys())
print(qe_cards["atomic_species"])
print(qe_cards["atomic_positions"])
print(qe_cards["k_points"])


qe_ntyp = qe_namelists["system"]["ntyp"]
qe_nat = qe_namelists["system"]["nat"]
print("ntype = ", qe_ntyp)
print("nat = ", qe_nat)

# read the structure from the espresso input using ase.io
atoms = read(filename, format="espresso-in")

print(atoms.get_cell())
print(atoms.get_scaled_positions())
print(atoms.get_atomic_numbers())
print(atoms.get_initial_magnetic_moments())
# print(atoms.get_magnetic_moments())

# magnetization
magmom_table, angle1_table, angle2_table = parse_noncolmag_table(filename, qe_ntyp)
print(magmom_table)
print(angle1_table)
print(angle2_table)

