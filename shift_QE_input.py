import numpy as np
from pymatgen.core import Structure

from spglib import get_symmetry_dataset, get_spacegroup_type
from spglib import get_magnetic_symmetry_dataset, get_magnetic_spacegroup_type  # develop branch
from spinspg import get_spin_symmetry

from ase import Atoms
from ase.io import read
import re
import sys


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
                current_namelist = current_namelist.lower()
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
                        # value = value[1:-1]
                        pass
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
            else:
                raise ValueError(f'Invalid line in ATOMIC_POSITIONS card: {line}')
        elif current_card == 'K_POINTS':
            if 'points' not in k_points:
                k_points['points'] = []
            k_points['points'].append(list(map(int, line.split())))

    return {
        'atomic_species': atomic_species,
        'atomic_positions': atomic_positions,
        'k_points': k_points
    }

def construct_ityp_table(qe_namelists, qe_cards):

    # check ntyp and nat
    qe_ntyp = qe_namelists["system"]["ntyp"]
    if(qe_ntyp != len(qe_cards["atomic_species"])):
        raise ValueError("ntyp is not consistent with ATOMIC_SPECIES card")
    qe_nat = qe_namelists["system"]["nat"]
    if(qe_nat != len(qe_cards["atomic_positions"])): 
        raise ValueError("nat is not consistent with ATOMIC_POSITIONS card")
    
    # construct ityp table
    ityp_table = np.zeros(qe_nat, dtype=int)
    for i in range(qe_nat):
        element = qe_cards["atomic_positions"][i]["element"]
        for j in range(qe_ntyp):
            if element == qe_cards["atomic_species"][j][0]:
                ityp_table[i] = j
                break
    return ityp_table
    

def get_noncol_magmoms(qe_namelists, qe_cards):

    # get magnetic moments for each atomic species
    mag_table = np.zeros(qe_ntyp)
    angle1_table = np.zeros(qe_ntyp)
    angle2_table = np.zeros(qe_ntyp)

    pattern_mags = re.compile(rf'starting_magnetization\((\d+)\)')
    pattern_angle1 = re.compile(rf'angle1\((\d+)\)')
    pattern_angle2 = re.compile(rf'angle2\((\d+)\)')

    for key, val in qe_namelists["system"].items():
        match = pattern_mags.match(key)
        if match:
            index = int(match.group(1))
            mag_table[index-1] = float(val)
        match = pattern_angle1.match(key)
        if match:
            index = int(match.group(1))
            angle1_table[index-1] = float(val)
        match = pattern_angle2.match(key)
        if match:
            index = int(match.group(1))
            angle2_table[index-1] = float(val)

    # construct magmom_table
    magmom_table = np.zeros([qe_nat, 3])
    for i in range(qe_ntyp):
        mag = mag_table[i]
        theta = angle1_table[i]/180.0*np.pi
        phi = angle2_table[i]/180.0*np.pi
        magmom_table[i,0] = mag*np.sin(theta)*np.cos(phi)
        magmom_table[i,1] = mag*np.sin(theta)*np.sin(phi)
        magmom_table[i,2] = mag*np.cos(theta)

    # construct ityp table
    itype_table = construct_ityp_table(qe_namelists, qe_cards)

    # get magnetic moments for each atom
    magmoms = np.zeros([qe_nat, 3])
    for i in range(qe_nat):
        magmoms[i] = magmom_table[itype_table[i]]

    return magmoms

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

