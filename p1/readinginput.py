from functions import *

def read_file(file, key):
    with open(file, 'r') as file_read:
        for line in file_read:
            if key in line:
                value = line.split()
                value_len = len(value)
                value_return = value[2:value_len]
    return value_return

rigid_type_input = read_file('INPUT', 'rigid_type')[0]
# rigid_type = rigid_input_mapping(rigid_type_input)
rigid_type = read_file('INPUT', 'rigid_type')[0]
compound_rigid = read_file('INPUT', 'compound_rigid')
compound_single_atom = read_file('INPUT', 'compound_single_atom')
ratio = int(read_file('INPUT', 'ratio')[0])
sym_no = int(read_file('INPUT', 'sym_no')[0])
rigid_no = int(read_file('INPUT', 'rigid_no')[0])
density = float(read_file('INPUT', 'density')[0])
box_type = read_file('INPUT', 'box_type')[0]
bond = float(read_file('INPUT', 'bond')[0])
dis_limit = float(read_file('INPUT', 'dis_limit')[0])
star = float(read_file('INPUT', 'step_ini')[0])
end = float(read_file('INPUT', 'step_fin')[0])
temp_star = float(read_file('INPUT', 'temp_ini')[0])
step_update = read_file('INPUT', 'step_algorithm')[0]
temp_update = read_file('INPUT', 'temp_algorithm')[0]
step_number = int(read_file('INPUT', 'step_number')[0])
rigid_type_input = read_file('INPUT', 'rigid_type')[0]
single_sym_restricted = int(read_file('INPUT', 'single_sym_restricted')[0])
internal_circulation = float(read_file('INPUT', 'internal_circulation')[0])
debug = int(read_file('INPUT', 'debug')[0])

cwd = os.getcwd()
sym_path = cwd + '/sym_information'
isfile1 = os.path.isfile(sym_path)
isfile = str(isfile1)
if isfile == 'True':
    wyckoff_rigid = int(read_file('sym_information', 'WyckoffRigid')[0])
    rigid_muti = int(read_file('sym_information', 'RigidMuti')[0])
    rigid_atom = int(read_file('sym_information', 'RigidAtom')[0])
    single_Atom = len(read_file('sym_information', 'SingleAtomSymmetry'))
    single_atom_symmetry = read_file('sym_information', 'SingleAtomSymmetry')
else:
    pass
