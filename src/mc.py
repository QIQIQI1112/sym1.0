from ase.spacegroup import crystal
from readinginput import *
from randommove import *
from symdata import *

root_path = os.getcwd()
trail_path = root_path +'/trail'

# output file
output_path = root_path +'/output'
isfile1 = os.path.isdir(output_path)
isfile = str(isfile1)
if isfile == 'False':
    os.mkdir('output')

# log
log = WriteLog(output_path)
script_count = log.write_script_count()
big_loop = log.write_big_loop(script_count, internal_circulation)

# read POSCAR
poscar = ReadPOSCAR('POSCAR', trail_path, root_path)
cell_poscar = poscar.lattice

# get rigid and single atom array
total_atom = rigid_atom + single_Atom
system_initial = read_system_initial(trail_path, total_atom)
rigid_array = system_initial[0:rigid_atom, :]
single_array = system_initial[rigid_atom:total_atom, :]

# generate new system
single_atom_symmetry_len = len(single_atom_symmetry)
wyckoff_position_single = []
for i in range(single_atom_symmetry_len):
    wyckoff_position_single.append(sg[f'sg_{sym_no}'][single_atom_symmetry[i]][2])
compound = [*compound_rigid, *(compound_single_atom * single_Atom)]

# wyckoff position
wyckoff_position_rigid = sg[f'sg_{sym_no}'][f's{wyckoff_rigid}'][2]
wyckoff_element_rigid = sg[f'sg_{sym_no}'][f's{wyckoff_rigid}'][1]

# current temp
print(f'temp_update:{temp_update}')
temp = temp_algorithm(temp_update, temp_star, big_loop)

# generate a rigid class for center method
rigid_type_class = rigid_select(rigid_type)(1.0, 1.0)

close_dis_check = 'no'
print('looking for new atom positions...')
while close_dis_check != 'yes':
    rigid_after_move1 = AtomRandomMove(rigid_array).symmetry_restricted('rigid', step_update, temp, wyckoff_position_rigid, rigid_type)
    # image select
    rigid_after_move = rigid_after_move1.copy()
    for i in range(1, rigid_atom):
        image_select = AtomImage(rigid_after_move1[i], cell_poscar).close_image_position(rigid_after_move1[0])
        rigid_after_move[i] = image_select
    
    # rotation
    rigid_center = rigid_type_class.center(rigid_after_move)
    rigid_after_rotation = rigid_rotation(rigid_after_move, rigid_center).random_rotation(sym_no, wyckoff_element_rigid, 15)
    # single atom move
    single_after_move = AtomRandomMove(single_array).symmetry_restricted('single', step_update, temp, wyckoff_position_single)

    system_after_move = np.concatenate((rigid_after_rotation, single_after_move))
    # CtoD
    system_after_move_dir = direct_cartesian_transform(system_after_move, cell_poscar, 'CtoD')
    system_after_move_tuple = []
    for k in range(total_atom):
        system_after_move_tuple.append(tuple(system_after_move_dir[k]))
    # transform
    system_transform = crystal(compound, system_after_move_tuple, spacegroup=sym_no, cell=cell_poscar)
    close_dis_check = dis_check(system_transform.positions, cell_poscar, dis_limit)

    # write poscar and system_initial
    os.chdir(trail_path)
    write_poscar(system_transform.positions, compound, rigid_muti, ratio, cell_poscar)
    write_system_initial(system_after_move)
    os.chdir(root_path)
