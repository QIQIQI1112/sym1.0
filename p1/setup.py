# initial_cell function only has cube and hex, universal vol 
"""
Class SymCases has the following directory format
The example is with tetrahedron and symmetry number 15)
{'1': [5, [['s1', 's6'], ['s2', 's6'], ['s3', 's6'], ['s4', 's6'], ['s5', 's6']]],
'2': [5, [['s1', 's2', 's3'], ['s1', 's2', 's4'], ['s1', 's2', 's5'], ['s1', 's3', 's4'], ['s1', 's3', 's5'],
['s1', 's4', 's5'], ['s2', 's3', 's4'], ['s2', 's3', 's5'], ['s2', 's4', 's5'], ['s3', 's4', 's5']]]}

'1' is the case index means the 1st case is with the rigid body at 5th Wyckoff position
['s1', 's6'] means single atoms are at s1 and s6

Add new rigid body: Add new rigid body information in rigid.py and rigid_generate function in functions.py
"""

from ase.spacegroup import crystal
from readinginput import *
from rigid import *
from casecount import *

rigid = rigid_select(rigid_type)
wyckoff_positions = ['x', 'y', 'z']
# trail_len = 0 means there is no possible combination of atoms
my_cwd = os.getcwd()
# make the test file where all the cases will be put in
path = my_cwd + '/SymTest'
isdir1 = os.path.isdir(path)
isdir = str(isdir1)
if isdir == 'False':
    os.mkdir('SymTest')

os.chdir(my_cwd + '/SymTest')
path_case = my_cwd + '/SymTest/' + f'sym{sym_no}'
isdir1 = os.path.isdir(path_case)
isdir = str(isdir1)
if isdir == 'False':
    os.mkdir(f'sym{sym_no}')
else:
    shutil.rmtree(f'{os.getcwd()}/sym{sym_no}', ignore_errors=True)
    # os.rmdir(f'sym{sym_no}')
    os.mkdir(f'sym{sym_no}')
my_cwd_cases = os.getcwd()
single_rigid = 5
total_rigid = single_rigid * rigid_no
total_single_atom = rigid_no * ratio
total_atom = total_rigid + total_single_atom

# generate single atoms
# total number of combination

os.chdir(my_cwd_cases + f'/sym{sym_no}')
os.mkdir('trail')
os.mkdir('relax')

# symmetry information on this case
# sym_information(compound_rigid, sym_no, wyckoff_no_list, rigid_muti, wyckoff_sym, sym_dic[f'{i + 1}'][1][n])
os.chdir(my_cwd_cases)

# generate the initial system combine rigid and single atoms
single_no = rigid_no * ratio
dis_transform_return = 'no'
system_initial = np.zeros((1, 3))
system_transform = np.zeros((1, 3))

density_adjust = density

while dis_transform_return != 'yes':
    density_adjust = density_adjust - 0.002
    count = 0
    while dis_transform_return != 'yes':
        # sim box
        sim_box = initial_cell(box_type, total_atom, density_adjust)
        # print(f'sim_box:{sim_box}')
        # generate rigid body
        rigid_generate_inter = rigid(sym_no, bond).atoms_position() @ rigid(sym_no, bond).rotaiton_matrix(['1']) @ setup_random_rotation(sym_no, ['1'], 360.0, sim_box)
        rigid_generate = move_rigid_to_wyckoff(rigid_generate_inter, wyckoff_positions, rigid_type, sim_box)
        if rigid_no > 1:
            for i in range(rigid_no - 1):
                rigid_generate_2_inter = rigid(sym_no, bond).atoms_position() @ rigid(sym_no, bond).rotaiton_matrix(['1']) @ setup_random_rotation(sym_no, ['1'], 360.0, sim_box)
                rigid_generate_2 = move_rigid_to_wyckoff(rigid_generate_2_inter, wyckoff_positions, rigid_type, sim_box)
                # print(rigid_generate)
                # print('--------------')
                # print(rigid_generate_2)
                rigid_generate = np.concatenate((rigid_generate, rigid_generate_2))

        # generate single atoms
        single_initial = generate_single_atoms_p1(single_no, sim_box)

        system_initial = np.concatenate((rigid_generate, single_initial))

        # atom distance check
        dis_transform_return = dis_check(system_initial, sim_box, dis_limit)
        # print(system_initial)
        # print('-----------------')
        count += 1
        if count == 100:
            break
    print(f'count:{count}')
    print(f'density_adjust:{density_adjust}')
    # print(f'sim_box:{sim_box}')
    # print(f'system_initial_dir:{system_initial_dir}')

# put rigid atom in order in POSCAR
system_final = rigid_poscar_transform_PS4(system_initial, rigid_no, 'RtoP')

# setup all input file 
# write poscar and atom position before transform
os.chdir(my_cwd)
write_poscar_p1(system_final, rigid_no, sim_box)

# copy file to individual trail
setup_file_copy_p1(my_cwd + '/SymTest', sym_no)
