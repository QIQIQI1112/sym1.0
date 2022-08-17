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
sym = SymCases(rigid_type, np.array([1, 3]), 1, 3)
sym_dic = sym.case_list(sym_no, 0)

trail_len = len(sym_dic)
# trail_len = 0 means there is no possible combination of atoms
if trail_len != 0:
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

    for i in range(trail_len):
        os.chdir(my_cwd_cases + f'/sym{sym_no}')
        os.mkdir(f'wyckoff{i + 1}')
        os.chdir('..')
        # wyckoff_no is where the rigid body gonna be
        wyckoff_no = sym_dic[f'{i + 1}'][0]
        # get the coordinates of the wyckoff position
        wyckoff_tuple = sg[f'sg_{sym_no}'][f's{wyckoff_no}'][2]
        wyckoff_rigid = list(wyckoff_tuple)
        # get the symmetry of the wyckoff position
        wyckoff_sym = sg[f'sg_{sym_no}'][f's{wyckoff_no}'][1]
        # get the total number of atoms and setup the sim box
        single_rigid = rigid(sym_no, bond).atom_no
        rigid_muti = sg[f'sg_{sym_no}'][f's{wyckoff_no}'][0]
        total_rigid = single_rigid * rigid_muti
        total_single_atom = rigid_muti * ratio
        total_atom = total_rigid + total_single_atom

        # generate single atoms
        # total number of combination
        single_case = len(sym_dic[f'{i + 1}'][1])
        for n in range(single_case):
            os.chdir(my_cwd_cases + f'/sym{sym_no}' + f'/wyckoff{i + 1}')
            os.mkdir(f'case{n + 1}')
            os.chdir(my_cwd_cases + f'/sym{sym_no}' + f'/wyckoff{i + 1}' + f'/case{n + 1}')
            os.mkdir('trail')
            os.mkdir('relax')
        
        # symmetry information on this case
        sym_information(compound_rigid, sym_no, wyckoff_no, rigid_muti, wyckoff_sym, sym_dic[f'{i + 1}'][1][n])
        os.chdir(my_cwd_cases)

        # generate the initial system combine rigid and single atoms
        single_no = len(sym_dic[f'{i + 1}'][1][0])
        dis_transform_return = 'no'
        system_initial = np.zeros((1, 3))
        system_transform = np.zeros((1, 3))
        compound = [*compound_rigid, *(compound_single_atom * single_no)]
        density_adjust = density

        while dis_transform_return != 'yes':
            density_adjust = density_adjust - 0.002
            count = 0
            while dis_transform_return != 'yes':
                # sim box
                sim_box = initial_cell(box_type, total_atom, density_adjust)
                # print(f'sim_box:{sim_box}')
                # generate rigid body
                rigid_generate = rigid(sym_no, bond).atoms_position() @ rigid(sym_no, bond).rotaiton_matrix(wyckoff_sym) @ setup_random_rotation(sym_no, wyckoff_sym, 360.0)
                # print(f'position:{rigid(sym_no, bond).atoms_position()}')
                # print(f'rotation:{rigid(sym_no, bond).rotaiton_matrix(wyckoff_sym)}')
                # print(f'setup_rotation:{setup_random_rotation(sym_no, wyckoff_sym, 360.0)}')
                # print(f'rigid_generate:{rigid_generate}')
                # print('-----------------------')
                # move rigid to wyckoff position
                rigid_initial = move_rigid_to_wyckoff(rigid_generate, wyckoff_rigid, rigid_type, sim_box)
                # print(f'rigid_initial:{rigid_initial}')
                # generate single atoms
                wyckoff_case = []
                for j in range(single_no):
                    wyckoff_single = sg[f'sg_{sym_no}'][sym_dic[f'{i + 1}'][1][n][j]][2]
                    wyckoff_case.append(wyckoff_single)
                single_initial = generate_single_atoms(single_no, wyckoff_case, sim_box)

                system_initial = np.concatenate((rigid_initial, single_initial))
                system_initial_dir = direct_cartesian_transform(system_initial, sim_box, 'CtoD')
                
                # space group transform
                system_initial_tuple = []
                for k in range(single_rigid + single_no):
                    system_initial_tuple.append(tuple(system_initial_dir[k]))
                system_transform = crystal(compound, system_initial_tuple, spacegroup=sym_no, cell=sim_box)
                # atom distance check
                dis_transform_return = dis_check(system_transform.positions, sim_box, dis_limit)
                count += 1
                if count == 100:
                    break
            print(f'count:{count}')
            print(f'density_adjust:{density_adjust}')
            # print(f'sim_box:{sim_box}')
            # print(f'system_initial_dir:{system_initial_dir}')

        # setup all input file 
        # write poscar and atom position before transform
        os.chdir('..')
        write_poscar(system_transform.positions, compound, rigid_muti, ratio, sim_box)
        write_system_initial(system_initial)

        # copy file to individual trail
        setup_file_copy(my_cwd + '/SymTest', sym_no, i, n)

else:
    print('no case found')