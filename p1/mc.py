from ase.spacegroup import crystal
from readinginput import *
from randommove import *
from symdata import *

root_path = os.getcwd()
trail_path = root_path +'/trail'

wyckoff_single = [['x', 'y', 'z']]
wyckoff_position_rigid = ['x', 'y', 'z']
wyckoff_element_rigid = ['1']
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

# positions array
# print(poscar.positions)
# print('================')
position_rigid = rigid_poscar_transform_PS4(poscar.positions, rigid_no, 'PtoR')

# current temp
print(f'temp_update:{temp_update}')
temp = temp_algorithm(temp_update, temp_star, big_loop)

# generate a rigid class for center method
rigid_type_class = rigid_select(rigid_type)(1.0, 1.0)

close_dis_check = 'no'
print('looking for new atom positions...')
print(f'big loop : {big_loop}')
while close_dis_check != 'yes':
    for n in range(rigid_no):
        rigid_array = position_rigid[n*5: (n+1)*5]
        rigid_after_move1 = AtomRandomMove(rigid_array, cell_poscar).symmetry_restricted('rigid', step_update, temp, wyckoff_position_rigid, rigid_type)
        # image select
        rigid_after_move = rigid_after_move1.copy()
        for i in range(1, 5):
            image_select = AtomImage(rigid_after_move1[i], cell_poscar).close_image_position(rigid_after_move1[0])
            rigid_after_move[i] = image_select
        
        # rotation
        rigid_center = rigid_type_class.center(rigid_after_move)
        rigid_after_rotation = rigid_rotation(rigid_after_move, rigid_center).random_rotation(sym_no, wyckoff_element_rigid, 350, cell_poscar)

        if n == 0:
            rigid_array_final = rigid_after_rotation
        else:
            rigid_array_final = np.append(rigid_array_final, rigid_after_rotation, axis = 0)
    # single atom move
    # print(rigid_array_final)
    # a = dis_check(rigid_array_final, cell_poscar, dis_limit)
    # print(a)
    # if a == 'no':
    #     print(rigid_array_final)
    # print('--------------')
    single_array = position_rigid[rigid_no*5:]
    # print(single_array)
    # print('--------------')
    wyckoff_position_single = wyckoff_single*len(single_array)
    single_after_move = AtomRandomMove(single_array, cell_poscar).symmetry_restricted('single', step_update, temp, wyckoff_position_single)

    system_after_move = np.concatenate((rigid_array_final, single_after_move))
    system_after_move_inorder = rigid_poscar_transform_PS4(system_after_move, rigid_no, 'RtoP')
    close_dis_check = dis_check(system_after_move_inorder, cell_poscar, dis_limit)
    # print(system_after_move_inorder)
    # print('-------------------------')
# CtoD

# write poscar and system_initial
os.chdir(trail_path)
# print(system_after_move_inorder)
# print('------------------------')
write_poscar_p1(system_after_move_inorder, rigid_no, cell_poscar)
os.chdir(root_path)
