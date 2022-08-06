import math
import os
import random
import shutil
from tkinter import Image
import numpy as np
from collections import Counter


def step_algorithm(algorithm, temp):
    shift_current = 0
    if algorithm == 'normal':
        s = np.random.normal(0, temp, 1)
        shift_current_1 = float(s[0])
        shift_current = shift_current_1
    return shift_current


def temp_algorithm(algorithm, temp_ini, script_count):
    current_temp = 0
    if algorithm == 'normal':
        current_temp = temp_ini / math.log(script_count)
    if algorithm == 'fast':
        current_temp = temp_ini * (0.99 ** script_count)
    return current_temp


def initial_cell(cell_type, atom_amount, density):
    if cell_type == 'cube':
            cube_len = round((atom_amount / density) ** (1./3), 10)
            cell = np.zeros((3, 3))
            for i in range(3):
                cell[i][i] = cube_len
    elif cell_type == 'hex':
        hex_len = round(((atom_amount / density) ** (1 / 3.0)) / math.sqrt(3.0), 10)
        cell = np.zeros((3, 3))
        cell[0][0], cell[2][2] = hex_len, hex_len
        cell[1][0] = - (hex_len / 2.0)
        cell[1][1] = (math.sqrt(3.0) / 2.0) * hex_len
    return cell


def direct_cartesian_transform(atom_positions, cell, transform_direction):
    atom_positions_final = atom_positions.copy()
    a1, a2, a3 = cell[0], cell[1], cell[2]
    # Build the matrix of lattice vectors stored column-wise
    # and get its inverse
    a_array = np.vstack([a1, a2, a3]).T
    a_array_inv = np.linalg.inv(a_array)
    if transform_direction == 'CtoD':
        atom_positions_final = np.matmul(a_array_inv, atom_positions.T).T
    elif transform_direction == 'DtoC':
        atom_positions_final = np.matmul(a_array, atom_positions.T).T
    return atom_positions_final


# atom positions are car coordinates
def move_atoms_into_box(atom_positions, cell):
    atom_amount = len(atom_positions)
    atom_positions_dir = direct_cartesian_transform(atom_positions, cell, 'CtoD')
    for i in range(atom_amount):
        for j in range(3):
            if atom_positions_dir[i][j] >= 1.0:
                while atom_positions_dir[i][j] >= 1.0:
                    atom_positions_dir[i][j] -= 1.0
            elif atom_positions_dir[i][j] < 0.0:
                while atom_positions_dir[i][j] < 0.0:
                    atom_positions_dir[i][j] += 1.0
    atom_positions_car = direct_cartesian_transform(atom_positions_dir, cell, 'DtoC')
    return atom_positions_car


def rotation_matrix(axis, angle):
    # make sure axis is unit vector
    if axis[0]**2 + axis[1]**2 + axis[2]**2 != 1.0:
        x = math.sqrt(axis[0]**2 + axis[1]**2 + axis[2]**2)
        axis_after_x = round(axis[0] / x, 6)
        axis_after_y = round(axis[1] / x, 6)
        axis_after_z = round(axis[2] / x, 6)
        axis_unit = [axis_after_x, axis_after_y, axis_after_z]
    else:
        axis_unit = axis
    
    y = float(angle / 360) * (2.0 * math.pi)
    r11 = math.cos(y) + (axis_unit[0]**2) * (1 - math.cos(y))
    r12 = axis_unit[0] * axis_unit[1] * (1 - math.cos(y)) - (axis_unit[2] * math.sin(y))
    r13 = axis_unit[0] * axis_unit[2] * (1 - math.cos(y)) + axis_unit[1] * math.sin(y)
    r21 = axis_unit[1] * axis_unit[0] * (1 - math.cos(y)) + axis_unit[2] * math.sin(y)
    r22 = math.cos(y) + axis_unit[1]**2 * (1 - math.cos(y))
    r23 = axis_unit[1] * axis_unit[2] * (1 - math.cos(y)) - axis_unit[0] * math.sin(y)
    r31 = axis_unit[2] * axis_unit[0] * (1 - math.cos(y)) - axis_unit[1] * math.sin(y)
    r32 = axis_unit[2] * axis_unit[1] * (1 - math.cos(y)) + axis_unit[0] * math.sin(y)
    r33 = math.cos(y) + axis_unit[2]**2 * (1 - math.cos(y))

    r_matrix = np.array([[r11, r12, r13],
                         [r21, r22, r23],
                         [r31, r32, r33]])
    return r_matrix


def sym_information(compound_rigid, sym_no, wyckoff_no, rigid_muti, wyckoff_sym, single_case):
    with open('sym_information', 'w') as sym_info:
        atom_in_rigid = len(compound_rigid)
        sym_info.write('RigidAtom : ' + f'{atom_in_rigid}' + '\n')
        sym_info.write('SymmetryNumber : ' + f'{sym_no}' + '\n')
        sym_info.write('WyckoffRigid : ' + f'{wyckoff_no}' + '\n')
        sym_info.write('RigidMuti : ' + f'{rigid_muti}' + '\n')
        wyckoff_sym_len = len(wyckoff_sym)
        sym_info.write('WyckoffSymmetry :')
        for i in range(wyckoff_sym_len):
            sym_info.write(' ' + f'{wyckoff_sym[i]}')
        sym_info.write('\n')
        single_case_len = len(single_case)
        sym_info.write('SingleAtomSymmetry :')
        for i in range(single_case_len):
            sym_info.write(' ' + f'{single_case[i]}')
        sym_info.write('\n')


# wyckoff_positions is dirc coordinates
def generate_single_atoms(atom_amount, wyckoff_positions, cell):
    atom_positions = np.zeros((atom_amount, 3))
    for i in range(atom_amount):
        for j in range(3):
            if wyckoff_positions[i][j] != 'x' and wyckoff_positions[i][j] != 'y' and wyckoff_positions[i][j] != 'z':
                atom_positions[i][j] = wyckoff_positions[i][j]
            else:
                atom_positions[i][j] = round(random.uniform(0, 1), 3)
    atom_positions = direct_cartesian_transform(atom_positions, cell, 'DtoC')
    return atom_positions


class ReadPOSCAR:
    def __init__(self, file, poscar_directory, root_directory):
        os.chdir(poscar_directory)
        with open(file, 'r') as ps:
            line = ps.readlines()
            lattice = np.zeros((3, 3))
            for i in range(3):
                for j in range(3):
                    lattice[i][j] = float(line[i + 2].split()[j])
            atom = line[6].split()
            atom_type = len(atom)
            atom_number = 0
            for i in range(atom_type):
                atom_number = atom_number + int(atom[i])
            atom_coordinate = np.zeros((atom_number, 3))
            for i in range(atom_number):
                for j in range(3):
                    atom_coordinate[i][j] = float(line[i + 8].split()[j])
        os.chdir(root_directory)
        self.lattice = lattice
        self.atom_number = atom_number
        self.positions = atom_coordinate


def read_system_initial(directory, total_atom):
    os.chdir(directory)
    with open('system_initial', 'r') as sys:
        line = sys.readlines()
        atom_coordinate = np.zeros((total_atom, 3))
        for i in range(total_atom):
            for j in range(3):
                atom_coordinate[i][j] = float(line[i].split()[j])
    return atom_coordinate


"""
use system_initial and poscar before_relax to find atom index in outcar after_relax
use atom number in rigid to find the right rigid image
"""
def sturcture_update(atom_in_rigid, system_initial, before_relax, after_relax, cell_before, cell_after):
    system_initial_inbox = move_atoms_into_box(system_initial, cell_before)
    before_relax_inbox = move_atoms_into_box(before_relax, cell_before)
    system_initial_atom_amount = len(system_initial)
    before_relax_atom_amount = len(before_relax)
    system_initial_after = np.zeros((system_initial_atom_amount, 3))
    right_image = np.zeros((system_initial_atom_amount, 3))
    atom_index = []
    for i in range(system_initial_atom_amount):
        for j in range(before_relax_atom_amount):
            atom_image = AtomImage(before_relax_inbox[j], cell_before)
            if atom_image.close_image_dis(system_initial_inbox[i]) < 0.001:
                atom_index.append(j)
    for i in range(system_initial_atom_amount):
        system_initial_after[i] = after_relax[atom_index[i]]
    system_initial_after_car = direct_cartesian_transform(system_initial_after, cell_after, 'DtoC')
    # find the right rigid image
    right_image[0] = system_initial_after_car[0]
    for i in range(1, atom_in_rigid):
        image = AtomImage(system_initial_after_car[i], cell_after)
        right_image[i] = image.close_image_position(right_image[0])
    for i in range(atom_in_rigid, system_initial_atom_amount):
        right_image[i] = system_initial_after_car[i]
    return right_image


def distance(atom1, atom2):
    x_dis = atom1[0] - atom2[0]
    y_dis = atom1[1] - atom2[1]
    z_dis = atom1[2] - atom2[2]
    dis = math.sqrt(x_dis ** 2 + y_dis ** 2 + z_dis ** 2)
    return dis


# atom_with_image is car coordinates
class AtomImage:
    def __init__(self, atom_with_image, cell):
        self.atom_with_image = atom_with_image
        self.cell = cell

    def image(self):    
        atom_with_image_dir = direct_cartesian_transform(self.atom_with_image, self.cell, 'CtoD')
        image = np.zeros((27, 3))
        for i in range(27):
            image[i] = atom_with_image_dir
        # images
        # change only 1 cord
        image[1][0] = atom_with_image_dir[0] + 1.0
        image[2][0] = atom_with_image_dir[0] - 1.0
        image[3][1] = atom_with_image_dir[1] + 1.0
        image[4][1] = atom_with_image_dir[1] - 1.0
        image[5][2] = atom_with_image_dir[2] + 1.0
        image[6][2] = atom_with_image_dir[2] - 1.0

        # change 2 cord (+x, +y) (+x, +z) (+y, +z)
        image[7][0] = atom_with_image_dir[0] + 1.0
        image[7][1] = atom_with_image_dir[1] + 1.0
        image[8][0] = atom_with_image_dir[0] + 1.0
        image[8][2] = atom_with_image_dir[2] + 1.0
        image[9][1] = atom_with_image_dir[1] + 1.0
        image[9][2] = atom_with_image_dir[2] + 1.0

        # change 2 cord (-x, -y) (-x, -z) (-y, -z)
        image[10][0] = atom_with_image_dir[0] - 1.0
        image[10][1] = atom_with_image_dir[1] - 1.0
        image[11][0] = atom_with_image_dir[0] - 1.0
        image[11][2] = atom_with_image_dir[2] - 1.0
        image[12][1] = atom_with_image_dir[1] - 1.0
        image[12][2] = atom_with_image_dir[2] - 1.0

        # change 2 cord (+x, -y) (+x, -z) (+y, -z)
        image[13][0] = atom_with_image_dir[0] + 1.0
        image[13][1] = atom_with_image_dir[1] - 1.0
        image[14][0] = atom_with_image_dir[0] + 1.0
        image[14][2] = atom_with_image_dir[2] - 1.0
        image[15][1] = atom_with_image_dir[1] + 1.0
        image[15][2] = atom_with_image_dir[2] - 1.0

        # change 2 cord (-x, +y) (-x, +z) (-y, +z)
        image[16][0] = atom_with_image_dir[0] - 1.0
        image[16][1] = atom_with_image_dir[1] + 1.0
        image[17][0] = atom_with_image_dir[0] - 1.0
        image[17][2] = atom_with_image_dir[2] + 1.0
        image[18][1] = atom_with_image_dir[1] - 1.0
        image[18][2] = atom_with_image_dir[2] + 1.0

        # change 3 cord (+x, +y, +z) (-x, -y, -z) (-x, +y, +z) (+x, -y, +z)
        #               (+x, +y, -z) (+x, -y, -z) (-x, +y, -z) (-x, -y, +z)
        image[19][0] = atom_with_image_dir[0] + 1.0
        image[19][1] = atom_with_image_dir[1] + 1.0
        image[19][2] = atom_with_image_dir[2] + 1.0

        image[20][0] = atom_with_image_dir[0] - 1.0
        image[20][1] = atom_with_image_dir[1] - 1.0
        image[20][2] = atom_with_image_dir[2] - 1.0

        image[21][0] = atom_with_image_dir[0] - 1.0
        image[21][1] = atom_with_image_dir[1] + 1.0
        image[21][2] = atom_with_image_dir[2] + 1.0

        image[22][0] = atom_with_image_dir[0] + 1.0
        image[22][1] = atom_with_image_dir[1] - 1.0
        image[22][2] = atom_with_image_dir[2] + 1.0

        image[23][0] = atom_with_image_dir[0] + 1.0
        image[23][1] = atom_with_image_dir[1] + 1.0
        image[23][2] = atom_with_image_dir[2] - 1.0

        image[24][0] = atom_with_image_dir[0] + 1.0
        image[24][1] = atom_with_image_dir[1] - 1.0
        image[24][2] = atom_with_image_dir[2] - 1.0

        image[25][0] = atom_with_image_dir[0] - 1.0
        image[25][1] = atom_with_image_dir[1] + 1.0
        image[25][2] = atom_with_image_dir[2] - 1.0

        image[26][0] = atom_with_image_dir[0] - 1.0
        image[26][1] = atom_with_image_dir[1] - 1.0
        image[26][2] = atom_with_image_dir[2] + 1.0
        return image

    def close_image_dis(self, single_atom):
        dis = np.zeros(27)
        image = AtomImage(self.atom_with_image, self.cell).image()
        image_car = direct_cartesian_transform(image, self.cell, 'DtoC')
        for d in range(27):
            dis[d] = distance(image_car[d], single_atom)
        close_image = np.amin(dis)
        return close_image
    
    # single atom being the first atom in the rigid body 
    def close_image_position(self, single_atom):
        dis = np.zeros(27)
        image = AtomImage(self.atom_with_image, self.cell).image()
        image_car = direct_cartesian_transform(image, self.cell, 'DtoC')
        for d in range(27):
            dis[d] = distance(image_car[d], single_atom)
        close_image_index = int(np.where(dis == np.amin(dis))[0][0])
        close_image_position = image_car[close_image_index]
        return close_image_position


def dis_check(atom_positions, cell, dis_limit):
    length = len(atom_positions)
    count2 = 0
    for i in range(1, length):
        count1 = 0
        for j in range(i):
            dis = AtomImage(atom_positions[i], cell).close_image_dis(atom_positions[j])
            if dis >= dis_limit:
                count1 += 1
            # else:
                # print(f'dis_check distance between atom{i+1} and atom{j+1} is {round(dis, 8)}')
        if count1 == i:
            count2 += 1
    if count2 == length - 1:
        return 'yes'
    else:
        return 'no'


def write_poscar(positions, compound, rigid_muti, ratio, sim_box):
    with open('POSCAR', 'w') as poscar:
        poscar.write('Good luck\n')
        poscar.write('1.0\n')
        for i in range(3):
            poscar.write(str(sim_box[i][0]) + '  ' + str(sim_box[i][1]) + '  ' + str(sim_box[i][2]) + '\n')
        compound_short = [compound[0]]
        compound_len = len(compound)
        for i in range(1, compound_len):
            check = 0
            for j in range(i):
                if compound[i] == compound[j]:
                    check = 1
            if check == 0:
                compound_short.append(compound[i])
        element_number = len(compound_short)
        for i in range(element_number):
            poscar.write(compound_short[i] + '  ')
        poscar.write('\n')
        compound_directory = dict(Counter(compound))
        compound_number = []
        for i in range(2):
            compound_number.append(compound_directory[compound_short[i]] * rigid_muti)
        compound_number.append(rigid_muti * ratio)
        for i in range(element_number):
            poscar.write(str(compound_number[i]) + '  ')
        poscar.write('\n')
        poscar.write('Car\n')
        positions_len = len(positions)
        for i in range(positions_len):
            x = round(positions[i][0], 10)
            y = round(positions[i][1], 10)
            z = round(positions[i][2], 10)
            poscar.write(str(x) + '  ' + str(y) + '  ' + str(z) + '\n')


def write_system_initial(system_initial):
    with open('system_initial', 'w') as sys_initial:
        sys_len = len(system_initial)
        for i in range(sys_len):
            x = str(round(system_initial[i][0], 10))
            y = str(round(system_initial[i][1], 10))
            z = str(round(system_initial[i][2], 10))
            sys_initial.write(x + '  ' + y + '  ' + z + '\n')


def read_system_initial(directory, total_atom):
    os.chdir(directory)
    with open('system_initial', 'r') as sys:
        line = sys.readlines()
        atom_coordinate = np.zeros((total_atom, 3))
        for i in range(total_atom):
            for j in range(3):
                atom_coordinate[i][j] = float(line[i].split()[j])
    return atom_coordinate


class WriteLog:
    def __init__(self, output_path):
        self.output_path = output_path

    def write_script_count(self):
        log_path = self.output_path + '/script_count'
        os.chdir(self.output_path)
        script1 = os.path.isfile(log_path)
        script = str(script1)
        if script == 'False':
            with open('script_count', 'w') as scrip:
                scrip.write('0')
        with open('script_count', 'r') as scrip_read:
            line_script = scrip_read.readlines()
            script_count = float(line_script[0].split()[0])
        script_count = script_count + 1
        with open('script_count', 'w') as scrip_write:
            scrip_write.truncate(0)
            scrip_write.write(f'{script_count}')
        return script_count
    
    def write_big_loop(self, script_count, internal_circulation):
        log_path = self.output_path + '/big_loop'
        os.chdir(self.output_path)
        script1 = os.path.isfile(log_path)
        script = str(script1)
        if script == 'False':
            with open('big_loop', 'w') as scrip:
                scrip.write('0')
        print(f'script_count:{script_count}')
        print(f'internal_circulation:{internal_circulation}')
        if script_count % internal_circulation == 0:
            with open('big_loop', 'r') as scrip_read:
                line_script = scrip_read.readlines()
                script_count = float(line_script[0].split()[0])
            script_count = script_count + 1
            with open('big_loop', 'w') as scrip_write:
                scrip_write.truncate(0)
                scrip_write.write(f'{script_count}')
        return script_count


def read_energy(keyword_file):
    keyword_text = 'sigma'
    with open(keyword_file, 'r') as outcar:
        lines = outcar.readlines()
        for line in reversed(lines):
            if keyword_text in line:
                energy = float(line.split()[6])
                break
    return energy


def metropolis_group_single(algorithm, energy, ref_energy, temp):
    accept = 0
    if algorithm == 'default':
        x = 1.0 / (1.0 + math.exp((energy - ref_energy) / temp))
        random_number = random.uniform(0, 1)
        if x > 0.5:
            accept = 1
        elif x == 1:
            accept = 2
        elif x > random_number:
            accept = 1
        else:
            accept = 0
    if algorithm == 'test':
        x = math.exp((ref_energy - energy) / temp)
        random_number = random.uniform(0, 1)
        if x > 1:
            accept = 1
        elif x == 1:
            accept = 2
        elif x > random_number:
            accept = 1
        else:
            accept = 0
    return accept


def setup_file_copy(my_cwd, sym_no, i, n):
    shutil.copy('system_initial', my_cwd + f'/sym{sym_no}' + f'/wyckoff{i + 1}' + f'/case{n + 1}' + '/trail')
    shutil.copy('INCAR', my_cwd + f'/sym{sym_no}' + f'/wyckoff{i + 1}' + f'/case{n + 1}' + '/trail')
    shutil.copy('POSCAR', my_cwd + f'/sym{sym_no}' + f'/wyckoff{i + 1}' + f'/case{n + 1}' + '/trail')
    shutil.copy('POTCAR', my_cwd + f'/sym{sym_no}' + f'/wyckoff{i + 1}' + f'/case{n + 1}' + '/trail')
    shutil.copy('KPOINTS', my_cwd + f'/sym{sym_no}' + f'/wyckoff{i + 1}' + f'/case{n + 1}' + '/trail')
    shutil.copy('vasp.sub', my_cwd + f'/sym{sym_no}' + f'/wyckoff{i + 1}' + f'/case{n + 1}' + '/trail')

    shutil.copy('INCAR_relax', my_cwd + f'/sym{sym_no}' + f'/wyckoff{i + 1}' + f'/case{n + 1}' + '/relax' + '/INCAR')
    shutil.copy('system_initial', my_cwd + f'/sym{sym_no}' + f'/wyckoff{i + 1}' + f'/case{n + 1}' + '/relax')
    shutil.copy('POSCAR', my_cwd + f'/sym{sym_no}' + f'/wyckoff{i + 1}' + f'/case{n + 1}' + '/relax')
    shutil.copy('POTCAR', my_cwd + f'/sym{sym_no}' + f'/wyckoff{i + 1}' + f'/case{n + 1}' + '/relax')
    shutil.copy('KPOINTS_relax', my_cwd + f'/sym{sym_no}' + f'/wyckoff{i + 1}' + f'/case{n + 1}' + '/relax' + '/KPOINTS')
    shutil.copy('vasp.sub', my_cwd + f'/sym{sym_no}' + f'/wyckoff{i + 1}' + f'/case{n + 1}' + '/relax')

    shutil.copy('functions.py', my_cwd + f'/sym{sym_no}' + f'/wyckoff{i + 1}' + f'/case{n + 1}')
    shutil.copy('interprocess.py', my_cwd + f'/sym{sym_no}' + f'/wyckoff{i + 1}' + f'/case{n + 1}')
    shutil.copy('mc.py', my_cwd + f'/sym{sym_no}' + f'/wyckoff{i + 1}' + f'/case{n + 1}')
    shutil.copy('metropolis.py', my_cwd + f'/sym{sym_no}' + f'/wyckoff{i + 1}' + f'/case{n + 1}')
    shutil.copy('rigid.py', my_cwd + f'/sym{sym_no}' + f'/wyckoff{i + 1}' + f'/case{n + 1}')
    shutil.copy('symdata.py', my_cwd + f'/sym{sym_no}' + f'/wyckoff{i + 1}' + f'/case{n + 1}')
    
    shutil.copy('INPUT', my_cwd + f'/sym{sym_no}' + f'/wyckoff{i + 1}' + f'/case{n + 1}')
    shutil.copy('structureupdate.py', my_cwd + f'/sym{sym_no}' + f'/wyckoff{i + 1}' + f'/case{n + 1}')
    shutil.copy('vasp2.sub', my_cwd + f'/sym{sym_no}' + f'/wyckoff{i + 1}' + f'/case{n + 1}')
    shutil.copy('readinginput.py', my_cwd + f'/sym{sym_no}' + f'/wyckoff{i + 1}' + f'/case{n + 1}')
    shutil.copy('randommove.py', my_cwd + f'/sym{sym_no}' + f'/wyckoff{i + 1}' + f'/case{n + 1}')

