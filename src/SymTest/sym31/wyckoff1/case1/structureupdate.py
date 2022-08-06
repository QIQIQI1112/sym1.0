from ase.spacegroup import crystal
from readinginput import *
from randommove import *

root_path = os.getcwd()
relax_path = root_path + '/relax'
trail_path = root_path + '/trail'
output_path = root_path + '/output'
good_ones_relaxed_path = output_path + '/good_ones_relaxed'

total_atom = rigid_atom + single_Atom
system_initial = read_system_initial(relax_path, total_atom)

poscar = ReadPOSCAR('POSCAR', relax_path, root_path)
contcar = ReadPOSCAR('CONTCAR', relax_path, root_path)
before_relax = poscar.positions
cell_before = poscar.lattice
after_relax = contcar.positions
cell_after = contcar.lattice

system_initial_after = sturcture_update(rigid_atom, system_initial, before_relax, after_relax, cell_before, cell_after)

compound = [*compound_rigid, *(compound_single_atom * single_Atom)]
system_initial_dir = direct_cartesian_transform(system_initial_after, cell_after, 'CtoD')
system_initial_tuple = []
for k in range(total_atom):
    system_initial_tuple.append(tuple(system_initial_dir[k]))
system_transform = crystal(compound, system_initial_tuple, spacegroup=sym_no, cell=cell_after)

os.chdir(trail_path)
write_poscar(system_transform.positions, compound, rigid_muti, ratio, cell_after)
write_system_initial(system_initial_after)

# save good structures
isdir1 = os.path.isdir(good_ones_relaxed_path)
isdir = str(isdir1)
if isdir == 'False':
    os.chdir(output_path)
    os.mkdir('good_ones_relaxed')
os.chdir(relax_path)
energy1 = read_energy('OUTCAR')
energy = round(float(energy1), 4)
list_file = os.listdir(good_ones_relaxed_path)
number_files = len(list_file)
if number_files < 11:
    shutil.copytree(relax_path, good_ones_relaxed_path + '/' + str(-energy))
elif number_files == 11:
    file = []
    for i in range(number_files):
        file.append(float(list_file[i]))
    file_sorted = sorted(file)
    shutil.rmtree(good_ones_relaxed_path + f'/{str(file_sorted[0])}')

os.chdir(root_path)