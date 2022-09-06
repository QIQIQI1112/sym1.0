from ase.spacegroup import crystal
from readinginput import *
from randommove import *

root_path = os.getcwd()
relax_path = root_path + '/relax'
trail_path = root_path + '/trail'
output_path = root_path + '/output'
good_ones_relaxed_path = output_path + '/good_ones_relaxed'

contcar = ReadPOSCAR('CONTCAR', relax_path, root_path)
after_relax = contcar.positions
cell_after = contcar.lattice

after_relax_car = direct_cartesian_transform(after_relax, cell_after, 'DtoC')

update = sturcture_update_p1(after_relax_car, rigid_no, cell_after)

os.chdir(trail_path)
write_poscar_p1(update, rigid_no, cell_after)

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

# relaxed vs script_count
os.chdir(output_path)
with open('script_count', 'r') as scrip_read:
    line_script = scrip_read.readlines()
    script_count = float(line_script[0].split()[0])
with open('energy_relax', 'a') as en_relax:
    en_relax.write(f'{script_count}  {energy}' + '\n')

os.chdir(root_path)