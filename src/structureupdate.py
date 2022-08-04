from ase.spacegroup import crystal
from readinginput import *
from randommove import *

root = os.getcwd()
os.chdir(root + '/trail')
trail_path = os.getcwd()
total_atom = rigid_atom + single_Atom
system_initial = read_system_initial(trail_path, total_atom)

os.chdir(root + '/relax')
relax_path = os.getcwd()
poscar = ReadPOSCAR('POSCAR', relax_path)
outcar = ReadPOSCAR('OUTCAR', relax_path)
before_relax = poscar.positions
cell_before = poscar.lattice
after_relax = outcar.positions
cell_after = outcar.lattice

system_initial_after = sturcture_update(system_initial, before_relax, after_relax, cell_before, cell_after)

compound = [*compound_rigid, *(compound_single_atom * single_Atom)]
system_initial_dir = direct_cartesian_transform(system_initial_after, cell_after, 'CtoD')
system_initial_tuple = []
for k in range(total_atom):
    system_initial_tuple.append(tuple(system_initial_dir[k]))
system_transform = crystal(compound, system_initial_tuple, spacegroup=sym_no, cell=cell_after)

os.chdir(trail_path)
write_poscar(system_transform.positions, compound, rigid_muti, ratio, cell_after)
write_system_initial(system_initial_after)

os.chdir(root)