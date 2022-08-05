from casecount import *
import matplotlib.pyplot as plt
from readinginput import *
from itertools import combinations, combinations_with_replacement

array = np.array([[1, 0, 1], [2, 2, 1]])
y = rigid_select(rigid_type)(1, 1.0).center(array)
z = Tetrahedron(1, 1.0).center(array)
print(y)