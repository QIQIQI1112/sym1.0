from casecount import *
import matplotlib.pyplot as plt
from itertools import combinations, combinations_with_replacement

# cell = np.array([[6.9336127435, 0, 0],[0, 6.9336127435, 0], [0, 0, 6.9336127435]])
# position = np.array([[0.0,  5.6786288369,  5.567691033],
# [1.6329939232,  6.6885531545, 5.0078808425],
# [0.0,  5.4603956831,  7.5557498905],
# [0.0,  3.8770133555,  4.6992525566],
# [-1.6329939232,  6.6885531545,  5.0078808425],
# [0.0,  3.3212005041,  1.6848678967],
# [5.6439607732,  4.7703255675,  6.4551934642]])

# aa = move_atoms_into_box(position, cell)
# print(aa)

# def aimfunction(x, t):
#     y = ((2.0 * math.pi * t)**(- 3.0 / 2.0)) * math.exp((-x**2) / (2.0 * (t**2)))
#     return y


# def a2(x, t):
#     y = ((2.0 * t * math.pi)**(- 3.0 / 2.0)) * math.exp((-x**2) / (2.0 * (t**2)))
#     return y

# def a3(x, t):
#     y = t / math.log(x)
#     return y

# # x = 0.5
# x = [i / 1.0 for i in range(1, 1000)]
# y = [1.0 for i in range(999)]
# for i in range(0, 999):
#     y[i] = a3(x[i], 0.5)

# plt.plot(x, y)
# plt.show()

#====================
# sym = SymCases('Tetrahedron', np.array([1, 2, 1]), 1, 3)
# aa = sym.case_count(1, 150, 0)
# print(aa)
# print('===================')
# for i in aa:
#     print(f'{i} {aa[i]}')

t = 1.0

print(1 % 50)

