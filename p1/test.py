from tkinter import X
from casecount2 import *
import matplotlib.pyplot as plt
from readinginput import *
from itertools import combinations, combinations_with_replacement

# temp = 0.2
# energy = -75.22656373
# ref_energy = -77.40847807

# x = 1.0 / (1.0 + math.exp((energy - ref_energy) / temp))
# y = step_algorithm('normal', temp)

# v = np.random.rand(3)
# z = np.linalg.norm(v)
# v_hat = v / np.linalg.norm(v)
# print(v)
# norm_v = (v[0]**2 + v[1]**2 + v[2]**2)**(0.5)
# print(f'norm_v:{norm_v}')
# print(z)
# print(v_hat)
# norm_v_hat = (v_hat[0]**2 + v_hat[1]**2 + v_hat[2]**2)**(0.5)
# print(f'norm_v_hat:{norm_v_hat}')
# print(x)
# print(y)

# x = SymCasesP1(rigid_type, np.array([1, 3]), 1, 3)
# y = x.case_list()
# print(y)

# a = np.array([[1, 2, 3],[4,5,6],[7,8,9]])
# b = np.append(a, [[4, 4, 4]], axis = 0)
# print(b)

a = [['x', 'y', 'z']]
b = a*4
print(b)