from tkinter import X
from casecount import *
import matplotlib.pyplot as plt
from readinginput import *
from itertools import combinations, combinations_with_replacement

temp = 0.9
energy = -75.22656373
ref_energy = -77.40847807

x = 1.0 / (1.0 + math.exp((energy - ref_energy) / temp))
y = step_algorithm('normal', temp)

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
print(x)
print(y)