from sympy.interactive.printing import init_printing

init_printing(use_unicode=False, wrap_line=False)
from sympy import pprint

import sympy.matrices as spm
from sympy.abc import s, n, a, t, m
from sympy import sin, cos, tan, pi
import numpy as np


deform = spm.Matrix([[1, -sin(pi/6), 0], [0, cos(pi/6), 0], [0, 0, 1]]) 
#translate = spm.Matrix([[1, 0, -a * (1 - sin(pi/6))], [0, 1, -a * cos(pi/6)], [0, 0, 1]])
#translate_back = spm.Matrix([[1, 0, a * (1 - sin(pi/6))], [0, 1, a * cos(pi/6)], [0, 0, 1]])
translate = spm.Matrix([[1, 0, -a], [0, 1, -a], [0, 0, 1]])
translate_back = spm.Matrix([[1, 0, a], [0, 1, a], [0, 0, 1]])
rot = spm.Matrix([[cos(t), -sin(t), 0], [sin(t), cos(t), 0], [0, 0, 1]])
flip = spm.Matrix([[((1 - m ** 2)/(1 + m ** 2)), ((2*m) / (1 + m ** 2)), 0], [((2*m)/ (1 + m **2)), ((m ** 2 - 1)/(1 + m ** 2)), 0], [0, 0, 1]])
#flip = spm.Matrix([[-1, 0, 0], [0, 1, 0], [0, 0, 1]])


main_func = spm.Matrix([s - (n * (s // n)), s // n, 1])
main_func_inv = spm.Matrix([1, n, 0]).transpose()

#full_transform = main_func_inv * deform.inv() * translate_back * flip * translate * deform * main_func 
full_transform = main_func_inv * translate_back * flip * translate * main_func 


angle = tan(-45 * pi / 180)
print(angle)
specific_transform = full_transform.subs(t, angle).subs(m, angle).subs(n, 5).subs(a, 2)
print(specific_transform[0, 0].subs(s, 14).evalf())

