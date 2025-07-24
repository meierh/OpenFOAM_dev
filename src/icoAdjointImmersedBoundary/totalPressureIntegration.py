import numpy as np
from scipy import integrate
U_m = 0.45
H = 0.41
f = lambda y,z : 0.5*(16*U_m*y*z*(H-y)*(H-z)/(pow(H,4)))**2
print(integrate.dblquad(f,0,H,0,H))

from sympy import*
from math import*
x=Symbol('x')
y=Symbol('y')
z=Symbol('z')
f = "0.5*(16*U_m*y*z*(H-y)*(H-z)/(pow(H,4)))**2"
inty=integrate(f,y)
intyz=integrate(inty,z)
y=H
z=H
f_HH = eval(intyz)
print(f_HH)
y=0
z=0
f_00 = eval(intyz)
print(f_00)
x=0
a=eval(abl)
x=6
b=eval(abl)
print(abl)
print(inty)
print(b-a)
