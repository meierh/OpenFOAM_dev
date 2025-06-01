import numpy as np
import matplotlib.pyplot as plt

def eval_fourier(a_0, list_ak, list_bk, phase, angle):
    val = 0
    val += a_0
    if len(list_ak)!=len(list_bk):
        raise Exception("Invalid argument")
    for c in range(len(list_ak)):
        val += list_ak[c]*np.cos((c+1)*angle+phase)
        val += list_bk[c]*np.sin((c+1)*angle+phase)
    return val

def eval_crossSec(a_0, list_ak, list_bk, phase, angle):
    r = eval_fourier(a_0, list_ak, list_bk, phase, angle)
    x = np.cos(angle)*r
    y = np.sin(angle)*r
    return x,y

a_0 = 0.05
list_ak = [0.025]
list_bk = [0]
phase = 0

angle = np.linspace(0,2*np.pi,100)
x_list = []
y_list = []
for a in angle:
    x,y = eval_crossSec(a_0,list_ak,list_bk,phase,a)
    x_list.append(x)
    y_list.append(y)
plt.plot(x_list,y_list)
plt.axis('equal')
plt.grid()
plt.show()
