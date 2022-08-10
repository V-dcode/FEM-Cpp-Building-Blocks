import numpy as np
from matplotlib import pyplot as plt

data_KE = np.loadtxt('KE_Truss.txt', delimiter=',')
data_Disp = np.loadtxt('Disp_Truss.txt', delimiter=',')
data_Conn = np.loadtxt('Conn_Truss.txt', delimiter=',', dtype=int)
Deform_x = data_Disp[:,0]+data_Disp[:,3]
Deform_y = data_Disp[:,1]+data_Disp[:,4]
# Deform_Config = [Deform_x, Deform_y]
plt.figure(1)
plt.plot(data_KE[:,0], data_KE[:,1])
plt.figure(2)
for elem in data_Conn:
    plt.plot([data_Disp[elem[0],0], data_Disp[elem[1],0]], [data_Disp[elem[0],1], data_Disp[elem[1],1]], 'k')
    plt.plot([Deform_x[elem[0]], Deform_x[elem[1]]], [Deform_y[elem[0]], Deform_y[elem[1]]], 'r')
plt.legend(['Undeformed', 'Deformed'])
plt.show()