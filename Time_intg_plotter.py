import numpy as np
from matplotlib import pyplot as plt

data_KE = np.loadtxt('KE.txt', delimiter=',')
dataDisp = np.loadtxt('disp.txt', delimiter=',')
plt.figure(1)
plt.plot(data_KE[:,0], data_KE[:,1])
plt.figure(2)
plt.plot(dataDisp[:,0], dataDisp[:,1],'-o')
plt.show()