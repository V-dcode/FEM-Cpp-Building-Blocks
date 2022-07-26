import numpy as np
from matplotlib import pyplot as plt

data_KE = np.loadtxt('KE.txt', delimiter=',')
plt.plot(data_KE[:,0], data_KE[:,1])
plt.show()