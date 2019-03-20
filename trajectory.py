#!/usr/bin/env python

import numpy as np
 
data = np.loadtxt("trajectory.txt", delimiter=' ', dtype=np.float)

#data = np.delete(data[:,[0,1]])
data = data.reshape(len(data), 5)[:,np.array([False, False, True, True, True])]

data = np.transpose(data)

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d 
# координаты точек
N = np.linspace(0,10,50)
x = data[0]
y = data[2]
z = data[1]
# построение и отображение графика
plt.gca(projection='3d')   # указываем тип
plt.xlabel('X')
plt.ylabel('Z')
plt.plot(x, z, y)                   # строим
plt.show()
