from matplotlib import pylab as plab
import numpy as np
import sys

f = open(sys.argv[1],'r+')
buf = f.readlines()

for i in range(len(buf)):
    buf[i] = buf[i].split()

rows = len(buf)
cols = len(buf[0])

z = np.zeros([rows,cols])

for i in range(rows):
    for j in range(cols):
        z[i][j] = float(buf[i][j].strip(','))

x = np.linspace(0,1,rows)
y = np.linspace(0,1,cols)
plab.pcolormesh(x,y,z)
plab.set_cmap('coolwarm')
plab.colorbar()
plab.show()
