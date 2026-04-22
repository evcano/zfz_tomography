import matplotlib.pyplot as plt
import numpy as np
import os
from glob import glob


files = glob("./*obs*.txt")

for fname in files:
    bname = os.path.basename(fname)
    X = np.loadtxt(fname)

    plt.plot(X[:,0], X[:,1], "-", label=bname)
    plt.plot(X[:,0], X[:,1], ".", c="k")

plt.xlim(1,14)
plt.ylim(0.5,5)
plt.legend()
plt.grid()
plt.show()
plt.close()
