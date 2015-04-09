import sys
import matplotlib
matplotlib.use("Qt4Agg")
matplotlib.interactive(True)
from matplotlib import pyplot as plt
from pygly2 import plot
from pygly2.io import glycoct

try:
    print sys.argv
    structure = glycoct.loads(sys.argv[1]).next()
    print(structure)
    plot.plot(structure)
    plt.show(block=True)
except Exception, e:
    print e
