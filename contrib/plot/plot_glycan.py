import sys
import matplotlib
matplotlib.use("Qt4Agg")
matplotlib.interactive(True)
from matplotlib import pyplot as plt
from glypy import plot
from glypy.io import glycoct

try:
    print sys.argv
    structure = glycoct.loads(sys.argv[1])
    print(structure)
    plot.plot(structure)
    plt.show(block=True)
except Exception, e:
    print e
