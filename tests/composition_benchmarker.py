from pygly2 import Composition

from common import load

import pstats
import cProfile

glyc = load("broad_n_glycan")

outfile = "test_reports/" + (Composition.__name__) + "_fragmentation_profile.stats"

cProfile.run("list(glyc.fragments('ABXY', max_cleavages=2))",
             filename="test_reports/" + (Composition.__name__) + "_fragmentation_profile.stats")


stat = pstats.Stats(outfile)
stat.strip_dirs().sort_stats("cumulative").print_stats()
