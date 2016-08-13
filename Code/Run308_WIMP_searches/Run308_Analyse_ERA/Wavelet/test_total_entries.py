
#Quick code to check entries of the trace_trigged_files
from ROOT import *
import PyROOTPlots as PyRPl
import glob

list_file = sorted(glob.glob("./ROOT_files/*root*"))

tg, fg = PyRPl.open_ROOT_object(list_file[0], "out")
th, fh = PyRPl.open_ROOT_object(list_file[1], "out")
ti, fi = PyRPl.open_ROOT_object(list_file[2], "out")
tj, fj = PyRPl.open_ROOT_object(list_file[3], "out")
tk, fk = PyRPl.open_ROOT_object(list_file[4], "out")


list_trees = [tg, th, ti, tj, tk]
sum_entries = 0


for tree in list_trees:
	sum_entries+=tree.GetEntries()
	print tree, tree.GetEntries()

print "total entries", sum_entries