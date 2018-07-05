import sys
sys.path.append("/usr/local/visit/current/linux-x86_64/lib/site-packages")
import visit

fname = "viz01:/bg01/homescinet/p/peltier/hsalehip/bgq-scratch/holm/Re-6000-Ri-016-Pr-8/pp/socppholm.nek5000"
visit.OpenDatabase(fname, 0)
visit.AddPlot("Pseudocolor", "temperature", 1, 1)
visit.AddOperator("Slice", 1)
visit.DrawPlots()
