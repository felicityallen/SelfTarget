import io, os, sys, shutil, csv

from selftarget.plot import setFigType
from selftarget.util import setPlotDir
from selftarget.data import setHighDataDir

sys.path.append('i1')
from plot_i1_summaries import runAnalysis as plotI1
sys.path.append('kl_comparisons')
from plot_kl_analysis import runAnalysis as plotKLCmp
sys.path.append('microhomology')
from plot_mh_analysis import runAnalysis as plotMH
sys.path.append('microhomology_mismatch')
from plot_mh_mismatch_frequencies import runAnalysis as plotMhMismatch
sys.path.append('scaffold_compare')
from plot_old_new import runAnalysis as plotKLOldNew
sys.path.append('indel_details')
from plot_pie_indel_summaries import runAnalysis as plotIndelDetails
sys.path.append('../indel_prediction/model_testing')
from plot_old_new_predictions import runAnalysis as plotPred

setFigType('png')
setPlotDir('/results/plots')
setHighDataDir('/data/summary_data')

plotI1()
plotKLCmp()
plotMH()
plotMhMismatch()
plotKLOldNew()
plotIndelDetails()
plotPred()


        


