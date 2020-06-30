import argparse
from pandas import *
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# Getting command-line options back
parser = argparse.ArgumentParser()
parser.add_argument('-d', dest='DESeq2', required=True) # DESeq2 csv result file
parser.add_argument('-e', dest='edgeR', required=True) # edgeR csv result file
parser.add_argument('-c1', dest='cond1', required=True) # first condition
parser.add_argument('-c2', dest='cond2', required=True) # second condition
options = parser.parse_args()

## lecture des fichiers
data01 = read_table(options.DESeq2,sep = ',')
data01 = data01[data01.padj <= 0.05]
data02 = read_table(options.edgeR,sep = ',')
data02 = data02[data02.FDR <= 0.05]

### sets au bon format pour venn2 et venn3
set1 = data01.index
set2 = data02.index

#legendes
label1 = 'DESeq2 (' + str(len(set1)) + ' DE genes)'
label2 = 'edgeR (' + str(len(set2)) + ' DE genes)'

# Plotting
plt.rcParams['figure.figsize'] = (10, 7)
venn2([set1, set2], (label1, label2))
#plt.title("Genes differentially expressed\n(WT-untagged Input vs WT-untagged IP)")
plt.title('Genes differentially expressed\n(' + options.cond1 + ' vs ' + options.cond2  + ')')
plt.show()