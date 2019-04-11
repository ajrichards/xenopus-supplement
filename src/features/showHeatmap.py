#!/usr/bin/env python
"""
create the manuscript figure summarizing features

volcano plot - showing sig genes/isoforms
clustering heatmap - showing how genes/isoforms cluster


"""

import sys,os,re,Image,csv,cPickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib as mpl
from htsint.tools import Heatmap
from matplotlib import rc
mpl.rcParams['text.usetex'] = True

## variables
homeDir = os.path.join(os.path.expanduser("~"),"sequencing","xenopus")
fontSize = 11
fontName = 'sans-serif'
useColor = True
colorList = ['#1122FF','#FF5500']
if useColor == True:
    myCmap = mpl.cm.gist_heat
else:
    myCmap = mpl.cm.gray

## load the dif exp data
deseqColumns = ['baseMean' 'log2FoldChange' 'lfcSE' 'stat' 'pvalue' 'padj']
dfeColumns = ["A","B","C","D","E","F","G","H"]
threshold = 0.05
_threshold = re.sub("\.","",str(threshold))
results = {}

for assembly in ['dn','gg','ref']:
    outfid = open("../filtered/filtered-%s-%s.pickle"%(assembly,_threshold),'r')
    results[assembly] = cPickle.load(outfid)
    outfid.close()
    print assembly, results[assembly]['dfeMat'].shape

## remove gene versions
for k,gene in enumerate(results['dn']['id']):
    results['dn']['id'][k] = re.sub("\.\d+$","",gene)

## only report full non-repeated gene symbols for visualization
usedInds = []
withRepeat = []
for g,gene in enumerate(np.unique(results['dn']['id'])):
    indx = np.where(results['dn']['id'] == gene)[0][0]
    if re.search("LOC",gene):
       pass
    else:
        usedInds.append(indx)
    
filter1 = np.array(usedInds)
origSize = results['dn']['id'].size
for key,item in results['dn'].iteritems():
    if key in ['dfeMat','deseqMat']:
        results['dn'][key] = results['dn'][key][filter1,:]
    else:
        results['dn'][key] = results['dn'][key][filter1]

print('keeping %s/%s genes from the dn results -- locus filter and multitranscript filter'%(len(usedInds),origSize))

colLabels = ['A','B','C','D','E','F','G','H']
_endurant = ['E','A','D','G']
_nonendurant = ['B','F','C','H']
nonendurant = [colLabels.index(lab) for lab in _nonendurant]
endurant = [colLabels.index(lab) for lab in _endurant]

## create the heatmap
for assembly in ['dn','gg','ref']:
    heatmapPath= os.path.join(".","heatmap-%s.png"%(assembly))

    if assembly == 'dn':
        rowLabels = results[assembly]['id']
        dfeMat = results[assembly]['dfeMat']
        width = 7
        height = 7
        rowFont = 6
        fSize = fontSize
    else:
        rowLabels = results[assembly]['id']
        dfeMat = results[assembly]['dfeMat']
        width = 7
        height = 8
        rowFont = 13
        fSize = fontSize + 4
    
    mpl.rcParams['text.usetex'] = False
    hm = Heatmap(dfeMat,rowLabels,colLabels,fontName=fontName,fontSize=fSize,
                 hpad=0.19,width=width,height=height)
    hm.draw_heatmap(cmap='uy',clabels=True,rlabels=True,rowFont=rowFont)
    endurantInds = [np.where(np.array(hm.indx['1']) == i)[0][0] for i in endurant]
    nonendurantInds = [np.where(np.array(hm.indx['1']) == i)[0][0] for i in nonendurant]
                        
    for t,txt in enumerate(hm.ax3.get_xticklabels()):
        if t in nonendurantInds:
            txt.set_color(colorList[1])
        if t in endurantInds:
            txt.set_color(colorList[0])

    ## format color bar labels
    for t in hm.ax4.get_xticklabels():
        t.set_fontsize(fSize-2)
        t.set_fontname(fontName)
    for t in hm.ax4.get_yticklabels():
        t.set_fontsize(fSize-2)
        t.set_fontname(fontName)
            
    hm.save(heatmapPath,dpi=500)
