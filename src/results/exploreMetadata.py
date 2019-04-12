#!/usr/local/bin/python
"""
explore the covariates (experiment meta-data)

(1) use unsupervised learning to cluster the data
(2) determine categories for the data
(3) use categories for supervised learing methods

avg-v --- average swim speed
vmax --- maximum swim speed
amax --- maximum acceleration
tmax-stam1 --- maximum stamina
dmax-stam1 --- maximum distance jumped 

"""

import os,sys,csv
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from matplotlib.transforms import offset_copy
from sklearn import datasets, svm
from sklearn.svm import SVR
from sklearn.mixture import GMM
from sklearn.preprocessing import scale
from sklearn.feature_selection import SelectPercentile, f_classif
from htsint.stats import get_silhouette_values

## control fonts
fontSize = 11
fontName = 'sans-serif'
colors = ['black','#1122FF','#FF5500']

## append the meta-data dir
metaDataDir = os.path.join(os.path.expanduser("~"),"sequencing","xenopus","meta-data")
from metadatalib import read_meta_data
mdata,X,T,features,targets = read_meta_data()
sampleList = ["A","B","C","D","E","F","G","H"]
#T = np.log10(T)
#X = np.log10(X)
T = scale(T)
X = scale(X)

## check for correlation among target variables
fig1 = plt.figure(figsize=(8,4))
comparisons = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
plotNum = 0
buff = 0.05
for comp in comparisons:
    plotNum += 1
    ax = fig1.add_subplot(2,3,plotNum)
    indX,indY = comp
    ax.plot(T[:,indX],T[:,indY],'bo')
    ax.set_xlabel(targets[indX],fontsize=fontSize)
    ax.set_ylabel(targets[indY],fontname=fontName)
    ax.set_yticks([])
    ax.set_xticks([])
    xbuff = buff * (T[:,indX].max() - T[:,indX].min())
    ybuff = buff * (T[:,indY].max() - T[:,indY].min())
    ax.set_ylim([T[:,indY].min() - ybuff, T[:,indY].max() + ybuff])
    ax.set_xlim([T[:,indX].min() - xbuff, T[:,indX].max() + xbuff])

    pcc,pvalue = pearsonr(T[:,indX],T[:,indY])
    ax.set_title(r"r=%s"%round(pcc,2),fontsize=fontSize,fontname=fontName)
    ax.set_aspect(1./ax.get_data_ratio())
    fig1.subplots_adjust(hspace=0.0005,wspace=0.1)

plt.tight_layout()
plt.savefig('target_variable_correlation.png',bbox_inches='tight',dpi=300)

## we thus remove vmax and dmax-stam1
#toRemove = ['vmax','dmax-stam1']
#indsToRemove = [targets.index(tr) for tr in toRemove]
#indsToKeep = list(set(range(len(targets))).difference(set(indsToRemove)))
#for tr in toRemove:
#    targets.remove(tr)
#T = T[:,indsToKeep]

def create_gmm_subplot(fig,k,labels,plotNum,gmm):
    comparisons = [(3,0),(3,1),(3,2),(0,1),(0,2),(2,1)]
    buff = 0.15
    for comp in comparisons:
        plotNum += 1
        ax = fig.add_subplot(2,6,plotNum)
        indX,indY = comp
        letterAdj = 0.01 * (T[:,indY].max() - T[:,indY].min())
        for c in range(k):
            clusterInds = np.where(labels==c)[0]
            
            for cind in clusterInds:
                color = colors[labels[cind]]
                m = sampleList[cind]
                ax.plot(T[cind,indX],T[cind,indY],'wo',ms=10)
                ax.text(T[cind,indX],T[cind,indY]-letterAdj,m,
                        color=color,ha="center",va="center",
                        fontsize=fontSize-3)

        ax.set_xlabel(targets[indX],fontsize=fontSize-2,fontname=fontName)
        ax.set_ylabel(targets[indY],fontsize=fontSize-2,fontname=fontName)
        ax.set_yticks([])
        ax.set_xticks([])
        xbuff = buff * (T[:,indX].max() - T[:,indX].min())
        ybuff = buff * (T[:,indY].max() - T[:,indY].min())
        ax.set_ylim([T[:,indY].min() - ybuff, T[:,indY].max() + ybuff])
        ax.set_xlim([T[:,indX].min() - xbuff, T[:,indX].max() + xbuff])
        ax.set_aspect(1./ax.get_data_ratio())
        
        ## calculate silhouete value
        if plotNum in [2,8]:
            silValues = get_silhouette_values([T],[labels],resultsType='raw')
            avgSilVal = silValues['0'].mean()
            ax.set_title("k=%s,silhouette-value=%s"%(k,round(avgSilVal,2)),ha='right',fontsize=fontSize,fontname=fontName)

## use a GMM to fit the target variables into classes
## 'spherical' (.49,.56), 'diag', 'tied', 'full']
print("Running GMM on %s"%targets)
plt.clf()
fig = plt.figure(figsize=(7,3))
for k in [2]:  ## missing statistical power to do more than 2
    if k == 2:
        plotNum = 0
    elif k == 3:
        plotNum = 6

    mm = GMM(n_components=k,covariance_type='tied', init_params='wc', n_iter=20)
    mm.fit(T)
    labels = mm.predict(T)
    probas = mm.predict_proba(T)
    create_gmm_subplot(fig,k,labels,plotNum,mm)

#fig.suptitle('blah')
fig.subplots_adjust(hspace=0.01,wspace=0.3)
plt.savefig('gmm_on_targets.png',bbox_inches='tight',dpi=300)

## plot the mixture model classes
plt.clf()
fig = plt.figure()
samples = mm.sample(1000)
slabels = mm.predict(samples)
comparisons = [(3,0),(3,1),(3,2),(0,1),(0,2),(2,1)]

plotNum = 0
for comp in comparisons:
    plotNum += 1
    ax = fig.add_subplot(1,6,plotNum)
    indX,indY = comp
    for c in np.sort(np.unique(slabels)):
        clusterInds = np.where(slabels==c)[0]
        ax.scatter(samples[clusterInds,indX],samples[clusterInds,indY],c=colors[c],s=1,edgecolor='none')
    ax.set_xlabel(targets[indX],fontsize=fontSize)
    ax.set_ylabel(targets[indY],fontname=fontName)
    ax.set_yticks([])
    ax.set_xticks([])
    xbuff = buff * (samples[:,indX].max() - samples[:,indX].min())
    ybuff = buff * (samples[:,indY].max() - samples[:,indY].min())
    ax.set_ylim([samples[:,indY].min() - ybuff, samples[:,indY].max() + ybuff])
    ax.set_xlim([samples[:,indX].min() - xbuff, samples[:,indX].max() + xbuff])

    ax.set_aspect(1./ax.get_data_ratio())
    fig.subplots_adjust(hspace=0.05,wspace=0.3)

plt.tight_layout()
plt.savefig('gmm_classes.png',bbox_inches='tight',dpi=300)

## print out the classes
#classes = np.array(['non-endurant-1','non-endurant-2','endurant'])
#means = np.array([T[np.where(labels==i)[0],3].mean() for i in [0,1,2]])
classes = np.array(['non-endurant','endurant'])
means = np.array([T[np.where(labels==i)[0],3].mean() for i in [0,1]])
sortedInds = np.argsort(means)

for i,sample in enumerate(sampleList):
    print("   * ``%s`` = **%s** (%s)"%(sample,classes[np.where(sortedInds==labels[i])[0]],
                                       round(probas[i].max(),4)))
print("\n")

## feature selection
plt.clf()
fig = plt.figure(figsize=(8,4))
ax = fig.add_subplot(111)

indicesX = np.arange(X.shape[-1])
y = labels

## univariate feature selection
width = 0.3
selector = SelectPercentile(f_classif, percentile=25)
selector.fit(X, y)
scores = -np.log10(selector.pvalues_)
scores /= scores.max()
ax.bar(indicesX, scores, width=width,
       label=r'Univariate score ($-Log(p_{value})$)', color='g')

#clf_selected = svm.SVC(kernel='linear')
#clf_selected.fit(selector.transform(X), y)

#svm_weights_selected = (clf_selected.coef_ ** 2).sum(axis=0)
#svm_weights_selected /= svm_weights_selected.max()

#ax.bar(indicesX[selector.get_support()] + width, svm_weights_selected, width=width,
#       label='SVM weights after selection', color='b')

clf = svm.SVC(kernel='linear')
clf.fit(X,y)
svm_weights = (clf.coef_ ** 2).sum(axis=0)
svm_weights /= svm_weights.max()
ax.bar(indicesX + width, svm_weights, width=width, label='SVM weight', color='r')

ax.set_title("Feature selection")
ax.set_xticks(indicesX+width)
ax.set_xticklabels(features,rotation=45,fontsize=7)

ax.set_xlabel('Feature')
ax.set_yticks(())
ax.axis('tight')
ax.legend(loc='upper right')

#print svm_weights.tolist()
plt.savefig('feature_selection.png',bbox_inches='tight',dpi=300)

combined = svm_weights + scores
combinedInds = np.argsort(combined)

rank = 0
for i in combinedInds[::-1]:
    rank += 1
    print "   %s."%rank, "``%s``"%features[i], 
    print "%s,%s,%s"%(round(scores[i],2), round(svm_weights[i],2), round(combined[i],2))


# use only the first three features to classify the frogs
## we thus remove vmax and dmax-stam1
keptFeatures = ['head-height','hand','radius']
keptInds = [features.index(kf) for kf in keptFeatures]
selectedX = X[:,keptInds]

selectedClf = SVR(kernel='rbf',degree=3)
selectedClf.fit(X, y)
newLabels = selectedClf.predict(X)
print("actual labels %s"% labels)
print("predicted labels %s"% [int(round(l)) for l in newLabels])


#create_subplot(ax,xaxis,t,t_predict7,'SVM-Regression')
#fig.savefig('via-ml-svmreg.png',dpi=200)


#toRemove = ['vmax','dmax-stam1']
#indsToRemove = [targets.index(tr) for tr in toRemove]
#indsToKeep = list(set(range(len(targets))).difference(set(indsToRemove)))
#for tr in toRemove:
#    targets.remove(tr)


