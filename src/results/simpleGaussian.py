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
metaDataDir = os.path.join(os.path.expanduser("~"),"sequencing","xenopus-data","meta-data")
from metadatalib import read_meta_data
mdata,X,T,features,targets = read_meta_data()
sampleList = np.array(["A","B","C","D","E","F","G","H"])
T = T[:,3]

## print previous classes
print("\nNon-endurant: C,D,F,G")
print("Endurant: A,B,E,H")

def run_mm(T,transform=False,standardize=False):
    y = T.copy()
    if transform == True:
        y = np.log10(y)
    if standardize == True:
        y = scale(y)

    mm = GMM(n_components=3,covariance_type='tied', init_params='wc', n_iter=20)
    mm.fit(y)
    labels = mm.predict(y)
    means = np.array([y[np.where(labels==i)[0]].mean() for i in [0,1,2]] )
    classes = ['low','medium','high']
    sortedInds = np.argsort(means)

    return labels,classes,sortedInds

## print the results
print "\nuntransformed, unstandardized"
labels,classes,sortedInds = run_mm(T,transform=False,standardize=False)
for i,s in enumerate(sortedInds):
    print classes[i], sampleList[np.where(labels==s)[0]]

print "\nuntransformed, standardized"
labels,classes,sortedInds = run_mm(T,transform=False,standardize=True)
for i,s in enumerate(sortedInds):
    print classes[i], sampleList[np.where(labels==s)[0]]

print "\ntransformed, unstandardized"
labels,classes,sortedInds = run_mm(T,transform=True,standardize=False)
for i,s in enumerate(sortedInds):
    print classes[i], sampleList[np.where(labels==s)[0]]

print "\ntransformed, standardized"
labels,classes,sortedInds = run_mm(T,transform=True,standardize=True)
for i,s in enumerate(sortedInds):
    print classes[i], sampleList[np.where(labels==s)[0]]


### untransformed, standardized
#print "untransformed, unstandardized"
#for i,s in enumerate(sortedInds):
#    print classes[i], sampleList[np.where(labels==s)[0]]


#T = np.log10(T)
#X = np.log10(X)
#T = scale(T)
#X = scale(X)




