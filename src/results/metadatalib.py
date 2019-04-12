#!/usr/local/bin/python
"""
explore the covariates (experiment meta-data)
"""

import os,sys,csv
import numpy as np
metaDataDir = os.path.join(os.path.expanduser("~"),"sequencing","xenopus","meta-data")

def read_meta_data():
    """
    return a data structure for the meta-data
    """

    filePath = os.path.join(metaDataDir,"rna_individual_measures.csv")
    fid = open(filePath,'rU')
    reader = csv.reader(fid)
    header = reader.next()

    mdata = {'id':[],'mobility':[],'svl':[],'mass':[],'head-length':[],'head-width':[], 
             'head-height':[],'jaw-length':[],'femur':[],'tibia':[],'foot':[], 
             'toe':[],'hindlimb-length':[],'humerus':[],'radius':[],'hand':[],'finger':[], 
             'fll':[],'ilium-length':[],'ilium-width':[],'avg-v':[],'vmax':[],'amax':[], 
             'tmax-stam1':[],'dmax-stam1':[]
             }

    for linja in reader:
        mdata['id'].append(linja[0])
        mdata['mobility'].append(linja[1])
        mdata['svl'].append(linja[2])
        mdata['mass'].append(linja[3])
        mdata['head-length'].append(linja[4])
        mdata['head-width'].append(linja[5])
        mdata['head-height'].append(linja[6])
        mdata['jaw-length'].append(linja[7])
        mdata['femur'].append(linja[8])
        mdata['tibia'].append(linja[9])
        mdata['foot'].append(linja[10])
        mdata['toe'].append(linja[11])
        mdata['hindlimb-length'].append(linja[12])
        mdata['humerus'].append(linja[13])
        mdata['radius'].append(linja[14])
        mdata['hand'].append(linja[15])
        mdata['finger'].append(linja[16])
        mdata['fll'].append(linja[17])
        mdata['ilium-length'].append(linja[18])
        mdata['ilium-width'].append(linja[19])
        mdata['avg-v'].append(linja[20])
        mdata['vmax'].append(linja[21])
        mdata['amax'].append(linja[22])
        mdata['tmax-stam1'].append(linja[23])
        mdata['dmax-stam1'].append(linja[24])
        
    for key,value in mdata.iteritems():
        if key not in ['id','mobility']:
            mdata[key] = [float(i) for i in value]
        mdata[key] = np.array((mdata[key],)).transpose()

    predictors = ['svl','mass','head-length','head-width','head-height',
                  'jaw-length','femur','tibia','foot','toe','hindlimb-length',
                  'humerus','radius','hand','finger','fll','ilium-length','ilium-width']

    X = None
    for p in predictors:
        if X == None:
            X = mdata[p]
        else:
            X = np.hstack([X,mdata[p]])

    print X.shape

    T = None

    ## ignoring 'avg-v'
    targets = ['vmax','amax','tmax-stam1','dmax-stam1']
    for t in targets:
        if T == None:
            T = mdata[t]
        else:
            T = np.hstack([T,mdata[t]])
  
    fid.close()
    return mdata,X,T,predictors,targets
