#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 16:58:12 2018

@author: onerva

This file contains parameters used for the bipartite analysis of the Climate 
Companions organization.
"""
import matplotlib.pylab as plt
import matplotlib.colors
import numpy as np
from matplotlib import cm

def categorical_cmap(nc, nsc, cmap, continuous=False):
    """
    A function for creating a tailored colormap. Original code from Stack Overflow: 
    https://stackoverflow.com/questions/47222585/matplotlib-generic-colormap-from-tab10.
    """
    if nc > plt.get_cmap(cmap).N:
        raise ValueError("Too many categories for colormap.")
    if continuous:
        ccolors = plt.get_cmap(cmap)(np.linspace(0,1,nc))
    else:
        ccolors = plt.get_cmap(cmap)(np.arange(nc, dtype=int))
    cols = np.zeros((nc*nsc, 3))
    for i, c in enumerate(ccolors):
        chsv = matplotlib.colors.rgb_to_hsv(c[:3])
        arhsv = np.tile(chsv,nsc).reshape(nsc,3)
        arhsv[:,1] = np.linspace(chsv[1],0.25,nsc)
        arhsv[:,2] = np.linspace(chsv[2],1,nsc)
        rgb = matplotlib.colors.hsv_to_rgb(arhsv)
        cols[i*nsc:(i+1)*nsc,:] = rgb       
    cmap = matplotlib.colors.ListedColormap(cols)
    return cmap

# input
companyInputPath = '/media/onerva/KINGSTON/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Members_alias_271018.csv'
eventInputPath = '/media/onerva/KINGSTON/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Events_alias_271018.csv'
linkInputPaths = ['/media/onerva/KINGSTON/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Links_all_271018.csv']#,
#                  '/media/onerva/KINGSTON/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Links_2013_2014.csv',
#                  '/media/onerva/KINGSTON/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Links_2015_2016.csv',
                  
#                  '/media/onerva/KINGSTON/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Links_all_271018.csv']
        
        #'/media/onerva/KINGSTON/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Links_all_271018.csv']
        
        #'/media/onerva/KINGSTON/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Links_2015_2016.csv',
         #         '/media/onerva/KINGSTON/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Links_2017_2018.csv']
        
#        '/media/onerva/KINGSTON/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Links_2011.csv',
#                  '/media/onerva/KINGSTON/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Links_2012.csv',
#                  '/media/onerva/KINGSTON/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Links_2013.csv',
#                  '/media/onerva/KINGSTON/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Links_2014.csv',
#                  '/media/onerva/KINGSTON/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Links_2015.csv',
#                  '/media/onerva/KINGSTON/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Links_2016.csv',
#                  '/media/onerva/KINGSTON/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Links_2017.csv',
#                  '/media/onerva/KINGSTON/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Links_2018.csv',
#                  '/media/onerva/KINGSTON/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Links_2011_2012.csv',
#                  '/media/onerva/KINGSTON/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Links_2013_2014.csv',
#                  '/media/onerva/KINGSTON/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Links_2015_2016.csv',
#                  '/media/onerva/KINGSTON/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Links_2017_2018.csv',
#                  '/media/onerva/KINGSTON/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Links_all_271018.csv']
                  
                  
companyColumnNames = ['Alias:', 'Member:']
eventColumnNames = ['Alias:', 'Event:']
linkColumnNames = ['Event:', 'Participant:']

years = ['all']#['11-12','13-14','15-16','17-18','all']

# distributions and binning
nDegreeBins = 20
cliqueHeatmapBottomBins = np.arange(0.5,16.5,1)
cliqueHeatmapTopBins = np.arange(0.5,68.5,1)
nRichnessBins = 5

lowDegreePercentile = 25
highDegreePercentile = 75

# comparison against random

nRandomIterations = 1000
nRandomBins = 20

# visualization
cmap = 'cool'
colors = cm.get_cmap(cmap, lut=2)
topColor = colors(0)
bottomColor = colors(1)

tags = ['C','CB','Con','CS','E','ENGO','F','H','HKI','II','In','L','Me','N','R','S','SD','Se','T','V' ]

#networkCMap = 'tab20'
networkCMap = categorical_cmap(8,3,cmap='Set1')
networkColors = cm.get_cmap('Set1', lut=len(tags))
networkBottomColor = 'k'
nodeSize = 50
edgeWidth = 0.5

cliqueTopColor = 'r'
nonCliqueColor = 'k'
nonCliqueAlpha = 0.5

cliqueHeatmapCmap = 'cool'
cliqueHetamapTopTicks = [0,9,19,29,39,49,59]
cliqueHeatmapBottomTicks = range(0,15)
cliqueHeatmapTopLabels = ['1','10','20','30','40','50','60']
cliqueHeatmapBottomLabels = [str(tick+1) for tick in cliqueHeatmapBottomTicks]

identityLineStyle = '--'
scatterMarker = '*'

dataColor = 'r'
dataMarker = '*'
dataLineWidth = 3.
randomColor = 'k'
randomMarker = '.'
randomAlpha = 0.2




# save paths
savePathBase = '/media/onerva/KINGSTON/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/'
degreeSaveName = 'degree-distributions'
degreeNodeDictionarySaveName = 'companies-per-degree'
networkSaveName = 'network'
cliqueSaveName = 'network-cliques'
cliqueHeatmapSaveName = 'clique-heatmap'
diversitySaveName = 'richness-vs-diversity'
comparisonVsRandomSaveName = 'comparison-vs-random'
relativeDiversitySaveName = 'relative-diversity'
diversityVsBottomIndexSaveName = 'diversity-vs-n-events'
diversityVsTopIndexSaveName = 'diveristy-vs-n-companies'