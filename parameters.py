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
companyInputPath = '/media/onerva/KINGSTON/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Members_alias.csv'
eventInputPath = '/media/onerva/KINGSTON/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Events_alias.csv'
linkInputPath = '/media/onerva/KINGSTON/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Linkit_2016.csv'
companyColumnNames = ['Alias:', 'Member:']
eventColumnNames = ['Alias:', 'Event:', 'Time:','Other information:']
linkColumnNames = ['Event:', 'Participant:']

# distributions and binning
nDegreeBins = 50
cliqueHeatmapBottomBins = np.arange(0.5,5.5,1)
cliqueHeatmapTopBins = np.arange(0.5,51.5,1)
nRichnessBins = 5

# comparison against random

nRandomIterations = 1000
nRandomBins = 20

# visualization
cmap = 'cool'
colors = cm.get_cmap(cmap, lut=2)
topColor = colors(0)
bottomColor = colors(1)

tags = ['R','In','Se','CB','Me','SD','L','Con','N','E','T','CT','CS','S','HKI','C','ENGO','H','II','F']
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
cliqueHetamapTopTicks = [0,9,19,29,39,49]
cliqueHeatmapBottomTicks = [0,1,2,3]
cliqueHeatmapTopLabels = ['1','10','20','30','40','50']
cliqueHeatmapBottomLabels = ['1','2','3','4']

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
degreeSaveName = 'degree-distributions.pdf'
networkSaveName = 'network.pdf'
cliqueSaveName = 'network-cliques'
cliqueHeatmapSaveName = 'clique-heatmap.pdf'
diversitySaveName = 'richness-vs-diversity.pdf'
comparisonVsRandomSaveName = 'comparison-vs-random'
relativeDiversitySaveName = 'relative-diversity.pdf'