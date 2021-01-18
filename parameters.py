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
import os
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

def getLinkInputPaths(years, linkInputStem, extension='.csv', manualLinkInputPaths=[]):
    """
    A function for generating a list of link input pahts. They're tried
    to generate from a file name stem and year, but if the file corresponding
    to the generated name doesn't exist, the name is picked from the manually
    entered list instead.
    
    Parameters:
    ----------
    years: list of strs
        the years, the data of which should be used as input
    linkInputStem: str
        a file name stem that, combined with the year and extension, should
        give the actual file names
    extension: str
        the file extension
    manualLinkInputPaths: list of strs
        elements of manualLinkInputPaths should be file names. If a file 
        is not found at the generated path, the corresponding element
        of manualLinkInputPaths is used instead. For years, for which the 
        generated file name should work, the corresponding element of 
        manualLinkInputPaths can be arbitrary, e.g. ''. If the list is empty,
        generated file paths are returned without checking for file existance.
        
    Returns:
    --------
    linkInputPaths: list of strs
        the generated link input file paths
    """
    linkInputPaths = []
    if len(manualLinkInputPaths) == 0:
        for year in years:
            linkInputPath = linkInputStem + year + extension
            linkInputPaths.append(linkInputPath)
    else:
        for i, year in enumerate(years):
            linkInputPath = linkInputStem + year + extension
            if os.path.isfile(linkInputPath):
                linkInputPaths.append(linkInputPath)
            else:
                linkInputPaths.append(manualLinkInputPaths[i])
    return linkInputPaths

# input
companyInputPath = '/media/onerva/0012-D687/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin_uusi/Member_alias.csv'
eventInputPath = '/media/onerva/0012-D687/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin_uusi/Events_alias_uusi.csv'
manualLinkInputPaths = []
                  #['/media/onerva/0012-D687/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Links_all_271018.csv',
                  #'/media/onerva/0012-D687/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Links_2011.csv']
                 #['/media/onerva/0012-D687/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Links_2011.csv',
                  #'/media/onerva/0012-D687/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Links_2012.csv',
                  #'/media/onerva/0012-D687/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Links_2013.csv',
                  #'/media/onerva/0012-D687/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Links_2014.csv',
                  #'/media/onerva/0012-D687/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Links_2015.csv',
                  #'/media/onerva/0012-D687/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Links_2016.csv',
                  #'/media/onerva/0012-D687/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Links_2017.csv',
                  #'/media/onerva/0012-D687/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Links_2018.csv',
                  #'/media/onerva/0012-D687/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Links_2011_2012.csv',
                  #'/media/onerva/0012-D687/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Links_2013_2014.csv',
                  #'/media/onerva/0012-D687/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Links_2015_2016.csv',
                  #'/media/onerva/0012-D687/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Links_2017_2018.csv',
                  #'/media/onerva/0012-D687/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin/Links_all_271018.csv']
years = ['2015_2016']#['11-12','13-14','15-16','17-18','all']
linkInputStem = '/media/onerva/0012-D687/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin_uusi/Links_'
linkInputExtension = '_uusi.csv'
linkInputPaths = getLinkInputPaths(years, linkInputStem, linkInputExtension, manualLinkInputPaths)

staticNetworkInputPath = '/media/onerva/0012-D687/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin_uusi/Links_all_uusi.csv'
                  
companyColumnNames = ['Alias:', 'Member:']
eventColumnNames = ['Alias:', 'Event:']
linkColumnNames = ['Event:', 'Participant:','Weight:']

csvSeparator = ';'

# distributions and binning
nTopDegreeBins = 20
nBottomDegreeBins = 20
cliqueHeatmapBottomBins = np.arange(0.5,16.5,1)
cliqueHeatmapTopBins = np.arange(0.5,68.5,1)
nRichnessBins = 5

lowDegreePercentile = 25
highDegreePercentile = 75

separateClasses = True # should top degree distribution be drawn separately for different membership classes?

analyzeZeroDegreeFields = True # should the field histogram be plotted separately for zero-degree nodes?

# comparison against random
nRandomIterations = 1000
nRandomBins = 20
ignoreNonMembers = True

# visualization
tags = ['C','CM','Co','E','I','II','LT','S','SD','RE','NGO','PS','AG','FPS','OC','HKI_1','HKI_2']
membershipClasses = ['BM','OM','NM','HKI']
nonMemberClass = 'NM'

cmap = 'cool'
colors = cm.get_cmap(cmap, lut=len(membershipClasses)+2)
topColor = colors(0)
classColors = [colors(i) for i in range(1,len(membershipClasses)+1)]
bottomColor = colors(-1)

#networkCMap = 'tab20'
networkCMap = categorical_cmap(8,3,cmap='Set1')
networkColors = cm.get_cmap('Set1', lut=len(tags))
networkBottomColor = 'k'
nodeSize = 50
businessMemberNodeShape = 'o'
otherMemberNodeShape = 's'
nonMemberNodeShape = 'd'
HKINodeShape = '^'
nodeShapes = [businessMemberNodeShape, otherMemberNodeShape, nonMemberNodeShape, HKINodeShape]
bottomShape = 'o'
edgeWidth = 'weight'
edgeAlpha = 0.5

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

histWidth = 0.75
fieldHistWidth = 0.2


# save paths
savePathBase = '/media/onerva/0012-D687/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/tulokset_jan_2021/'
degreeSaveName = 'degree-distributions'
degreeHistogramSaveName = '/degree_histograms'
degreeNodeDictionarySaveName = 'companies-per-degree'
networkSaveName = 'network'
cliqueSaveName = 'network-cliques'
cliqueHeatmapSaveName = 'clique-heatmap'
diversitySaveName = 'richness-vs-diversity'
comparisonVsRandomSaveName = 'comparison-vs-random'
relativeDiversitySaveName = 'relative-diversity'
diversityVsBottomIndexSaveName = 'diversity-vs-n-events'
diversityVsTopIndexSaveName = 'diveristy-vs-n-companies'
fieldHistogramClassesSaveName = 'field-histogram-classes'
fieldHistogramSaveName = 'field_histogram'