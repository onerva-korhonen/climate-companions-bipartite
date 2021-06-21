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
    A function for generating a list of link input paths. They're tried
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
# NOTE on link input paths: Link input paths are expected to follow a shared structure:
# linkInputStem + year + linkInputExtension
# If some of the paths don't follow this structure, they can be given using the manualLinkInputPaths
# that contains an element for each input path; for those paths that follow the shared structure, the elements
# can be arbitrary (e.g. ''). The staticNetworkInputPath (to the full network) should be given separately since there are some 
# analysis that are performed only for the static network. However, to get all analysis performed on both the full network
# and networks obtained in time windows, the static network path should be included in the linkInputPaths as well (probably
# it need to be given through manualLinkInputPaths)
companyInputPath = '/media/onerva/0012-D687/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin_uusi/Member_alias.csv'
eventInputPath = '/media/onerva/0012-D687/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin_uusi/Events_alias_uusi.csv'
manualLinkInputPaths = ['','','','','/media/onerva/0012-D687/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin_uusi/Links_all.csv']
years = ['2011_2012','2013_2014','2015_2016','2017_2018','2011_2018'] # the last 4 characters of the year should form an integer
linkInputStem = '/media/onerva/0012-D687/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin_uusi/Links_'
linkInputExtension = '_uusi.csv'
linkInputPaths = getLinkInputPaths(years, linkInputStem, linkInputExtension, manualLinkInputPaths)

staticNetworkInputPath = '/media/onerva/0012-D687/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Analyysiin_uusi/Links_all.csv'
                  
companyColumnNames = ['Alias:', 'Member:','Index:','Max degree:','Index slope:','Mask:', 'Joined on:']
eventColumnNames = ['Alias:', 'Event:']
linkColumnNames = ['Event:', 'Participant:','Weight:']

indexKey = 'Index:'
degreeNormalizationKey = 'Max degree:'
indexChangeKey = 'Index slope:'
maskKey = 'Mask:'
joiningYearKey = 'Joined on:'

csvSeparator = ';'

# distributions and binning
nTopDegreeBins = 10
nBottomDegreeBins = 20
cliqueHeatmapBottomBins = np.arange(0.5,16.5,1)
cliqueHeatmapTopBins = np.arange(0.5,68.5,1)
nRichnessBins = 5
nIndexBins = 20

lowDegreePercentile = 25
highDegreePercentile = 75

separateClasses = True # should top degree distribution be drawn separately for different membership classes?

analyzeZeroDegreeFields = True # should the field histogram be plotted separately for zero-degree nodes?

normalizeDegreeInScatter = True # should degree be normalized by its theoretical max value in the degree-index scatter?

# comparison against random
nRandomIterations = 1000
nRandomBins = 20
ignoreNonMembers = True

starnessXLims = (0,0.6)
starnessYLims = (0,0.3)
richnessXLims = (3,14)
richnessYLims = (0,0.25)
relativeDivXLims = (0.6,0.95)
relativeDivYLims = (0,0.18)

fieldMeanDegreesSigLimit = 0.01

# visualization
tags = ['C','CM','Co','E','I','II','LT','S','SD','RE','NGO','PS','AG','FPS','OC','HKI_1','HKI_2']
membershipClasses = ['BM','OM','NM','HKI']
nonMemberClasses = ['NM']
nonMemberClassesForStarness = ['NM','HKI']
nodesToExcludeFromScatter = ['OM_PS3','HKI_1','HKI_2']
# This is for excluding some outlier nodes from the degree - index scatter (OM_PS3 is an error and HKI_1 and HKI_2 instances of the city of Helsinki)
nodesToExcludeFromDegrees = ['HKI_1','HKI_2'] 
# We use this to exclude the instances of the City of Helsinki when calculating the average number of participants per event and event per participants

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
scatterMarker = '.'

classMarkers = ['.','^','','o']
markerAlpha = 0.7

dataColor = 'r'
dataMarker = '*'
dataLineWidth = 3.
randomColor = 'k'
randomMarker = '.'
randomAlpha = 0.2

histWidth = 0.75
fieldHistWidth = 0.2

richnessLineStyle = '-'
diversityLineStyle = '--'

indexPercentile = 0.5
indexPercentileLineStyle = '--'
indexPercentileColor = 'k'
indexPercentileAlpha = 0.5


# save paths
savePathBase = '/media/onerva/0012-D687/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/tulokset_apr_2021/'
degreeSaveName = 'degree-distributions'
degreeHistogramSaveName = '/degree_histograms'
degreeNodeDictionarySaveName = 'companies-per-degree'
networkSaveName = 'network'
cliquesSaveName = 'network-cliques'
cliqueHeatmapSaveName = 'clique-heatmap'
diversitySaveName = 'richness-vs-diversity'
comparisonVsRandomSaveName = 'comparison-vs-random'
relativeDiversitySaveName = 'relative-diversity'
diversityVsBottomIndexSaveName = 'diversity-vs-n-events'
diversityVsTopIndexSaveName = 'diveristy-vs-n-companies'
fieldHistogramClassesSaveName = 'field-histogram-classes'
fieldHistogramSaveName = 'field_histogram'
degreeIndexScatterSaveName = 'degree_index_scatter'
degreeIndexHeatmapSaveName = 'degree_index_heatmap'
sortedCliquesSaveName = 'cliques_sorted_by_events'