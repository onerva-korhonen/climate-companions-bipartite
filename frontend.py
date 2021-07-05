#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 16:57:52 2018

@author: onerva

A frontend script for analysing the Climate Companions organization as a bipartite
network
"""
import parameters as pms
import functions

import numpy as np

years = pms.years
linkInputPaths = pms.linkInputPaths

topNodes = []
topNodesWithoutNonmembers = [] # this is needed for calculating Jaccards without non-members

densities = []
nCliques = []
starnesses = []
meanRichnesses = []
meanDiversities = []
meanRelativeDiversities = []

numbersOnly = False

if numbersOnly:
    cfg = {}

    cfg['companyInputPath'] = pms.companyInputPath
    cfg['eventInputPath'] = pms.eventInputPath
    cfg['linkInputPath'] = pms.staticNetworkInputPath
    cfg['companyColumnNames'] = pms.companyColumnNames
    cfg['eventColumnNames'] = pms.eventColumnNames
    cfg['linkColumnNames'] = pms.linkColumnNames
    cfg['tags'] = pms.tags
    cfg['classes'] = pms.membershipClasses
    cfg['csvSeparator'] = pms.csvSeparator
    
    cfg['ignoreNonMembers'] = pms.ignoreNonMembers
    cfg['nonMemberClasses'] = pms.nonMemberClasses
    nodesToExcludeFromDegrees = pms.nodesToExcludeFromDegrees
    
    cfg['topColor'] = pms.topColor
    cfg['networkColors'] = pms.networkColors
    cfg['networkBottomColor'] = pms.networkBottomColor
    cfg['nodeShapes'] = pms.nodeShapes
    cfg['bottomShape'] = pms.bottomShape

    for year, linkInputPath in zip(years,linkInputPaths):
        
        cfg['linkInputPath'] = linkInputPath
        
        cfg['savePathBase'] = pms.savePathBase
        cfg['degreeSaveName'] = pms.degreeSaveName + '_' + year + '.pdf'
        cfg['degreeHistogramSaveName'] = pms.degreeHistogramSaveName + '_' + year + '.pdf'
        cfg['degreeNodeDictionarySaveName'] = pms.degreeNodeDictionarySaveName + '_' + year + '.csv'
        cfg['networkSaveName'] = pms.networkSaveName + '_' + year + '.pdf'
        cfg['cliquesSaveName'] = pms.cliquesSaveName + '_' + year
        cfg['cliqueHeatmapSaveName'] = pms.cliqueHeatmapSaveName + '_' + year + '.pdf'
        cfg['diversitySaveName'] = pms.diversitySaveName + '_' + year + '.pdf' 
        cfg['comparisonVsRandomSaveName'] = pms.comparisonVsRandomSaveName + '_' + year
        cfg['relativeDiversitySaveName'] = pms.relativeDiversitySaveName + '_' + year
        cfg['diversityVsBottomIndexSaveName'] = pms.diversityVsBottomIndexSaveName + '_' + year + '.pdf'
        cfg['diversityVsTopIndexSaveName'] = pms.diversityVsTopIndexSaveName + '_' + year + '.pdf'
        cfg['fieldHistogramSaveName'] = pms.fieldHistogramSaveName + '_' + year + '.pdf'
    
        bnet = functions.createBipartite(cfg)
        
        # Counting the number of members and non-members participating in events in this window
        assert year[-4::].isdigit(),'Please give the years in a form where the last 4 characters form an integer (e.g. "2011_2012")'
        endYear = int(year[-4::])
        top,_ = functions.getTopAndBottom(bnet)
        tlist = list(top)
        joiningYears = [bnet.nodes(data=True)[node][pms.joiningYearKey] for node in tlist]
        nodeClasses = [bnet.nodes(data=True)[node]['class'] for node in tlist]
        degrees = [bnet.degree(node) for node in tlist]
        mclasses = list(cfg['classes'])
        for nonMemberClass in cfg['nonMemberClasses']:
            mclasses.remove(nonMemberClass)
        for mclass in mclasses:
            partCount = 0
            nonpartCount = 0
            for joiningYear,nodeClass,degree in zip(joiningYears,nodeClasses,degrees):
                if nodeClass == mclass:
                    if degree > 0:
                        partCount += 1
                    elif joiningYear <= endYear:
                        nonpartCount += 1
            print 'In year(s) ' + year + ', in class ' + mclass + ' ' + str(partCount) + ' participants, ' + str(nonpartCount) + ' non-participants'
        partCount = 0
        for joiningYear,nodeClass,degree in zip(joiningYears,nodeClasses,degrees):
            if nodeClass in cfg['nonMemberClasses'] and degree > 0:
                partCount += 1
        nonMemberStr = ', '.join(cfg['nonMemberClasses'])
        print 'In year(s) ' + year + ', in classes ' + nonMemberStr + ' ' + str(partCount) + ' participants'
                
        degreeDict = functions.getDegreeNodeDictionary(bnet,cfg)
        bnet, nZeroDegree = functions.pruneBipartite(bnet)
        
        density = functions.getDensity(bnet)
        densityWithoutNonMembers = functions.getDensity(bnet,excludeNonMembers=True,nonMemberClasses=cfg['nonMemberClasses'])
        print 'Density in year(s) ' + year + ': ' + str(density) + ', ' + str(densityWithoutNonMembers) + ' excluding non-members'
        densities.append(density)
     
        top,bottom = functions.getTopAndBottom(bnet)
        topNodes.append(top)
        
        print str(len(bottom)) + ' events organized in year(s) ' + year
        
        cnet = bnet.copy() # calculating the events/participant and participants/event requires excluding some nodes, so let's use a copy of the network
        for node in top:
            if node in nodesToExcludeFromDegrees: # first removing only the nodes particularly listed for being removed (in our case, the instances of the city of Helsinki that participate in all events and articially increase bottom degree)
                cnet.remove_node(node)
        
        ctop, cbottom = functions.getTopAndBottom(cnet)
        bottomDegrees = [cnet.degree(node) for node in cbottom]
        print 'In year(s) ' + year + ', on average ' + str(np.mean(bottomDegrees)) + ' participants per event, min ' + str(min(bottomDegrees)) + ', max ' + str(max(bottomDegrees))
        
        for node in ctop:
            if cnet.nodes(data=True)[node]['class'] in cfg['nonMemberClasses']: # now, removing all non-member nodes (that often participate only one event, thus artificially lowering mean events/participant)
                cnet.remove_node(node)
        
        ctop,_ = functions.getTopAndBottom(cnet)
        topDegrees = [cnet.degree(node) for node in ctop]
        print 'In year(s) ' + year + ', on average ' + str(np.mean(topDegrees)) + ' events per participant, min ' + str(min(topDegrees)) + ', max ' + str(max(topDegrees))
        
        ctopNodes = list(ctop) + nodesToExcludeFromDegrees # before calculating Jaccards without non-members, let's add back the member nodes that were excluded from degree calculations (instances of the city of Helsinki)
        topNodesWithoutNonmembers.append(ctopNodes)
        
        cliqueInfo = {'topNodes':top,'bottomNodes':bottom,'topIndex':len(top),'bottomIndex':len(bottom),'isBridge':False,'isStar':False}
    
        richness, cliqueDiversity, count, majorField = functions.getCliqueFieldDiversity(bnet,cliqueInfo)
    
        print 'Diversity of full network in year(s) ' + year + ': richness: ' + str(richness) + ', effective diversity: ' + str(cliqueDiversity) + ', relative diversity: ' + str(cliqueDiversity/richness)
        print str(nZeroDegree) + ' zero-degree nodes found'
    
    jaccards = []
    jaccardsWithoutNonmembers = []
        
    for i in range(len(years)-1):
        jaccard = functions.getJaccardIndex(topNodes[i],topNodes[i+1])
        jaccards.append(jaccard)
        jaccardWithoutNonmembers = functions.getJaccardIndex(topNodesWithoutNonmembers[i],topNodesWithoutNonmembers[i+1])
        jaccardsWithoutNonmembers.append(jaccardWithoutNonmembers)
        print 'Jaccard index, years ' + years[i] + ', ' + years[i+1] + ': ' + str(jaccard) + ', without non-members ' + str(jaccardWithoutNonmembers)
        
else:   
    cfg = {}
    
    cfg['companyInputPath'] = pms.companyInputPath
    cfg['eventInputPath'] = pms.eventInputPath
    cfg['linkInputPath'] = pms.staticNetworkInputPath
    cfg['companyColumnNames'] = pms.companyColumnNames
    cfg['eventColumnNames'] = pms.eventColumnNames
    cfg['linkColumnNames'] = pms.linkColumnNames
    cfg['tags'] = pms.tags
    cfg['classes'] = pms.membershipClasses
    cfg['csvSeparator'] = pms.csvSeparator
    cfg['indexKey'] = pms.indexKey
    cfg['degreeNormalizationKey'] = pms.degreeNormalizationKey
    cfg['indexChangeKey'] = pms.indexChangeKey
    cfg['maskKey'] = pms.maskKey
    cfg['nodesToExcludeFromScatter'] = pms.nodesToExcludeFromScatter
    
    cfg['ignoreNonMembers'] = pms.ignoreNonMembers
    cfg['nonMemberClasses'] = pms.nonMemberClasses
    
    cfg['nTopDegreeBins'] = pms.nTopDegreeBins
    cfg['nBottomDegreeBins'] = pms.nBottomDegreeBins
    cfg['cliqueHeatmapTopBins'] = pms.cliqueHeatmapTopBins
    cfg['cliqueHeatmapBottomBins'] = pms.cliqueHeatmapBottomBins
    cfg['nRichnessBins'] = pms.nRichnessBins
    cfg['nIndexBins'] = pms.nIndexBins
    
    cfg['separateClasses'] = pms.separateClasses
    cfg['analyzeZeroDegreeFields'] = pms.analyzeZeroDegreeFields
    cfg['normalizeDegreeInScatter'] = pms.normalizeDegreeInScatter
    
    cfg['nRandomIterations'] = pms.nRandomIterations
    cfg['nRandomBins'] = pms.nRandomBins
    cfg['fieldMeanDegreesSigLimit'] = pms.fieldMeanDegreesSigLimit
    
    cfg['topColor'] = pms.topColor
    cfg['bottomColor'] = pms.bottomColor
    cfg['classColors'] = pms.classColors
    cfg['networkColors'] = pms.networkColors
    cfg['nonUniqueColorTags'] = pms.nonUniqueColorTags
    cfg['networkNonUniqueColor'] = pms.networkNonUniqueColor
    cfg['networkBottomColor'] = pms.networkBottomColor
    cfg['nodeSize'] = pms.nodeSize
    cfg['nodeShapes'] = pms.nodeShapes
    cfg['bottomShape'] = pms.bottomShape
    cfg['edgeWidth'] = pms.edgeWidth
    cfg['edgeAlpha'] = pms.edgeAlpha
    cfg['cliqueTopColor'] = pms.cliqueTopColor
    cfg['nonCliqueColor'] = pms.nonCliqueColor
    cfg['nonCliqueAlpha'] = pms.nonCliqueAlpha
    cfg['cliqueHeatmapCmap'] = pms.cliqueHeatmapCmap
    cfg['cliqueHeatmapTopTicks'] = pms.cliqueHetamapTopTicks
    cfg['cliqueHeatmapBottomTicks'] = pms.cliqueHeatmapBottomTicks
    cfg['cliqueHeatmapTopLabels'] = pms.cliqueHeatmapTopLabels
    cfg['cliqueHeatmapBottomLabels'] = pms.cliqueHeatmapBottomLabels
    cfg['identityLineStyle'] = pms.identityLineStyle
    cfg['scatterMarker'] = pms.scatterMarker
    cfg['classMarkers'] = pms.classMarkers
    cfg['markerAlpha'] = pms.markerAlpha
    cfg['randomColor'] = pms.randomColor
    cfg['randomMarker'] = pms.randomMarker
    cfg['randomAlpha'] = pms.randomAlpha
    cfg['dataColor'] = pms.dataColor
    cfg['dataMarker'] = pms.dataMarker
    cfg['dataLineWidth'] = pms.dataLineWidth
    cfg['histWidth'] = pms.histWidth
    cfg['fieldHistWidth'] = pms.fieldHistWidth
    cfg['starnessXLims'] = pms.starnessXLims
    cfg['starnessYLims'] = pms.starnessYLims
    cfg['richnessXLims'] = pms.richnessXLims
    cfg['richnessYLims'] = pms.richnessYLims
    cfg['relativeDivXLims'] = pms.relativeDivXLims
    cfg['relativeDivYLims'] = pms.relativeDivYLims
    cfg['richnessLineStyle'] = pms.richnessLineStyle
    cfg['diversityLineStyle'] = pms.diversityLineStyle
    cfg['indexPercentile'] = pms.indexPercentile
    cfg['indexPercentileLineStyle'] = pms.indexPercentileLineStyle
    cfg['indexPercentileColor'] = pms.indexPercentileColor
    cfg['indexPercentileAlpha'] = pms.indexPercentileAlpha
    
    cfg['savePathBase'] = pms.savePathBase
    cfg['degreeIndexScatterSaveName'] = pms.degreeIndexScatterSaveName + '_all_degree_normalized.pdf'
    cfg['degreeIndexHeatmapSaveName'] = pms.degreeIndexHeatmapSaveName + '_all_degree_normalized.pdf'
    
    bnet = functions.createBipartite(cfg)
    functions.createDegreeIndexScatter(bnet, cfg)
    functions.createDegreeIndexHeatmap(bnet, cfg)
        
    for year, linkInputPath in zip(years,linkInputPaths):
        cfg['linkInputPath'] = linkInputPath
        
        cfg['degreeSaveName'] = pms.degreeSaveName + '_' + year + '.pdf'
        cfg['degreeHistogramSaveName'] = pms.degreeHistogramSaveName + '_' + year + '.pdf'
        cfg['networkSaveName'] = pms.networkSaveName + '_' + year + '.pdf'
        cfg['cliquesSaveName'] = pms.cliquesSaveName + '_' + year
        cfg['cliqueHeatmapSaveName'] = pms.cliqueHeatmapSaveName + '_' + year + '.pdf'
        cfg['diversitySaveName'] = pms.diversitySaveName + '_' + year + '.pdf' 
        cfg['comparisonVsRandomSaveName'] = pms.comparisonVsRandomSaveName + year
        cfg['relativeDiversitySaveName'] = pms.relativeDiversitySaveName + year + '.pdf'
        cfg['diversityVsBottomIndexSaveName'] = pms.diversityVsBottomIndexSaveName + '_' + year + '.pdf'
        cfg['diversityVsTopIndexSaveName'] = pms.diversityVsTopIndexSaveName + '_' + year + '.pdf'
        cfg['fieldHistogramClassesSaveName'] = pms.fieldHistogramClassesSaveName + '_' + year + '.pdf'
        cfg['fieldHistogramSaveName'] = pms.fieldHistogramSaveName + '_' + year + '.pdf'
        sortedCliquesSaveName = pms.savePathBase + pms.sortedCliquesSaveName + '_' + year + '.txt'
        
        cfg['skipNonMembersInVisualization'] = False
        cfg['nonMemberClasses'] = pms.nonMemberClasses
        
        bnet = functions.createBipartite(cfg)
        functions.getDegreeHistogram(bnet, cfg)
        functions.getFieldHistogram(bnet,cfg)
        bnet, nZeroDegree = functions.pruneBipartite(bnet)
    #    print str(nZeroDegree) + ' zero degree nodes removed'
        lowPercentileNodes,highPercentileNodes = functions.findTopNodesInPercentile(bnet,pms.lowDegreePercentile,pms.highDegreePercentile,cfg)
        density = functions.getDensity(bnet)
        print 'Density: ' + str(density)
        densities.append(density)
    #    #functions.getDegreeDistributions(bnet, cfg)
        functions.drawNetwork(bnet,cfg)
        cfg['skipNonMembersInVisualization'] = True
        cfg['networkSaveName'] = pms.networkSaveName + '_' + year + '_without_nonmembers.pdf'
        functions.drawNetwork(bnet,cfg) # repeat network drawing without non-member nodes
        cliques, cliqueInfo = functions.findBicliques(bnet)
        nCliques.append(len(cliques))
        functions.saveBicliques(cliqueInfo,sortedCliquesSaveName)
        richnesses, diversities, counts, majorFields = functions.getCliqueFieldDiversityWrapper(bnet,cliqueInfo)
        meanRichnesses.append(np.mean(np.array(richnesses)))
        meanDiversities.append(np.mean(np.array(diversities)))
        relativeDiversity = np.array(diversities)/np.array(richnesses)
        meanRelativeDiversities.append(np.mean(relativeDiversity))
        functions.plotRichnessVsDiversity(richnesses,diversities,cfg) # Note: for some reason, this command does not work nicely in Spyder (only one plot is saved). It works as it should when run from the terminal.
        functions.plotRelativeDiversity(cliques,richnesses,diversities,cfg)
        functions.plotDiversityVsIndices(cliqueInfo,richnesses,diversities,cfg)
        functions.createCliqueIndexHeatmap(cliqueInfo, cfg)        
        
        measures = {'richness':richnesses,'diversity':diversities, 'cliques':cliques}
        
        functions.compareAgainstRandom(bnet,cfg,measures)
        
        cfg['nonMemberClasses'] = pms.nonMemberClassesForStarness
        # for analyzing starness, let's remove from stars the nodes that participate also in other cliques
        cliques, cliqueInfo = functions.pruneStars(bnet,cliques,cliqueInfo,ignoreNonMembers=cfg['ignoreNonMembers'],nonMemberClasses=cfg['nonMemberClasses'])
        functions.visualizeBicliques(bnet,cliqueInfo,cfg)
    #    #functions.createCliqueIndexHeatmap(cliqueInfo, cfg)
        starness = functions.getStarness(bnet,cliqueInfo,ignoreNonMembers=cfg['ignoreNonMembers'],nonMemberClasses=cfg['nonMemberClasses'])
        print 'Starness: ' + str(starness)
        starnesses.append(starness)
        
        measures = {'starness':starness}
        
        functions.compareAgainstRandom(bnet,cfg,measures)
        print year + ' OK'
        
    print 'Densities:'
    print densities
    print 'Starnesses:'
    print starnesses
    print 'Mean richnesses:'
    print meanRichnesses
    print 'Mean diversities:'
    print meanDiversities
    print 'Mean relative diversities:'
    print meanRelativeDiversities