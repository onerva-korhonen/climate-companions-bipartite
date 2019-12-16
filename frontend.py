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

densities = []
nCliques = []
starnesses = []
meanRichnesses = []
meanDiversities = []
meanRelativeDiversities = []

numbersOnly = True

if numbersOnly:
    cfg = {}

    cfg['companyInputPath'] = pms.companyInputPath
    cfg['eventInputPath'] = pms.eventInputPath
    cfg['linkInputPath'] = pms.staticNetworkInputPath
    cfg['companyColumnNames'] = pms.companyColumnNames
    cfg['eventColumnNames'] = pms.eventColumnNames
    cfg['linkColumnNames'] = pms.linkColumnNames
    cfg['tags'] = pms.tags
    
    cfg['topColor'] = pms.topColor
    cfg['networkColors'] = pms.networkColors
    cfg['networkBottomColor'] = pms.networkBottomColor
     
    bnet = functions.createBipartite(cfg)
    bnet, nZeroDegree = functions.pruneBipartite(bnet)
     
    top,bottom = functions.getTopAndBottom(bnet)
    cliqueInfo = {'topNodes':top,'bottomNodes':bottom,'topIndex':len(top),'bottomIndex':len(bottom),'isBridge':False,'isStar':False}
    
    richness, cliqueDiversity, count, majorField = functions.getCliqueFieldDiversity(bnet,cliqueInfo)
    
    print 'Diversity of full network: richness: ' + str(richness) + ', effective diversity: ' + str(cliqueDiversity) + ', relative diversity: ' + str(cliqueDiversity/richness)
    print str(nZeroDegree) + ' zero-degree nodes found'
    
    topNodes = []
    
    for linkInputPath in linkInputPaths:
        cfg['linkInputPath'] = linkInputPath
        bnet = functions.createBipartite(cfg)
        bnet, nZeroDegree = functions.pruneBipartite(bnet)
        
        top,bottom = functions.getTopAndBottom(bnet)
        topNodes.append(top)
    
    for i,yeari in enumerate(years):
        for j,yearj in enumerate(years):
            jaccard = functions.getJaccardIndex(topNodes[i],topNodes[j])
            print 'Jaccard index, years ' + yeari + ', ' + yearj + ': ' + str(jaccard)
    
    
else:   
    for year, linkInputPath in zip(years,linkInputPaths):
        cfg = {}
    
        cfg['companyInputPath'] = pms.companyInputPath
        cfg['eventInputPath'] = pms.eventInputPath
        cfg['linkInputPath'] = linkInputPath
        cfg['companyColumnNames'] = pms.companyColumnNames
        cfg['eventColumnNames'] = pms.eventColumnNames
        cfg['linkColumnNames'] = pms.linkColumnNames
        cfg['tags'] = pms.tags
        
        cfg['nDegreeBins'] = pms.nDegreeBins
        cfg['cliqueHeatmapTopBins'] = pms.cliqueHeatmapTopBins
        cfg['cliqueHeatmapBottomBins'] = pms.cliqueHeatmapBottomBins
        cfg['nRichnessBins'] = pms.nRichnessBins
        
        cfg['nRandomIterations'] = pms.nRandomIterations
        cfg['nRandomBins'] = pms.nRandomBins
        
        cfg['topColor'] = pms.topColor
        cfg['bottomColor'] = pms.bottomColor
        cfg['networkColors'] = pms.networkColors
        cfg['networkBottomColor'] = pms.networkBottomColor
        cfg['nodeSize'] = pms.nodeSize
        cfg['edgeWidth'] = pms.edgeWidth
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
        cfg['randomColor'] = pms.randomColor
        cfg['randomMarker'] = pms.randomMarker
        cfg['randomAlpha'] = pms.randomAlpha
        cfg['dataColor'] = pms.dataColor
        cfg['dataMarker'] = pms.dataMarker
        cfg['dataLineWidth'] = pms.dataLineWidth
        
        cfg['savePathBase'] = pms.savePathBase
        cfg['degreeSaveName'] = pms.degreeSaveName + '_' + year + '.pdf'
        cfg['networkSaveName'] = pms.networkSaveName + '_' + year + '.pdf'
        cfg['cliquesSaveName'] = pms.cliqueSaveName + '_' + year
        cfg['cliqueHeatmapSaveName'] = pms.cliqueHeatmapSaveName + '_' + year + '.pdf'
        cfg['diversitySaveName'] = pms.diversitySaveName + '_' + year + '.pdf' 
        cfg['comparisonVsRandomSaveName'] = pms.comparisonVsRandomSaveName + year
        cfg['relativeDiversitySaveName'] = pms.relativeDiversitySaveName + year
        
        bnet = functions.createBipartite(cfg)
        bnet, nZeroDegree = functions.pruneBipartite(bnet)
    #    print str(nZeroDegree) + ' zero degree nodes removed'
        density = functions.getDensity(bnet)
        print 'Density: ' + str(density)
        densities.append(density)
    #    #functions.getDegreeDistributions(bnet, cfg)
        functions.drawNetwork(bnet,cfg)
        cliques, cliqueInfo = functions.findBicliques(bnet)
        nCliques.append(len(cliques))
        richnesses, diversities, counts, majorFields = functions.getCliqueFieldDiversityWrapper(bnet,cliqueInfo)
        meanRichnesses.append(np.mean(np.array(richnesses)))
        meanDiversities.append(np.mean(np.array(diversities)))
        relativeDiversity = np.array(diversities)/np.array(richnesses)
        meanRelativeDiversities.append(np.mean(relativeDiversity))
        functions.plotRichnessVsDiversity(richnesses,diversities,cfg) # Note: for some reason, this command does not work nicely in Spyder (only one plot is saved). It works as it should when run from the terminal.
        functions.plotRelativeDiversity(cliques,richnesses,diversities,cfg)
        
        measures = {'richness':richnesses,'diversity':diversities, 'cliques':cliques}
        
        functions.compareAgainstRandom(bnet,cfg,measures)
        
        # for analyzing starness, let's remove from stars the nodes that participate also in other cliques
        cliques, cliqueInfo = functions.pruneStars(bnet,cliques,cliqueInfo)
        functions.visualizeBicliques(bnet,cliqueInfo,cfg)
    #    #functions.createCliqueIndexHeatmap(cliqueInfo, cfg)
        starness = functions.getStarness(bnet,cliqueInfo)
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