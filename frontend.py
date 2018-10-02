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

cfg = {}

cfg['companyInputPath'] = pms.companyInputPath
cfg['eventInputPath'] = pms.eventInputPath
cfg['linkInputPath'] = pms.linkInputPath
cfg['companyColumnNames'] = pms.companyColumnNames
cfg['eventColumnNames'] = pms.eventColumnNames
cfg['linkColumnNames'] = pms.linkColumnNames
cfg['tags'] = pms.tags

cfg['nDegreeBins'] = pms.nDegreeBins
cfg['cliqueHeatmapTopBins'] = pms.cliqueHeatmapTopBins
cfg['cliqueHeatmapBottomBins'] = pms.cliqueHeatmapBottomBins

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

cfg['savePathBase'] = pms.savePathBase
cfg['degreeSaveName'] = pms.degreeSaveName
cfg['networkSaveName'] = pms.networkSaveName
cfg['cliquesSaveName'] = pms.cliqueSaveName
cfg['cliqueHeatmapSaveName'] = pms.cliqueHeatmapSaveName

bnet = functions.createBipartite(cfg)
bnet = functions.pruneBipartite(bnet)
#for node, d in bnet.nodes(data=True):
#    print node, d['bipartite']
density = functions.getDensity(bnet)
print 'Density: ' + str(density)
#functions.getDegreeDistributions(bnet, cfg)
#functions.drawNetwork(bnet,cfg)
cliques, cliqueInfo = functions.findBicliques(bnet)
cliques, cliqueInfo = functions.pruneStars(bnet,cliques,cliqueInfo)
#functions.visualizeBicliques(bnet,cliqueInfo,cfg)
functions.createCliqueIndexHeatmap(cliqueInfo, cfg)