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

cfg['topColor'] = pms.topColor
cfg['bottomColor'] = pms.bottomColor
cfg['networkColors'] = pms.networkColors
cfg['networkBottomColor'] = pms.networkBottomColor
cfg['nodeSize'] = pms.nodeSize
cfg['edgeWidth'] = pms.edgeWidth

cfg['savePathBase'] = pms.savePathBase
cfg['degreeSaveName'] = pms.degreeSaveName
cfg['networkSaveName'] = pms.networkSaveName



bnet = functions.createBipartite(cfg)
bnet = functions.pruneBipartite(bnet)
for node, d in bnet.nodes(data=True):
    print node, d['bipartite']
density = functions.getDensity(bnet)
print 'Density: ' + str(density)
#functions.getDegreeDistributions(bnet, cfg)
functions.drawNetwork(bnet,cfg)
