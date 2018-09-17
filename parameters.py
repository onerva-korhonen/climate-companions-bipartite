#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 16:58:12 2018

@author: onerva

This file contains parameters used for the bipartite analysis of the Climate 
Companions organization.
"""
from matplotlib import cm

# input
companyInputPath = '/media/onerva/KINGSTON/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Members_alias.csv'
eventInputPath = '/media/onerva/KINGSTON/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Events_alias.csv'
linkInputPath = '/media/onerva/KINGSTON/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/Linkit_2016.csv'
companyColumnNames = ['Alias:', 'Member:']
eventColumnNames = ['Alias:', 'Event:', 'Time:','Other information:']
linkColumnNames = ['Event:', 'Participant:']

# distributions and binning
nDegreeBins = 50

# visualization
cmap = 'cool'
colors = cm.get_cmap(cmap, lut=2)
topColor = colors(0)
bottomColor = colors(1)

# save paths
savePathBase = '/media/onerva/KINGSTON/aallon-tyokoneelta/lappari/misc_projects/ilmastokumppanit/'
degreeSaveName = 'degree-distributions.pdf'