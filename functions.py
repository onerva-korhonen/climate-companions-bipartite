# -*- coding: utf-8 -*-
"""
Created on Fri May  4 12:24:58 2018

@author: aokorhon

This file contains functions related to the Ilmastokumppanit bipartite project. To call, either write a 
frontend script (see XX for example) or call interactively from ipython.
"""
import pandas as pd
import networkx as nx
import numpy as np
import matplotlib.pylab as plt
import random
from networkx.algorithms import bipartite
from scipy.stats import binned_statistic, binned_statistic_2d
from itertools import combinations

# Functions for reading metadata

def readNodes(cfg):
    """
    Reads information about network nodes and their attributes. To this end, the information
    must be stored in a csv file with named columns. By default, column names are set for reading
    company information ('Alias' (anonymized alias used to refer to the
    node in the network), 'Name' (real name of the company), 'FoB' (field of business), and 
    'Membership class' (BM = business member, OM = other member, NM = non-member)).
    Other sets of column names can be given as parameters. Data of the first column will be used
    as node labels and must therefore be anonymized.
    
    Parameters:
    -----------
    cfg: dict, contains:
        inputPath: str, path to the input csv file
        columnNames: list of strs, column names in the csv file (default ['Alias', 'Name', 'FoB', 'Membership class'])
        
    Returns:
    --------
    nodeInfo: 2D list, first element contains node aliases and further elements additional attributes
    columnNames: list of strs, column names of the csv file (and therefore titles of the lists in nodeInfo)
    
    """
    inputPath = cfg['inputPath']
    if 'columnNames' in cfg.keys():
        columnNames= cfg['columnNames']
    else:
        columnNames = ['Alias', 'Name', 'FoB', 'Membership class']
    
    nodeData = pd.read_csv(inputPath)
    nodeInfo = []
    
    #import pdb; pdb.set_trace()
    
    for columnIndex, columnName in enumerate(columnNames):
        dlist = []
        data = nodeData[columnName]
        for i in data:
            dlist.append(i)
        nodeInfo.append(dlist)
        
    return nodeInfo, columnNames
    
def readLinks(cfg):
    """
    Reads the link information from a csv file. The default column names are ['Source', 'Target'] but other
    names can be given as parameters. Note that the identifiers of source and target nodes in the file must
    equal to nodes' aliases.
    
    Parameters:
    -----------
    cfg: dict, contains:
        inputPath: str, path to the input csv file
        
    Returns:
    --------
    links: list of (source,target) tuples
    """
    inputPath = cfg['inputPath']
    if 'columnNames' in cfg.keys():
        columnNames= cfg['columnNames']
    else:
        columnNames = ['Source', 'Target']
    
    linkData = pd.read_csv(inputPath)
    links = []
    
    sources = linkData[columnNames[0]]
    targets = linkData[columnNames[1]]
    
    for source, target in zip(sources,targets):
        links.append((source,target))
    
    return links

    
    
    
    
    
# Functions for network construction
    
def addNodes(cfg, net, bipartiteClass):
    """
    Based on the given metadata, adds a set of nodes to a bipartite graph.
    
    Parameters:
    -----------
    cfg: dict, contains:
        inputPath: str
            path to the metadata input csv file
        columnNames: list of strs 
            column names in the csv file
        net: nx.Graph()
            the network where to add nodes
        bipartite: int (0 or 1)
            bipartite node class (top or bottom)
        tags: list of strs
            the field of business tags
        classes: list of strs
            the possible membership classes
        topColors: bln
            shall the node color attributes be set according to the field of 
            business indicates in the alias? If False, a separately given 
            networkBottomColor is used instead.
        classShapes: bln
            shall the node shape attributes be set according to the membership class? If false, 'o' is used as the shape
            of all nodes.
        networkBottomColor: str
            the color for bottom (event) nodes
        nodeShapes: list of strs
            shapes associated with the different membership classes
        bottomShape: str
            the shape of bottom (event) nodes
    Returns:
    --------
        net: nx.Graph()
            the network with given nodes added
    """
    topColors = cfg['topColors']
    classShapes = cfg['classShapes']
    networkColors = cfg['networkColors']
    nodeShapes = cfg['nodeShapes']
    tags = cfg['tags']
    classes = cfg['classes']
    # Reading node information
    nodeInfo, columnNames = readNodes(cfg)
    aliases = nodeInfo[0] # first column is used as node alias
    attributes = nodeInfo[1::] # additional columns contain attributes
    if 'FoB' in columnNames:
        FoBs = attributes[columnNames.index('FoB')-1]
    if 'Membership class' in columnNames:
        membershipClasses = attributes[columnNames.index('Membership class')-1]
    
    # Adding nodes and updating their attributes
    for i, node in enumerate(aliases):
        net.add_node(node, bipartite=bipartiteClass)
        attr_dic = {}
        for columnName, attribute in zip(columnNames[1::], attributes):
            attr_dic[columnName] = attribute[i]
        if topColors:
            if 'FoB' in columnNames:
                attr_dic['nodeColor'] = networkColors(tags.index(FoBs[i]))
                attr_dic['tag'] = FoBs[i]
            else: # assuming that tag is included in the alias
                foundTags = []
                for tag in tags:
                    if tag in node:
                        foundTags.append(tag)
                foundLengths = [len(found) for found in foundTags]
                foundTag = foundTags[foundLengths.index(max(foundLengths))]
                attr_dic['nodeColor'] = networkColors(tags.index(foundTag))
                attr_dic['tag'] = foundTag
            if classShapes:
                if 'Membership class' in columnNames:
                    attr_dic['nodeShape'] = nodeShapes[classes.index(membershipClasses[i])]
                    attr_dic['class'] = membershipClasses[i]
                else: # assuming that membership class is included in the alias
                    foundClasses = []
                    for mclass in classes:
                        teststr = mclass + '_'
                        if teststr in node:
                            foundClasses.append(mclass)
                    foundLengths = [len(found) for found in foundClasses]
                    foundClass = foundClasses[foundLengths.index(max(foundLengths))]
                    attr_dic['nodeShape'] = nodeShapes[classes.index(foundClass)]
                    attr_dic['class'] = foundClass
            else:
                attr_dic['nodeShape'] = cfg['bottomShape']
                attr_dic['class'] = 'bottom (event)'
        elif classShapes and not topColors:
            if 'Membership class' in columnNames:
                attr_dic['nodeShape'] = nodeShapes[classes.index(membershipClasses[i])]
                attr_dic['class'] = membershipClasses[i]
            else: # assuming that membership class is included in the alias
                foundClasses = []
                for mclass in classes:
                    teststr = mclass + '_'
                    if teststr in node:
                        foundClasses.append(mclass)
                foundLengths = [len(found) for found in foundClasses]
                foundClass = foundClasses[foundLengths.index(max(foundLengths))]
                attr_dic['nodeShape'] = nodeShapes[classes.index(foundClass)]
                attr_dic['class'] = foundClass
                
            attr_dic['nodeColor'] = cfg['networkBottomColor']
            attr_dic['tag'] = 'bottom (event)'
        else:
            attr_dic['nodeColor'] = cfg['networkBottomColor']
            attr_dic['tag'] = 'bottom (event)'
        net.nodes[node].update(attr_dic)
        
    return net
    
def createBipartite(cfg):
    """
    Creates a bipartite network with two types of nodes (companies and events), node attributes (real names, fields of business, ...), 
    and links between different types of nodes.
    
    Parameters:
    -----------
    cfg: dict, contains:
        companyInputPath: str
            path to the company metadata input csv file
        eventInputPath: str
            path to the event metadata input csv file
        linkInputPath: str
            path to the link input csv file
        companyColumnNames: list of strs
            column names in the csv file containing the company information (default ['Alias', 'Name', 'FoB'])
        eventColumnNames: list of strs
            column names in the csv file containing the company information (default ['Alias', 'Name'])
        linkColumnNames: list of strs
            column names in the link csv file (default ['Source', 'Target'])
        topColors: bln
            shall the node color attributes be set according to the field of business indicated in the alias? If False, 
            a separately given networkBottomColor is used instead.
        classShapes: bln
            shall the node shape attributes be set according to the membership class? If false, a separately given
            bottomShape is used as the shape
            of all nodes.
        networkColors: matplotlib.cmap
            colors associated with different fields of business
        nodeShapes: list of strs
            shapes associated with the different membership classes
        tags: list of strs
            the possible field of business tags
        classes: list of strs
            the possible membership classes
        networkBottomColor: str
            the color for bottom (event) nodes
        bottomShape: str
            the shape of bottom (event) nodes
        
    Returns:
    --------
    bnet: networkx.Graph(), a bipartite
    """
    bnet = nx.Graph()
    
    # Reading company information, adding company nodes
    cfg['inputPath'] = cfg['companyInputPath']
    cfg['topColors'] = True
    cfg['classShapes'] = True
    if 'companyColumnNames' in cfg.keys():
        cfg['columnNames'] = cfg['companyColumnNames']
    else:
        cfg['columnNames'] = ['Alias', 'Name', 'FoB', 'Membership class']
    bnet = addNodes(cfg, bnet, bipartiteClass=0)   
    
    # Reading event information, adding event nodes
    cfg['inputPath'] = cfg['eventInputPath']
    cfg['topColors'] = False
    cfg['classShapes'] = False
    if 'eventColumnNames' in cfg.keys():
        cfg['columnNames'] = cfg['eventColumnNames']
    else:
        cfg['columnNames'] = ['Alias', 'Name']
    bnet = addNodes(cfg, bnet, bipartiteClass=1)
    
    # Reading link data, adding links
    cfg['inputPath'] = cfg['linkInputPath']
    if 'linkColumnNames' in cfg.keys():
        cfg['columnNames'] = cfg['linkColumnNames']
    else:
        cfg['columnNames'] = ['Source','Target']
    links = readLinks(cfg)
    for link in links:
        source = link[0]
        target = link[1]
        if not source in bnet.nodes or not target in bnet.nodes:
            print 'Error in link'
            print link
    bnet.add_edges_from(links)
    
    return bnet
    
def pruneBipartite(bnet):
    """
    Removes from a (bipartite) network all nodes with degree 0
    
    Parameters:
    ----------
    bnet: networkx.Graph(), a bipartite
    
    Returns:
    --------
    bnet: networkx.Graph(), the input network without nodes with degree 0
    nZeroDegree: int, number of the removed nodes
    """
    nodesToRemove = []
    
    for node in bnet.nodes():
        if bnet.degree(node) == 0:
            nodesToRemove.append(node)
    
    bnet.remove_nodes_from(nodesToRemove)
    nZeroDegree = len(nodesToRemove)
    
    return bnet, nZeroDegree
    
# Basic analysis
def getDegreeDistributions(bnet, cfg):
    """
    Calculates the top and bottom (companies and events) degree distributions of
    a bipartite network, plots the distributions, and saves the plot.
    
    Parameters:
    -----------
    bnet: networkx.Graph(), bipartite
    cfg: dict, contains:
        savePathBase: str, a base path (e.g. to a shared folder) for saving figures
        degreeSaveName: str, name of the file where to save the degree distribution plots
        nTopDegreeBins: int, number of bins used to calculate the top degree distribution
        nBottomDegreeBins: int, number of bins used to calculate the bottom degree distribution
        topColor: str, color for plotting the top degree distribution
        bottomColor: str, color for plotting the bottom degree distribution
        
    Returns:
    --------
    topDegrees: dict, degrees of the topNodes (companies)
    topDistribution: np.array, degree distribution (pdf) of top nodes
    topBinCenters: np.array, centers of bins for top nodes
    bottomDegree: dict, degrees of the bottomNodes
    bottomPdf: np.array, degree distribution (pdf) of bottom nodes
    bottomBinCenters: np.array, centers of bins for bottom nodes
    
    Further, saves the degree distributions in a file
    """
    nTopBins = cfg['nTopDegreeBins']
    nBottomBins = cfg['nBottomDegreeBins']
    topColor = cfg['topColor']
    bottomColor = cfg['bottomColor']  
    savePath = cfg['savePathBase'] + cfg['degreeSaveName']
    
    top, bottom = getTopAndBottom(bnet)
    topDegrees = nx.degree(bnet,top)
    topDegrees = dict(topDegrees) # bipartite.degrees returns a DegreeView so this is required for accessing values
    bottomDegrees = nx.degree(bnet, bottom)
    bottomDegrees = dict(bottomDegrees)
    topDegree = topDegrees.values() 
    bottomDegree = bottomDegrees.values()
    topPdf, topBinCenters = getDistribution(topDegree, nTopBins)
    bottomPdf, bottomBinCenters = getDistribution(bottomDegree, nBottomBins)
    
    fig = plt.figure()
    ax = fig.add_subplot(121)
    ax.plot(topBinCenters, topPdf, color=topColor, label='Companies (top nodes)')
    ax.set_xlabel('Degree')
    ax.set_ylabel('PDF')
    ax.set_title('Companies (top nodes)')
    #ax.legend()
    ax = fig.add_subplot(122)
    ax.plot(bottomBinCenters, bottomPdf, color=bottomColor, label='Events (bottom nodes)')
    ax.set_xlabel('Degree')
    ax.set_ylabel('PDF')
    ax.set_title('Events (bottom nodes)')
    #ax.legend()
    
    plt.tight_layout()
    plt.savefig(savePath,format='pdf',bbox_inches='tight')
    
    return topDegrees, topPdf, topBinCenters, bottomDegrees, bottomPdf, bottomBinCenters

def getDensity(bnet):
    """
    Returns density of the bipartite graph.
    
    Parameters:
    -----------
    bnet: networkx.Graph(), bipartite
    
    Returns:
    d: float, density of bnet
    """
    top, _ = getTopAndBottom(bnet)
    d = bipartite.density(bnet,top)
    return d

def findTopNodesInPercentile(bnet,lowPercentile,highPercentile,cfg):
    """
    Finds the top nodes that belong to the given lowest and highest percentiles
    of the degree distribution (low and high percentiles don't need to be equal).
    
    Parameters:
    -----------
    bnet: networkx.Graph(), bipartite
    lowPercentile: int, the percentile for the bottom (left tail) of the distribution
                   (must be between 0 and 100)
    highPercentile: int, the percentile for the top (right tail) of the distribution
                   (must be between 0 and 100)
    cfg: dict, contains:
         savePathBase: str, a base path (e.g. to a shared folder) for saving figures
         degreeSaveName: str, name of the file where to save the degree distribution plots
         nDegreeBins: int, number of bins used to calculate the distributions
         topColor: str, color for plotting the top degree distribution
                   
    Returns:
    --------
    lowPercentileNodes: list, tags of nodes in the low percentile
    highPercentileNodes: list, tags of nodes in the high percentile
    
    Further, saves the top degree distribution with the percentiles marked
    with dashed lines.
    """
    topDegrees,topPdf,topBinCenters,_,_,_ = getDegreeDistributions(bnet,cfg)
    lowPercentileValue = np.percentile(topDegrees.values(),lowPercentile)
    highPercentileValue = np.percentile(topDegrees.values(),highPercentile)
    lowPercentileNodes = []
    highPercentileNodes = []
    for node, degree in topDegrees.items():
        if degree <= lowPercentileValue:
            lowPercentileNodes.append(node)
        elif degree > highPercentileValue:
            highPercentileNodes.append(node)
            
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.plot(topBinCenters, topPdf, color=cfg['topColor'], label='Companies (top nodes)')
    plt.plot([lowPercentileValue,lowPercentileValue],[min(topPdf),max(topPdf)],ls='--',color='k')
    ax.plot([highPercentileValue,highPercentileValue],[min(topPdf),max(topPdf)],ls='--',color='k')
    
    ax.set_xlabel('Degree')
    ax.set_ylabel('PDF')
    ax.set_title('Companies (top nodes)')
    
    savePath = cfg['savePathBase'] + cfg['degreeSaveName'][0:-4] + '_lowPercentile_' + str(lowPercentile) + '_highPercentile_' + str(highPercentile) + cfg['degreeSaveName'][-4::]
    
    plt.tight_layout()
    plt.savefig(savePath,format='pdf',bbox_inches='tight')
    
    return lowPercentileNodes,highPercentileNodes

def getDegreeNodeDictionary(bnet,cfg,nameKey='Member:'):
    """
    Creates a dictionary where keys are degree values of top nodes (companies)
    and values are lists of companies with the given degree. This dictionary
    is also saved as a .csv file.
    
    Parameters:
    -----------
    bnet: nx.Graph(), a bipartite
    nameKey: str, the key corresponding the name of the company in the node
             attributes dictionary (default: 'Member:')
    cfg: dict, contains (at least):
         savePathBase: str, a base path (e.g. to a shared folder) for saving figures
         degreeNodeDictionarySaveName: str, name of the file to which save the dictionary
    
    
    Returns:
    --------
    degreeDict: dict where keys are degrees of top nodes and values are lists of
                companies with the key degree
    """
    top, _ = getTopAndBottom(bnet)
    topDegrees = nx.degree(bnet,top)
    topDegrees = dict(topDegrees) # bipartite.degrees returns a DegreeView so this is required for accessing values
    degreeKeys = np.unique(topDegrees.values())
    degreeDict = {degreeKey:[] for degreeKey in degreeKeys}
    for node, degree in topDegrees.items():
        degreeDict[degree].append(bnet.nodes()[node]['Member:'])
    saveDict = {'Degree:':[],'Number of companies:':[],'Companies:':[]}
    for degree, companies in degreeDict.items():
        saveDict['Degree:'].append(degree)
        saveDict['Number of companies:'].append(len(companies))
        saveDict['Companies:'].append(companies)
    df = pd.DataFrame(saveDict,columns=['Degree:','Number of companies:','Companies:'])
    savePath = cfg['savePathBase'] + cfg['degreeNodeDictionarySaveName']
    df.to_csv(savePath,index=None,header=True)
    
    return degreeDict
                

# Clique analysis
    
def stuffNetwork(bnet):
    """
    Adds links between nodes of same type so that top and bottom nodes from two
    complete subgraphs. Note that the output is no more a bipartite graph.
    
    Parameters:
    -----------
    bnet: networkx.Graph(), bipartite
    
    Returns:
    --------
    stuffedNet: networkx.Graph(), a general graph with same nodes as bnet
    """
    top, bottom = getTopAndBottom(bnet)
    topLinks = list(combinations(top,2))
    bottomLinks = list(combinations(bottom,2))
    stuffedNet = nx.Graph()
    stuffedNet.add_nodes_from(bnet.nodes(data=True))
    stuffedNet.add_edges_from(bnet.edges)
    stuffedNet.add_edges_from(topLinks)
    stuffedNet.add_edges_from(bottomLinks)
    
    return stuffedNet

def findBicliques(bnet):
    """
    Finds the maximal bicliques of the bipartite network bnet. Biclique is a complete subgraph 
    of bnet so that each member of a biclique is connected to each other member 
    of the biclique that are from a different node class. Maximal bicliques are
    not contained by any larger biclique. To find the bicliques, the network is
    first stuffed to a general network by adding links between all nodes of same class,
    and then the cliques of this general network are detected. The method is based on Makino &
    Uno 2004: New algorithms for enumerating all maximal cliques. Note that this
    algorithm does not scale nicely and may take long time to run for large networks.
    
    Parameters:
    -----------
    bnet: networkx.Graph(), a bipartite
    
    Returns:
    --------
    cliques: list of lists, cliques of bnet; each clique is represented as a 
             list containing the nodes of the clique
    cliqueInfo: list of dicts; each dict contains one clique separated to top and
                bottom nodes (keys: 'topNodes', 'bottomNodes')
    
    """
    snet = stuffNetwork(bnet)
    cliques = list(nx.find_cliques(snet)) # this finds all cliques, including those formed by top and bottom nodes
    top, bottom = getTopAndBottom(bnet) 
    for clique in cliques: # remocing the top and bottom cliques
        if set(clique) == top:
            cliques.remove(clique)
            print 'removed top:'
            print clique
            break
    for clique in cliques:
        if set(clique) == bottom:
            cliques.remove(clique)
            print 'removed bottom:'
            print clique
            break
    
    cliqueInfo = getCliqueIndices(bnet,cliques)
        
    return cliques, cliqueInfo

def getCliqueIndices(bnet, cliques):
    """
    Calculates the number of top and bottom nodes in each of the bicliques of a
    bipartite network.
    
    Parameters:
    -----------
    bnet: networkx.Graph(), a bipartite
    cliques: list of lists, cliques of bnet; each clique is represented as a 
             list containing the nodes of the clique
             
    Returns:
    --------
    cliqueInfo: list of dicts; each dict contains one clique separated to top and
                bottom nodes (keys: 'topNodes': top nodes (companies) of the clique, 
                                    'bottomNodes' : bottom nodes (events) of the clique,
                                    'topIndex': number of top nodes of the clique,
                                    'bottomIndex': number of bottom nodes of the clique,
                                    'isBridge': does the clique consist of one top node (company) connecting multiple bottom nodes,
                                    'isStar': does the clique consist of node bottom node (event) connecting multiple top nodes)
    """
    top, bottom = getTopAndBottom(bnet)
    
    cliqueInfo = []
    
    for i, clique in enumerate(cliques):
        info = {}
        cliqueTop = set(clique) & top
        cliqueBottom = set(clique) & bottom
        info['topNodes'] = cliqueTop
        info['bottomNodes'] = cliqueBottom
        info['topIndex'] = len(cliqueTop)
        info['bottomIndex'] = len(cliqueBottom)
        if len(cliqueTop) == 1:
            info['isBridge'] = True
        else:
            info['isBridge'] = False
        if len(cliqueBottom) == 1:
            info['isStar'] = True
        else:
            info['isStar'] = False
        cliqueInfo.append(info)
    
    return cliqueInfo
    
def pruneStars(bnet, cliques, cliqueInfo, nonMemberClass='NM'):
    """
    Removes from bistars all top nodes (companies) that participate in any other
    cliques besides the bistar, as well as all top nodes whose membership class
    matches the non-member class. Bistar is a clique where one bottom node (event)
    is surrounded by multiple top nodes (companies). If all top nodes of the bistar
    are removed, the whole star is removed.
    
    Parameters:
    -----------
    bnet: networx.Graph()
        a bipartite
    cliques: list of lists
        cliques of bnet; each clique is represented as a 
        list containing the nodes of the clique
    cliqueInfo: list of dicts; each dict contains one clique separated to top and
                bottom nodes (keys: 'topNodes': top nodes (companies) of the clique, 
                                    'bottomNodes' : bottom nodes (events) of the clique,
                                    'topIndex': number of top nodes of the clique,
                                    'bottomIndex': number of bottom nodes of the clique,
                                    'isBridge': does the clique consist of one top node (company) connecting multiple bottom nodes,
                                    'isStar': does the clique consist of node bottom node (event) connecting multiple top nodes)
    nonMemberClass: str
        the membership class (class attribute) of those nodes that should be
        removed from stars. Non-member companies or instances often
        participate only one event and therefore artificially increase starness.
    
    Returns:
    --------
    cliques: list of lists, the original set of cliques with bistars are replaced
             with the pruned versions
    cliqueInfo: list of dicts, updated info for the updated clique list
    """
    cliquesToRemove = []
    infosToRemove = []
    bottomOnlyCliques = []
    bottomOnlyInfos = []
    newCliques = []
    newInfos = []
    for clique, info in zip(cliques,cliqueInfo):
        if info['isStar'] == True:
            topToRemove = []
            for top in info['topNodes']:
                if nx.degree(bnet,top) > 1:
                    topToRemove.append(top)
                if bnet.nodes(data=True)[top]['class'] == nonMemberClass:
                    topToRemove.append(top)
            if len(topToRemove) == len(info['topNodes']): # if all the top nodes of a clique are to be removed, let's remove the whole clique
                bottomOnlyCliques.append(clique)
                bottomOnlyInfos.append(info)
            else:
                cliquesToRemove.append(clique)
                infosToRemove.append(info)
                newClique = list(clique)
                newInfo = dict(info)
                newInfo['topIndex'] = newInfo['topIndex'] - len(topToRemove)
                newInfo['topNodes'] = set(list(newInfo['topNodes']))
                for top in topToRemove:
                    newClique.remove(top)
                    newInfo['topNodes'].remove(top)
                newCliques.append(newClique)
                newInfos.append(newInfo)
    for clique, info, newClique, newInfo in zip(cliquesToRemove, infosToRemove, newCliques, newInfos):
        cliques.remove(clique)
        cliqueInfo.remove(info)
        cliques.append(newClique)
        cliqueInfo.append(newInfo)
    for bottomOnlyClique, bottomOnlyInfo in zip(bottomOnlyCliques,bottomOnlyInfos):
        cliques.remove(bottomOnlyClique)
        cliqueInfo.remove(bottomOnlyInfo)
    
    return cliques, cliqueInfo

def createCliqueIndexHeatmap(cliqueInfo, cfg):
    """
    Creates a 2D histogram of the number of bicliques with different numbers of
    top and bottom nodes and visualizes it as a heatmap.
    
    Parameters:
    -----------
    cliqueInfo: list of dicts; each dict contains one clique separated to top and
                bottom nodes (keys: 'topNodes': top nodes (companies) of the clique, 
                                    'bottomNodes' : bottom nodes (events) of the clique,
                                    'topIndex': number of top nodes of the clique,
                                    'bottomIndex': number of bottom nodes of the clique,
                                    'isBridge': does the clique consist of one top node (company) connecting multiple bottom nodes,
                                    'isStar': does the clique consist of node bottom node (event) connecting multiple top nodes)
    cfg: dict, containing:
        cliqueHeatmapCmap: str, colormap to be used for the heatmap
        cliqueHeatmapTopBins: list of ints, bin edges for the top (company) index
        cliqueHeatmapBottomBins: list of ints, bin edges for the bottom (event) index
        savePathBase: str, a base path (e.g. to a shared folder) for saving figures
        cliqueHeatmapSaveName: str, name of the file where to save the heatmap
        
    Returns:
    --------
    no direct output, saves the heatmap to the given path
    """    
    topIndices = []
    bottomIndices = []
    for clique in cliqueInfo:
        topIndices.append(clique['topIndex'])
        bottomIndices.append(clique['bottomIndex'])
        
    topBins = cfg['cliqueHeatmapTopBins']
    bottomBins = cfg['cliqueHeatmapBottomBins']
    count,xEdges,yEdges,_ = binned_statistic_2d(topIndices, bottomIndices, bottomIndices,statistic='count',bins=[topBins,bottomBins])
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    im = ax.imshow(count,interpolation='none',cmap=cfg['cliqueHeatmapCmap'],aspect='auto',origin='lower')
    ax.autoscale(False)
    ax.set_xlabel('Bottom index (number of events)')
    ax.set_ylabel('Top index (number of companies)')
    ax.set_xticks(cfg['cliqueHeatmapBottomTicks'])
    ax.set_yticks(cfg['cliqueHeatmapTopTicks'])
    ax.set_xticklabels(cfg['cliqueHeatmapBottomLabels'])
    ax.set_yticklabels(cfg['cliqueHeatmapTopLabels'])
    ax.tick_params(top='off',right='off')
    cbar = ax.figure.colorbar(im,ax=ax)
    cbar.ax.set_ylabel('Count', rotation=-90, va="bottom")
    
    #ax.add_colorbar()
    fig.tight_layout()
    
    savePath = cfg['savePathBase'] + cfg['cliqueHeatmapSaveName']
    plt.savefig(savePath,format='pdf',bbox_inches='tight')
        
def getStarness(bnet,cliqueInfo):
    """
    Calculates the starness of a bipartite graph. Starness is defined as the 
    fraction of top nodes (companies) belonging to bistars out of all top nodes
    of the graph. Bistar is a clique where one bottom node is surrounded by multiple top nodes.
    
    Parameters:
    -----------
    bnet: networkx.Graph(), a bipartite
    cliqueInfo: list of dicts; each dict contains one clique separated to top and
                bottom nodes (keys: 'topNodes': top nodes (companies) of the clique, 
                                    'bottomNodes' : bottom nodes (events) of the clique,
                                    'topIndex': number of top nodes of the clique,
                                    'bottomIndex': number of bottom nodes of the clique,
                                    'isBridge': does the clique consist of one top node (company) connecting multiple bottom nodes,
                                    'isStar': does the clique consist of node bottom node (event) connecting multiple top nodes)
    
    Returns:
    --------
    starness: double, starness of the bigraph
    """
    nStarNodes = 0
    top, _ = getTopAndBottom(bnet) 
    nTotal = len(top)
    for clique in cliqueInfo:
        if clique['isStar']:
            nStarNodes = nStarNodes + clique['topIndex']
    starness = nStarNodes/float(nTotal)
    return starness
    
def getCliqueFieldDiversity(bnet,cliqueInfo):
    """
    Calculates the diversity of a biclique in terms of the fields of business
    of the participating top nodes (companies). For obtaining the diversity, the
    Gini-Simpson index transformed to effective diversity is used:
    GS = sum_1^S(p_1^2), D = 1/(1-GS)
    
    Parameters:
    -----------
    bnet: networkx.Graph(), a bipartite
    cliqueInfo: dict, contains:
        topNodes: list of nodes, top nodes (companies) of the clique, 
        bottomNodes : list of nodes, bottom nodes (events) of the clique,
        topIndex: int, number of top nodes of the clique,
        bottomIndex: int, number of bottom nodes of the clique,
        isBridge: boolean, does the clique consist of one top node (company) connecting multiple bottom nodes,
        isStar: boolean, does the clique consist of node bottom node (event) connecting multiple top nodes
    
    Returns:
    --------
    richness: int, number of different fields in the clique
    cliqueDiversity: double, the effective diversity based on Gini-Simpson index
    count: dict, number of abundance of different fields in the clique
    majorField: str, the most common field tag in the clique
    """
    fieldsOfBusiness = []
    topNodes = cliqueInfo['topNodes']
    nTops = len(topNodes)
    allFields = nx.get_node_attributes(bnet,'tag')
    for top in topNodes:
        fieldsOfBusiness.append(allFields[top])
    uniqueFields = set(fieldsOfBusiness)
    richness = len(uniqueFields)
    count = {}
    for field in fieldsOfBusiness:
        count[field] = count.get(field,0) + 1
    fractions = []
    for field, abundance in count.iteritems():
        fractions.append(abundance/float(nTops))
    GS = 1 - sum(np.array(fractions)**2)
    cliqueDiversity = 1/(1-GS)
    majorField = max(count,key=count.get)
    nMajor = count[majorField]
    majorFraction = nMajor/float(len(topNodes))
    print 'For the present clique, effective diversity: ' + str(cliqueDiversity) + ', major field: ' + majorField + ', ' + str(nMajor) + ', ' + str(majorFraction) + ' of all (' + str(len(topNodes)) + ') companies'
    return richness, cliqueDiversity, count, majorField
    
def getCliqueFieldDiversityWrapper(bnet,cliqueInfo):
    """
    As a wrapper, handles the loop over cliques to obtain the biclique diversity
    using getCliqueFieldDiversity.
    
    Parameters:
    -----------
    bnet: networkx.Graph(), a bipartite
    cliqueInfo: list of dicts, each of them containing:
        topNodes: list of nodes, top nodes (companies) of the clique, 
        bottomNodes : list of nodes, bottom nodes (events) of the clique,
        topIndex: int, number of top nodes of the clique,
        bottomIndex: int, number of bottom nodes of the clique,
        isBridge: boolean, does the clique consist of one top node (company) connecting multiple bottom nodes,
        isStar: boolean, does the clique consist of node bottom node (event) connecting multiple top nodes)
    
    Returns:
    --------
    richenesses: list of ints, numbers of different fields in the cliques
    cliqueDiversities: list of doubles, the effective diversities based on Gini-Simpson index
    counts: list of dicts, numbers of abundance of different fields in the cliques
    majorFields: list of strs, the most common field tags in the cliques
    """
    richnesses = []
    cliqueDiversities = []
    counts = []
    majorFields = []
    for clique in cliqueInfo:
        richness, diversity, count, majorField = getCliqueFieldDiversity(bnet,clique)
        richnesses.append(richness)
        cliqueDiversities.append(diversity)
        counts.append(count)
        majorFields.append(majorField)
    return richnesses, cliqueDiversities, counts, majorFields
    
    
# Null models
    
def createRandomBipartite(bnet):
    """
    Creates a randomly wired bipartite network with the same number of top and
    bottom nodes and the same density as a given network.
    
    Parameters:
    -----------
    bnet: networkx.Graph(), a bipartite
    
    Returns:
    --------
    randNet: networkx.Graph(), a randomly wired bipartite
    """
    top, bottom = getTopAndBottom(bnet) 
    nTop = len(top)
    nBottom = len(bottom)
    nEdges = len(bnet.edges())
    randNet = nx.bipartite.gnmk_random_graph(nTop,nBottom,nEdges)
    return randNet

def shuffleFields(bnet):
    """
    Creates a version of bnet with the tags (fields) of top nodes (companies) randomly shuffled.
    The link configuration of randNet is the same as in bnet.
    
    Parameters:
    -----------
    bnet: networkx.Graph(), a bipartite
    
    Returns:
    --------
    randNet: networkx.Graph(), a shuffled version of bnet
    """
    top, bottom = getTopAndBottom(bnet)
    tags = nx.get_node_attributes(bnet,'tag')
    colors = nx.get_node_attributes(bnet,'nodeColor')
    topAttributes = []
    for topNode in top:
        topAttributes.append((tags[topNode],colors[topNode]))
    shuffledAttributes = list(topAttributes) # shuffling happens in-place so a copy of the list is needed
    random.shuffle(shuffledAttributes)
    newTags = {}
    newColors = {}
    for topNode, attribute in zip(top,shuffledAttributes):
        newTags[topNode] = attribute[0]
        newColors[topNode] = attribute[1]
    for bottomNode in bottom:
        newTags[bottomNode] = tags[bottomNode]
        newColors[bottomNode] = colors[bottomNode]
    randNet = nx.Graph()
    randNet.add_nodes_from(bnet.nodes(data=True))
    nx.set_node_attributes(randNet,newTags,'tag')
    nx.set_node_attributes(randNet,newColors,'color')
    randNet.add_edges_from(bnet.edges())
    
    return randNet
    
def compareAgainstRandom(bnet,cfg,measures):
    """
    Compares the measures obtained from bnet agains a suitable null model (edge-shuffled
    for starness, tag-shuffled for richness and diversity).
    
    Parameters:
    -----------
    bnet: networkx.Graph(), the original bipartite network
    cfg: dict, containing:
               nRandomIterations: int, number of null model instances to be used
               nRandomBins: int, number of bins for obtaining the random distribution
               nRichnessBins: int, number of bins for obtaining distributions of richness
               randomColor: str, color for visualizing the values obtained from random networks
               dataColor: str, color for visualizing the actual values
               randomMarker: str, marker for the values obtained from random networks
               dataMarker: str, marker for the actual values
               dataLineWidth: double, width of the lines presenting data (increased from default to increase data-random contrast)
               randomAlpha: double, transparency value for the data points from random networks
               identityLineStyle: str, line style for plotting the identity line
               savnRandomBins: int, number of bins for obtaining the random distributionePathBase: str, a base path (e.g. to a shared folder) for saving figures
               comparisonVsRandomSaveName: str, name of the file where to save visualizations of the comparison
    measures, dict, possible keys:
                    starness: double, starness of the bigraph
                    cliques: list of lists, cliques of the bipartite (each clique presented as a list of nodes)
                    richness: list of ints, numbers of different fields in the cliques
                    diversity: list of doubles, the effective diversities based on Gini-Simpson index
                    
    Returns:
    --------
    No direct output, saves the test output as figures to the given paths.
    
    """
    nIters = cfg['nRandomIterations']
    randColor = cfg['randomColor']
    dataColor = cfg['dataColor']
    randMarker = cfg['randomMarker']
    randAlpha = cfg['randomAlpha']
    dataMarker = cfg['dataMarker']
    dataLineWidth = cfg['dataLineWidth']
    nRandBins = cfg['nRandomBins']
    savePathBase = cfg['savePathBase']
    saveNameBase = cfg['comparisonVsRandomSaveName']  
    
    if 'starness' in measures.keys():
        starness = measures['starness']
        randStarness = []
        for i in range(nIters):
            randNet = createRandomBipartite(bnet)
            randNet,_ = pruneBipartite(randNet)
            cliques, cliqueInfo = findBicliques(randNet)
            cliques, cliqueInfo = pruneStars(randNet,cliques,cliqueInfo)
            randStarness.append(getStarness(bnet,cliqueInfo))
        randLabel = 'Starness of random networks, mean: ' + str(np.mean(randStarness))
        dataLabel = 'True starness: ' + str(starness)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        pdf, binCenters = getDistribution(randStarness,nRandBins)
        plt.plot(binCenters,pdf,color=randColor,alpha=randAlpha,label=randLabel)
        plt.plot([starness,starness],[0,max(pdf)],color=dataColor,label=dataLabel)
        ax.set_xlabel('Starness')
        ax.set_ylabel('PDF')
        ax.legend()
        plt.tight_layout()
        savePath = savePathBase + saveNameBase + '_starness.pdf'
        plt.savefig(savePath,format='pdf',bbox_inches='tight')
    if 'richness' in measures.keys():
        richness = measures['richness']
        nBins = cfg['nRichnessBins']
        trueDist, trueBinCenters = getDistribution(richness,nBins)
        sizes = [len(clique) for clique in measures['cliques']]
        randDists = []
        randBinCenters = []
        randRichnesses = []
        meanRichnesses = []
        if 'diversity' in measures.keys():
            diversity = measures['diversity']
            trueDiverDist, trueDiverBinCenters = getDistribution(diversity,nBins)
            randDiverDists = []
            randDiverBinCenters = []
            randDiversities = []
            randMeanDiversities = []
            randRelativeDiversities = []
        for i in range(nIters):
            randNet = shuffleFields(bnet)
            cliques,cliqueInfo = findBicliques(randNet)
            randRichness, randDiversity,_,_ = getCliqueFieldDiversityWrapper(randNet,cliqueInfo)
            randDist,binCenters = getDistribution(randRichness,nBins)
            randDists.append(randDist)
            randBinCenters.append(binCenters)
            randRichnesses.extend(randRichness)
            meanRichnesses.append(np.mean(np.array(randRichness)))
            if 'diversity' in measures.keys():
                randDiverDist,diverBinCenters = getDistribution(randDiversity,nBins)
                randDiverDists.append(randDiverDist)
                randDiverBinCenters.append(diverBinCenters)
                randDiversities.extend(randDiversity)
                randMeanDiversities.append(np.mean(np.array(randDiversity)))
                randRelativeDiversities.append(np.mean(np.array(randDiversity)/np.array(randRichness)))
                
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for dist, centers in zip(randDists,randBinCenters):
            plt.plot(centers,dist,color=randColor,alpha=randAlpha)
        plt.plot(trueBinCenters,trueDist,color=dataColor,linewidth=dataLineWidth,label='True richness')
        ax.set_xlabel('Richness')
        ax.set_ylabel('PDF')
        ax.set_title('Richness')
        ax.legend()
        plt.tight_layout()
        savePath = savePathBase + saveNameBase + '_richness_dist.pdf'
        plt.savefig(savePath,format='pdf',bbox_inches='tight')

        fig = plt.figure()
        ax = fig.add_subplot(111)
        randMeanDist, randMeanCenters = getDistribution(meanRichnesses,nRandBins)
        randLabel = 'mean richness in random networks, mean of means: ' + str(np.mean(meanRichnesses))
        dataLabel = 'mean richess in data: ' + str(np.mean(richness))
        plt.plot(randMeanCenters,randMeanDist,color=randColor,alpha=randAlpha,label=randLabel)
        plt.plot([np.mean(richness),np.mean(richness)],[0,max(randDist)],color=dataColor,label=dataLabel)
        ax.set_xlabel('Mean richness')
        ax.set_ylabel('PDF')
        ax.legend()
        plt.tight_layout()
        savePath = savePathBase + saveNameBase + '_mean_richness_dist.pdf'
        plt.savefig(savePath,format='pdf',bbox_inches='tight')
        
        if 'diversity' in measures.keys():
            fig = plt.figure()
            ax = fig.add_subplot(111)
            for dist, centers in zip(randDiverDists,randDiverBinCenters):
                plt.plot(centers,dist,color=randColor,alpha=randAlpha)
            plt.plot(trueDiverBinCenters,trueDiverDist,color=dataColor,linewidth=dataLineWidth,label='True diversity')
            ax.set_xlabel('Effective diversity')
            ax.set_ylabel('PDF')
            ax.set_title('Effective diversity')
            ax.legend()
            plt.tight_layout()
            savePath = savePathBase + saveNameBase + '_effective_diversity_dist.pdf'
            plt.savefig(savePath,format='pdf',bbox_inches='tight')
            
            fig = plt.figure()
            ax = fig.add_subplot(111)
            randDist,randCenters = getDistribution(randMeanDiversities,nRandBins)
            trueMean = np.mean(diversity)
            randomMean = np.mean(randMeanDiversities)
            randLabel = 'mean effective diversity in random networks, pooled: ' + str(randomMean)
            dataLabel = 'mean effective diversity in data: ' + str(trueMean)
            plt.plot(randCenters,randDist,color=randColor,alpha=randAlpha,label=randLabel)
            plt.plot([trueMean,trueMean],[0,max(randDist)],color=dataColor,label=dataLabel)
            ax.set_xlabel('Mean effective diversity')
            ax.set_ylabel('PDF')
            ax.legend()
            plt.tight_layout()
            savePath = savePathBase + saveNameBase + '_mean_effective_diversity_dist.pdf'
            plt.savefig(savePath,format='pdf',bbox_inches='tight')
            
            fig = plt.figure()
            ax = fig.add_subplot(111)
            randDist,randCenters = getDistribution(randRelativeDiversities,nRandBins)
            relativeDiversity = np.array(diversity)/np.array(richness)
            trueMean = np.mean(relativeDiversity)
            randomMean = np.mean(randRelativeDiversities)
            randLabel = 'mean relative diversity in random networks, pooled: ' + str(randomMean)
            dataLabel = 'mean relative diversity in data: ' + str(trueMean)
            plt.plot(randCenters,randDist,color=randColor,alpha=randAlpha,label=randLabel)
            plt.plot([trueMean,trueMean],[0,max(randDist)],color=dataColor,label=dataLabel)
            ax.set_xlabel('Mean relative diversity')
            ax.set_ylabel('PDF')
            ax.legend()
            plt.tight_layout()
            savePath = savePathBase + saveNameBase + '_relative_diversity_dist.pdf'
            plt.savefig(savePath,format='pdf',bbox_inches='tight')
            
            fig = plt.figure()
            ax = fig.add_subplot(111)
            identity = range(min(randRichnesses),max(randRichnesses))
            plt.plot(randRichnesses,randDiversities,color=randColor,marker=randMarker,alpha=randAlpha,ls='',label='random')
            plt.plot(identity,identity,color='k',ls=cfg['identityLineStyle'],label='x = y')
            plt.plot(richness,diversity,color=dataColor,marker=dataMarker,ls='',label='Data')
            ax.set_xlabel('Richness')
            ax.set_ylabel('Effective diversity')
            ax.legend()
            plt.tight_layout()
            savePath = savePathBase + saveNameBase + '_diversity_vs_richness.pdf'
            plt.savefig(savePath,format='pdf',bbox_inches='tight')
            
            fig = plt.figure()
            ax = fig.add_subplot(111)
            randRelativeDiversity = np.array(randDiversities)/np.array(randRichnesses)
            plt.plot(np.tile(sizes,nIters),randRelativeDiversity,color=randColor,marker=randMarker,alpha=randAlpha,ls='',label='random')
            plt.plot(sizes,relativeDiversity,color=dataColor,marker=dataMarker,ls='',label='data')
            ax.set_xlabel('Clique size')
            ax.set_ylabel('Relative diversity')
            ax.legend()
            plt.tight_layout()
            savePath = savePathBase + saveNameBase + '_relative_diversity_vs_size.pdf'
            plt.savefig(savePath,format='pdf',bbox_inches='tight')
        
    elif 'diversity' in measures.keys():
        diversity = measures['diversity']
        trueDiverDist, trueDiverBinCenters = getDistribution(diversity,nBins)
        randDiverDists = []
        randDiverBinCenters = []
        randDiversities = []
        for i in range(nIters):
            randNet = shuffleFields(bnet)
            cliques,cliqueInfo = findBicliques(randNet)
            _,randDiversity,_,_ = getCliqueFieldDiversityWrapper(randNet,cliqueInfo)
            randDiverDist,diverBinCentres = getDistribution(randDiversity,nBins)
            randDiverDists.append(randDiverDist)
            randDiverBinCenters.append(diverBinCenters)
            randDiversities.append(randDiversity)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for dist, centers in zip(randDiverDists,randDiverBinCenters):
                plt.plot(centers,dist,color=randColor,alpha=randAlpha)
        plt.plot(trueDiverBinCenters,trueDiverDist,color=dataColor,linewidth=dataLineWidth,label='True diversity')
        ax.set_xlabel('Effective diversity')
        ax.set_ylabel('PDF')
        ax.set_title('Effective diversity')
        ax.legend()
        savePath = savePathBase + saveNameBase + '_diversity_dist.pdf'
        plt.savefig(savePath,format='pdf',bbox_inches='tight')

        
            
        
            
    
    
# Visualization
    
def drawNetwork(bnet, cfg):
    """
    Visualizes the bipartite network and saves the figure as pdf.
    
    Parameters:
    -----------
    bnet: networkx.Graph(), the bipartite to be visualized
    cfg: dict, contains:
        tags: list of strs
            the field of business tags
        classes: list of strs
            the possible membership classes
        networkColors: matplotlib.cmap
            colors associated with different fields of business
        bottomNetworkColor: str
            color of the bottom nodes (events)
        nodeShapes: list of strs
            shapes associated with the different membership classes
        bottomShape: str
            the shape of bottom (event) nodes
        nodeSize: int
            size of nodes in the visualization
        edgeWidth: double
            width of network edges
        savePathBase: str
            a base path (e.g. to a shared folder) for saving figures
        networkSaveName: str
            name of the file where to save the network visualization
        
    Returns:
    --------
    no direct output, saves the network visualization as pdf
    """
    
    tags = cfg['tags']    
    classes = cfg['classes']
    networkColors = cfg['networkColors']
    bottomColor = cfg['networkBottomColor']
    nodeShapes = cfg['nodeShapes']
    bottomShape = cfg['bottomShape']
    nodeSize = cfg['nodeSize']
    edgeWidth = cfg['edgeWidth']
    
    top, bottom = getTopAndBottom(bnet)
    
    pos = nx.spring_layout(bnet)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    nx.draw_networkx_nodes(bnet,pos,ax=ax,nodelist=bottom,node_color=bottomColor,node_shape=bottomShape,
                           node_size=nodeSize,label='Events)')
    
    for i, tag in enumerate(tags):
        for j, mclass in enumerate(classes):
            taggedNodes = []
            for topn in top:
                if bnet.nodes(data=True)[topn]['tag'] == tag and bnet.nodes(data=True)[topn]['class'] == mclass:
                    taggedNodes.append(topn)
        nx.draw_networkx_nodes(bnet,pos,ax=ax,nodelist=taggedNodes,node_color=networkColors(i),node_shape=nodeShapes[j],
                               node_size=nodeSize,label=tag)
    
    nx.draw_networkx_edges(bnet,pos,width=edgeWidth)
    #ax.legend()
    plt.axis('off')    
    
    savePath = cfg['savePathBase'] + cfg['networkSaveName']
    plt.savefig(savePath,format='pdf',bbox_inches='tight')
    
def visualizeBicliques(bnet, cliqueInfo, cfg):
    """
    For each (bi)clique of a network, visualized the network so that members of
    the clique are colored in one color and rest of the network in another.
    
    Parameters:
    -----------
    bnet: networkx.Graph(), a bipartite
    cliques: list of lists, cliques of bnet; each clique is represented as a 
             list containing the nodes of the clique
    cfg: dict, contains:
        bottomNetworkColor: str, color of the bottom nodes (events) in the clique
        cliqueTopColor: str, color of the top nodes (companies) in the clique
        nonCliqueColor: str, color of the nodes not belonging to the clique
        nonCliqueAlpha: double, transparency of the non-clique nodes
        nodeSize: int, size of nodes in the visualization
        edgeWidth: double, width of network edges
        savePathBase: str, a base path (e.g. to a shared folder) for saving figures
        cliquesSaveName: str, name of the file where to save the clique visualization
        
    Returns:
    --------
    no direct output, saves the visualizations to given path
    """
    bottomColor = cfg['networkBottomColor']
    topColor = cfg['cliqueTopColor']
    nonCliqueColor = cfg['nonCliqueColor']
    nonCliqueAlpha = cfg['nonCliqueAlpha']
    nodeSize = cfg['nodeSize']
    edgeWidth = cfg['edgeWidth']
    
    saveName = cfg['cliquesSaveName']
    
    pos = nx.spring_layout(bnet)
    
    for i, clique in enumerate(cliqueInfo):
        topCliqueNodes = clique['topNodes']
        bottomCliqueNodes = clique['bottomNodes']
        nonCliqueNodes = set(bnet) - topCliqueNodes.intersection(bottomCliqueNodes)
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        nx.draw_networkx_nodes(bnet,pos=pos,nodelist=topCliqueNodes,node_size=nodeSize,node_color=topColor)
        nx.draw_networkx_nodes(bnet,pos=pos,nodelist=bottomCliqueNodes,node_size=nodeSize,node_color=bottomColor)
        nx.draw_networkx_nodes(bnet,pos=pos,nodelist=nonCliqueNodes,node_size=nodeSize,node_color=nonCliqueColor,alpha=nonCliqueAlpha)
        nx.draw_networkx_edges(bnet,pos=pos,width=edgeWidth)
        plt.axis('off')  
        
        savePath = cfg['savePathBase'] + saveName + '_' + str(i) + '.pdf'
        plt.savefig(savePath,format='pdf',bbox_inches='tight')
        
def plotRichnessVsDiversity(richnesses,diversities,cfg):
    """
    Plots the diversity of cliques as a function of the number of different fields
    in the cliques.
    
    Parameters:
    -----------
    richnesses: list of ints, numbers of different fields in the cliques
    diversities: list of doubles, effective diversities of the cliques
    cfg: a dictionary containing:
        savePathBase: str, a base path (e.g. to a shared folder) for saving figures
        diversitySaveName: str, name of the file where to save the visualization
        identityLineStyle: str, line style for plotting the identity line
        scatterMarker: str, point marker style
        
    Returns:
    --------
    no direct output, saves the figure to the given path
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    identity = range(min(richnesses),max(richnesses))
    plt.plot(identity,identity,ls=cfg['identityLineStyle'],label='x = y')
    plt.plot(richnesses,diversities,cfg['scatterMarker'],ls='',label='data')
    ax.set_xlabel('Number of fields')
    ax.set_ylabel('Effective diversity')
    ax.legend()
    fig.tight_layout()
    
    savePath = cfg['savePathBase'] + cfg['diversitySaveName']
    plt.savefig(savePath,format='pdf',bbox_inches='tight')
    
def plotRelativeDiversity(cliques,richnesses,diversities,cfg):
    """
    Plots the relative diversity as a function of clique size. The relative
    diversity is defined as effective diversity/richness and tells how far
    away a clique is from it's theoretical maximal diversity. It's plotted against
    clique size as richness strongly depends on it.
    
    Parameters:
    -----------
    cliques: list of lists, cliques of the network
    richnesses: list of ints, numbers of different fields in the cliques
    diversities: list of doubles, effective diversities of the cliques
    cfg: a dictionary containing:
        savePathBase: str, a base path (e.g. to a shared folder) for saving figures
        relativeDiversitySaveName: str, name of the file where to save the visualization
        identityLineStyle: str, line style for plotting the identity line
        scatterMarker: str, point marker style
    
    Returns:
    --------
    no direct output, saves the figure to the given path
    """
    relativeDiversity = np.array(diversities)/np.array(richnesses)
    sizes = [len(clique) for clique in cliques]
    identity = range(min(sizes),max(sizes))
    fig = plt.figure()
    ax = fig.add_subplot(121)
    plt.plot(identity,identity,ls=cfg['identityLineStyle'],label='x = y')
    plt.plot(sizes,richnesses,ls='',marker=cfg['scatterMarker'],label='richness')
    ax.set_xlabel('Clique size')
    ax.set_ylabel('Richness')
    ax.legend()
    fig.tight_layout()
   
    ax = fig.add_subplot(122)
    
    plt.plot(sizes,relativeDiversity,ls='',marker=cfg['scatterMarker'],label='relative diversity')
    ax.set_xlabel('Clique size')
    ax.set_ylabel('Relative diversity')
    ax.legend()
    fig.tight_layout()
    
    savePath = cfg['savePathBase'] + cfg['relativeDiversitySaveName']
    plt.savefig(savePath,format='pdf',bbox_inches='tight')

def plotDiversityVsIndices(cliqueInfo,richnesses,diversities,cfg):
    """
    Plots the clique diversities (richness, effective diversity, and effective
    diversity) as a function of the number of top and bottom nodes in the clique.
    
    Parameters:
    -----------
    cliqueInfo: list of dicts, each of them containing:
        topNodes: list of nodes, top nodes (companies) of the clique, 
        bottomNodes : list of nodes, bottom nodes (events) of the clique,
        topIndex: int, number of top nodes of the clique,
        bottomIndex: int, number of bottom nodes of the clique,
        isBridge: boolean, does the clique consist of one top node (company) connecting multiple bottom nodes,
        isStar: boolean, does the clique consist of node bottom node (event) connecting multiple top nodes)
    richnesses: list of ints, numbers of different fields in the cliques
    diversities: list of doubles, effective diversities of the cliques
    cfg: dict, contains:
        savePathBase: str, a base path (e.g. to a shared folder) for saving figures
        diversityVsBottomIndexSaveName: str, file name for saving the figure of diversity vs bottom index
        diversityVsTopIndexSaveName: str, file name for saving the figure of diveristy vs top index
        scatterMarker: str, point marker style
        
    Returns:
    --------
    No direct output, saves the plot to the given path
    """
    bottomSavePath = cfg['savePathBase'] + cfg['diversityVsBottomIndexSaveName']
    topSavePath = cfg['savePathBase'] + cfg['diversityVsTopIndexSaveName']
    scatterMarker = cfg['scatterMarker']
    bottomIndices = [info['bottomIndex'] for info in cliqueInfo]
    topIndices = [info['topIndex'] for info in cliqueInfo]
    cliqueSizes = [top + bottom for top,bottom in zip(topIndices,bottomIndices)]
    
    print 'Max number of events ' + str(max(bottomIndices))
    print 'Max number of companies ' + str(max(topIndices))
    
    bottomFig = plt.figure()
    bottomAx1 = bottomFig.add_subplot(221)
    bottomAx1.plot(bottomIndices,richnesses,marker=scatterMarker,ls='')
    bottomAx1.set_xlabel('Number of bottom nodes (events)')
    bottomAx1.set_ylabel('Richness')
    bottomAx2 = bottomFig.add_subplot(222)
    bottomAx2.plot(bottomIndices,diversities,marker=scatterMarker,ls='')
    bottomAx2.set_xlabel('Number of bottom nodes (events)')
    bottomAx2.set_ylabel('Effective diversity')
    bottomAx3 = bottomFig.add_subplot(223)
    relativeDiversity = np.array(diversities)/np.array(richnesses)
    bottomAx3.plot(bottomIndices,relativeDiversity,marker=scatterMarker,ls='')
    bottomAx3.set_xlabel('Number of bottom nodes (events)')
    bottomAx3.set_ylabel('Relative diversity')
    bottomAx4 = bottomFig.add_subplot(224)
    bottomAx4.plot(bottomIndices,cliqueSizes,marker=scatterMarker,ls='')
    bottomAx4.set_xlabel('Number of bottom nodes (events)')
    bottomAx4.set_ylabel('Clique size')
    bottomFig.tight_layout()
    plt.savefig(bottomSavePath,format='pdf',bbox_inches='tight')
    
    topFig = plt.figure()
    topAx1 = topFig.add_subplot(221)
    topAx1.plot(topIndices,richnesses,marker=scatterMarker,ls='')
    topAx1.set_xlabel('Number of top nodes (companies)')
    topAx1.set_ylabel('Richness')
    topAx2 = topFig.add_subplot(222)
    topAx2.plot(topIndices,diversities,marker=scatterMarker,ls='')
    topAx2.set_xlabel('Number of top nodes (companies)')
    topAx2.set_ylabel('Effective diversity')
    topAx3 = topFig.add_subplot(223)
    relativeDiversity = np.array(diversities)/np.array(richnesses)
    topAx3.plot(topIndices,relativeDiversity,marker=scatterMarker,ls='')
    topAx3.set_xlabel('Number of top nodes (companies)')
    topAx3.set_ylabel('Relative diversity')
    topAx4 = topFig.add_subplot(224)
    topAx4.plot(topIndices,cliqueSizes,marker=scatterMarker,ls='')
    topAx4.set_xlabel('Number of top nodes (companies)')
    topAx4.set_ylabel('Clique size')
    topFig.tight_layout()
    plt.savefig(topSavePath,format='pdf',bbox_inches='tight')
    
        
    
    
    
    
    
    
# Accessories:
    
def getDistribution(data, nBins):
    """
    Calculates the PDF of the given data
    
    Parameters:
    -----------
    data: a container of data points, e.g. list or np.array
    nBins: int, number of bins used to calculate the distribution
    
    Returns:
    --------
    pdf: np.array, PDF of the data
    binCenters: np.array, points where pdf has been calculated
    """
    count, binEdges, _ = binned_statistic(data, data, statistic='count', bins=nBins)
    pdf = count/float(np.sum(count))
    binCenters = 0.5*(binEdges[:-1]+binEdges[1:])
    
    return pdf, binCenters

def getTopAndBottom(bnet):
    """
    Finds the sets of top and bottom nodes of a bipartite network.
    
    Parameters:
    -----------
    bnet: networkx.Graph(), a bipartite
    
    Returns:
    --------
    top: set, top nodes of the bipartite
    bottom: set, bottom nodes of the bipartite
    """
    top = {n for n, d in bnet.nodes(data=True) if d['bipartite']==0}
    bottom = set(bnet) - top
    return top, bottom
    
def getJaccardIndex(a,b):
    """
    Calculates the Jaccard index of two lists (Jaccard(a,b) = |intersection(a,b)|/|union(a,b)|)
    
    Parameters:
    -----------
    a, b: two lists (or sets)
    
    Returns:
    J: Jaccard index of a and b
    """
    a = set(a)
    b = set(b)
    J = len(a.intersection(b))/float(len(a.union(b)))
    return J
    




    
        

    
    
