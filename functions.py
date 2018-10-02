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
from networkx.algorithms import bipartite
from scipy.stats import binned_statistic, binned_statistic_2d
from itertools import combinations

# Functions for reading metadata

def readNodes(cfg):
    """
    Reads information about network nodes and their attributes. To this end, the information
    must be stored in a csv file with named columns. By default, column names are set for reading
    company information ('Alias' (anonymized alias used to refer to the
    node in the network), 'Name' (real name of the company), and 'FoB' (field of business)).
    Other sets of column names can be given as parameters. Data of the first column will be used
    as node labels and must therefore be anonymized.
    
    Parameters:
    -----------
    cfg: dict, contains:
        inputPath: str, path to the input csv file
        columnNames: list of strs, column names in the csv file (default ['Alias', 'Name', 'FoB'])
        
    Returns:
    --------
    nodeInfo: 2D list, first element contains node aliases and further elements additional attributes
    columnNames: list of strs, column names of the csv file (and therefore titles of the lists in nodeInfo)
    
    """
    inputPath = cfg['inputPath']
    if 'columnNames' in cfg.keys():
        columnNames= cfg['columnNames']
    else:
        columnNames = ['Alias', 'Name', 'FoB']
    
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
        inputPath: str, path to the metadata input csv file
        columnNames: list of strs, column names in the csv file
        net: nx.Graph()
        bipartite: int (0 or 1), bipartite node class
        tags: list of strs, the field of business tags
        
    Returns:
    --------
        net: nx.Graph()
    """
    topColors = cfg['topColors']
    networkColors = cfg['networkColors']
    tags = cfg['tags']
    # Reading node information
    nodeInfo, columnNames = readNodes(cfg)
    aliases = nodeInfo[0] # first column is used as node alias
    attributes = nodeInfo[1::] # additional columns contain attributes
    
    # Adding nodes and updating their attributes
    for i, node in enumerate(aliases):
        print node
        net.add_node(node, bipartite=bipartiteClass)
        attr_dic = {}
        for columnName, attribute in zip(columnNames[1::], attributes):
            attr_dic[columnName] = attribute[i]
        if topColors:
            foundTags = []
            for tag in tags:
                if tag in node:
                    foundTags.append(tag)
            foundLengths = [len(found) for found in foundTags]
            foundTag = foundTags[foundLengths.index(max(foundLengths))]
            attr_dic['nodeColor'] = networkColors(tags.index(foundTag))
            attr_dic['tag'] = foundTag
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
        companyInputPath: str, path to the company metadata input csv file
        eventInputPath: str, path to the event metadata input csv file
        linkInputPath: str, path to the link input csv file
        companyColumnNames: list of strs, column names in the csv file containing the company information (default ['Alias', 'Name', 'FoB'])
        eventColumnNames: list of strs, column names in the csv file containing the company information (default ['Alias', 'Name'])
        linkColumnNames: list of strs, column names in the link csv file (default ['Source', 'Target'])
        topColors: bln, shall the node color attributes be set according to the field of business indicates in the alias? If False, a separately set networkBottomColor is used instead.
        networkColors: matplotlib.cmap, colors associated with different fields of business
        tags: list of strs, the field of business tags
        networkBottomColor: str, the color for bottom (event) nodes
        
    Returns:
    --------
    bnet: networkx.Graph(), a bipartite
    """
    bnet = nx.Graph()
    
    # Reading company information, adding company nodes
    cfg['inputPath'] = cfg['companyInputPath']
    cfg['topColors'] = True
    if 'companyColumnNames' in cfg.keys():
        cfg['columnNames'] = cfg['companyColumnNames']
    else:
        cfg['columnNames'] = ['Alias', 'Name', 'FoB']
    bnet = addNodes(cfg, bnet, bipartiteClass=0)   
    
    # Reading event information, adding event nodes
    cfg['inputPath'] = cfg['eventInputPath']
    cfg['topColors'] = False
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
    """
    nodesToRemove = []
    
    for node in bnet.nodes():
        if bnet.degree(node) == 0:
            nodesToRemove.append(node)
    
    bnet.remove_nodes_from(nodesToRemove)
    
    return bnet
    
# Basic analysis
def getDegreeDistributions(bnet, cfg):
    """
    Calculates the top and bottom (companies and events) degree distributions of
    a bipartite network, plots the distributions, and saves the plot.
    
    Parameters:
    -----------
    cfg: dict, contains:
        savePathBase: str, a base path (e.g. to a shared folder) for saving figures
        degreeSaveName: str, name of the file where to save the degree distribution plots
        nDegreeBins: int, number of bins used to calculate the distributions
        TODO: add a possibility for different number of bins in top and bottom
        topColor: str, color for plotting the top degree distribution
        bottomColor: str, color for plotting the bottom degree distribution
    bnet: networkx.Graph(), bipartite
        
    Returns:
    --------
    no direct output, saves the degree distributions in a file
    """
    nBins = cfg['nDegreeBins']
    topColor = cfg['topColor']
    bottomColor = cfg['bottomColor']   
    savePath = cfg['savePathBase'] + cfg['degreeSaveName']
    
    top = {n for n, d in bnet.nodes(data=True) if d['bipartite']==0}
    bottom = set(bnet) - top
    topDegree = nx.degree(bnet,top)
    bottomDegree = nx.degree(bnet, bottom)
    topDegree = dict(topDegree).values() # bipartite.degrees returns a DegreeView so this is required for accessing values
    bottomDegree = dict(bottomDegree).values()
    topPdf, topBinCenters = getDistribution(topDegree, nBins)
    bottomPdf, bottomBinCenters = getDistribution(bottomDegree, nBins)
    
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

def getDensity(bnet):
    """
    Returns density of the bipartite graph.
    
    Parameters:
    -----------
    bnet: networkx.Graph(), bipartite
    
    Returns:
    d: float, density of bnet
    """
    top = {n for n, d in bnet.nodes(data=True) if d['bipartite']==0}
    d = bipartite.density(bnet,top)
    return d

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
    top = {n for n, d in bnet.nodes(data=True) if d['bipartite']==0}
    bottom = set(bnet) - top
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
    top = {n for n, d in bnet.nodes(data=True) if d['bipartite']==0}
    bottom = set(bnet) - top
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
    top = {n for n, d in bnet.nodes(data=True) if d['bipartite']==0}
    bottom = set(bnet) - top
    
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
    
def pruneStars(bnet, cliques, cliqueInfo):
    """
    Removes from bistars all top nodes (companies) that participate in any other
    cliques besides the bistar. Bistar is a clique where one bottom node (event)
    is surrounded by multiple top nodes (companies).
    
    Parameters:
    -----------
    bnet: networx.Graph(), a bipartite
    cliques: list of lists, cliques of bnet; each clique is represented as a 
             list containing the nodes of the clique
    cliqueInfo: list of dicts; each dict contains one clique separated to top and
                bottom nodes (keys: 'topNodes': top nodes (companies) of the clique, 
                                    'bottomNodes' : bottom nodes (events) of the clique,
                                    'topIndex': number of top nodes of the clique,
                                    'bottomIndex': number of bottom nodes of the clique,
                                    'isBridge': does the clique consist of one top node (company) connecting multiple bottom nodes,
                                    'isStar': does the clique consist of node bottom node (event) connecting multiple top nodes)
    
    Returns:
    --------
    cliques: list of lists, the original set of cliques with bistars are replaced
             with the pruned versions
    cliqueInfo: list of dicts, updated info for the updated clique list
    """
    cliquesToRemove = []
    infosToRemove = []
    newCliques = []
    newInfos = []
    for clique, info in zip(cliques,cliqueInfo):
        if info['isStar'] == True:
            cliquesToRemove.append(clique)
            infosToRemove.append(info)
            topToRemove = []
            for top in info['topNodes']:
                if nx.degree(bnet,top) > 1:
                    topToRemove.append(top)
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
    print xEdges[0],xEdges[-1],yEdges[0],yEdges[-1]
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
        
            
            
            
    
    
# Visualization
    
def drawNetwork(bnet, cfg):
    """
    Visualizes the bipartite network and saves the figure as pdf.
    
    Parameters:
    -----------
    bnet: networkx.Graph(), the bipartite to be visualized
    cfg: dict, contains:
        tags: list of strs, the field of business tags
        networkColors: matplotlib.cmap, colors associated with different fields of business
        bottomNetworkColor: str, color of the bottom nodes (events)
        nodeSize: int, size of nodes in the visualization
        edgeWidth: double, width of network edges
        savePathBase: str, a base path (e.g. to a shared folder) for saving figures
        networkSaveName: str, name of the file where to save the network visualization
        
    Returns:
    --------
    no direct output, saves the network visualization as pdf
    """
    
    tags = cfg['tags']    
    networkColors = cfg['networkColors']
    bottomColor = cfg['networkBottomColor']
    nodeSize = cfg['nodeSize']
    edgeWidth = cfg['edgeWidth']
    
    top = {n for n, d in bnet.nodes(data=True) if d['bipartite']==0}
    bottom = set(bnet) - top
    
    pos = nx.spring_layout(bnet)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    nx.draw_networkx_nodes(bnet,pos,ax=ax,nodelist=bottom,node_color=bottomColor,node_size=nodeSize,label='Events)')
    
    for i, tag in enumerate(tags):
        taggedNodes = []
        for topn in top:
            if bnet.nodes(data=True)[topn]['tag'] == tag:
                taggedNodes.append(topn)
        nx.draw_networkx_nodes(bnet,pos,ax=ax,nodelist=taggedNodes,node_color=networkColors(i),node_size=nodeSize,label=tag)
    
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
    
    #top = {n for n, d in bnet.nodes(data=True) if d['bipartite']==0}
    #bottom = set(bnet) - top
    
    pos = nx.spring_layout(bnet)
    
    for i, clique in enumerate(cliqueInfo):
        topCliqueNodes = clique['topNodes']#set(clique) & top
        bottomCliqueNodes = clique['bottomNodes']#set(clique) & bottom
        nonCliqueNodes = set(bnet) - topCliqueNodes.intersection(bottomCliqueNodes)#set(clique)
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        nx.draw_networkx_nodes(bnet,pos=pos,nodelist=topCliqueNodes,node_size=nodeSize,node_color=topColor)
        nx.draw_networkx_nodes(bnet,pos=pos,nodelist=bottomCliqueNodes,node_size=nodeSize,node_color=bottomColor)
        nx.draw_networkx_nodes(bnet,pos=pos,nodelist=nonCliqueNodes,node_size=nodeSize,node_color=nonCliqueColor,alpha=nonCliqueAlpha)
        nx.draw_networkx_edges(bnet,pos=pos,width=edgeWidth)
        plt.axis('off')  
        
        savePath = cfg['savePathBase'] + saveName + str(i) + '.pdf'
        plt.savefig(savePath,format='pdf',bbox_inches='tight')
    
    
    
    
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

    
    




    
        

    
    
