# -*- coding: utf-8 -*-
"""
Created on Fri May  4 12:24:58 2018

@author: aokorhon

This file contains functions related to the Ilmastokumppanit bipartite project. To call, either write a 
frontend script (see XX for example) or call interactively from ipython.
"""
import pandas as pd
import networkx as nx
from networkx.algorithms import bipartite


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
    cfg: dic, contains:
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
    cfg: dic, contains:
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
    cfg: dic, contains:
        inputPath: str, path to the metadata input csv file
        columnNames: list of strs, column names in the csv file
        net: nx.Graph()
        bipartite: int (0 or 1), bipartite node class
        
    Returns:
    --------
        net: nx.Graph()
    """
    # Reading node information
    nodeInfo, columnNames = readNodes(cfg)
    aliases = nodeInfo[0] # first column is used as node alias
    attributes = nodeInfo[1::] # additional columns contain attributes
    
    # Adding nodes and updating their attributes
    for i, node in enumerate(aliases):
        net.add_node(node, bipartite=bipartiteClass)
        attr_dic = {}
        for columnName, attribute in zip(columnNames[1::], attributes):
            attr_dic[columnName] = attribute[i]
        net.nodes[node].update(attr_dic)
        
    return net
    
def createBipartite(cfg):
    """
    Creates a bipartite network with two types of nodes (companies and events), node attributes (real names, fields of business, ...), 
    and links between different types of nodes.
    
    Parameters:
    -----------
    cfg: dic, contains:
        companyInputPath: str, path to the company metadata input csv file
        eventInputPath: str, path to the event metadata input csv file
        linkInputPath: str, path to the link input csv file
        companyColumnNames: list of strs, column names in the csv file containing the company information (default ['Alias', 'Name', 'FoB'])
        eventColumnNames: list of strs, column names in the csv file containing the company information (default ['Alias', 'Name'])
        linkColumnNames: list of strs, column names in the link csv file (default ['Source, 'Target'])
        
    Returns:
    --------
    bnet: networkx.Graph(), a bipartite
    """
    bnet = nx.Graph()
    
    # Reading company information, adding company nodes
    cfg['inputPath'] = cfg['companyInputPath']
    if 'companyColumnNames' in cfg.keys():
        cfg['columnNames'] = cfg['companyColumnNames']
    else:
        cfg['columnNames'] = ['Alias', 'Name', 'FoB']
    bnet = addNodes(cfg, bnet, bipartiteClass=0)    
    
    # Reading event information, adding event nodes
    cfg['inputPath'] = cfg['eventInputPath']
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
    cfg: dic, contains:
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
    topDegree, bottomDegree = bipartite.degrees(bnet, top)
    topDegree = dict(topDegree).values() # bipartite.degrees returns a DegreeView so this is required for accessing values
    bottomDegree = dict(bottomDegree).values()
    topPdf, topBinCenters = getDistribution(topDegree, nBins)
    bottomPdf, bottomBinCenters = getDistribution(bottomDegree, nBins)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(topBinCenters, topPdf, color=topColor, label='Companies (top nodes)')
    ax.plot(bottomBinCenters, bottomPdf, color=bottomColor, label='Events (bottom nodes)')
    ax.set_xlabel('Degree')
    ax.set_ylabel('PDF')
    ax.legend()
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

    
    




    
        

    
    
