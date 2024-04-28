from functools import cache
import networkx as nx
import numpy as np
import csv
import itertools
from tqdm import tqdm
import igraph as ig

__author__ = "Giulio Rossetti"
__contact__ = "giulio.rossetti@gmail.com"
__github__ = "https://github.com/GiulioRossetti"


def median(lst):
    return np.median(np.array(lst))


def community_modularity(coms, g):
    if type(g) != nx.Graph:
        raise TypeError("Bad graph type, use only non directed graph")

    inc = dict([])
    deg = dict([])
    links = g.size(weight='weight')
    if links == 0:
        raise ValueError("A graph without link has an undefined modularity")

    for node in g:
        try:
            com = coms[node]
            deg[com] = deg.get(com, 0.) + g.degree(node, weight='weight')
            for neighbor, dt in g[node].items():
                weight = dt.get("weight", 1)
                if coms[neighbor] == com:
                    if neighbor == node:
                        inc[com] = inc.get(com, 0.) + float(weight)
                    else:
                        inc[com] = inc.get(com, 0.) + float(weight) / 2.
        except:
            pass

    res = 0.
    for com in set(coms.values()):
        res += (inc.get(com, 0.) / links) - (deg.get(com, 0.) / (2.*links))**2
    return res


def internal_edge_density(coms):
    ms = coms.number_of_edges()
    ns = coms.number_of_nodes()
    try:
        internal_density = float(ms) / (float(ns * (ns - 1)) / 2)
    except:
        return 0
    return internal_density


def edges_inside(coms):
    return coms.number_of_edges()


def average_internal_degree(coms):
    ms = coms.number_of_edges()
    ns = coms.number_of_nodes()
    try:
        avg_id = float(2*ms) / ns
    except:
        return 0
    return avg_id


def fraction_over_median_degree(coms):
    ns = coms.number_of_nodes()
    degs = coms.degree()

    med = median([d[1] for d in degs])
    above_med = len([d[0] for d in degs if d[1] > med])
    try:
        ratio = float(above_med) / ns
    except:
        return 0
    return ratio


def expansion(g, coms):
    ns = coms.number_of_nodes()
    cs = len(tuple(nx.edge_boundary(g, coms)))
    try:
        return float(cs) / ns
    except:
        return 0


def cut_ratio(g, coms):
    ns = coms.number_of_nodes()
    cs = len(tuple(nx.edge_boundary(g, coms)))
    try:
        ratio = float(cs) / (ns * (len(g) - ns))
    except:
        return 0
    return ratio


def conductance(g, coms):
    ms = coms.number_of_edges()
    cs = len(tuple(nx.edge_boundary(g, coms)))
    try:
        ratio = float(cs) / ((2 * ms) + cs)
    except:
        return 0
    return ratio


def normalized_cut(g, coms):
    ms = coms.number_of_edges()
    cs = len(tuple(nx.edge_boundary(g, coms)))
    try:
        ratio = (float(cs) / ((2 * ms) + cs)) + \
            float(cs) / (2 * (g.number_of_edges - ms) + cs)
    except:
        return 0

    return ratio


def __out_degree_fraction(g, coms):
    nds = []
    for n in coms:
        nds.append(g.degree(n) - coms.degree(n))
    return nds


def max_odf(g, coms):
    return max(__out_degree_fraction(g, coms))


def avg_odf(g, coms):
    return float(sum(__out_degree_fraction(g, coms)))/len(coms)


def flake_odf(g, coms):
    df = 0
    for n in coms:
        fr = coms.degree(n) - (g.degree(n) - coms.degree(n))
        if fr < 0:
            df += 1
    return float(df)/len(coms)


def triangle_participation_ratio(coms):
    cls = nx.triangles(coms)
    nc = [n for n in cls if cls[n] > 0]
    return float(len(nc))/len(coms)


def modularity(g, coms):
    part = {}
    ids = 0
    for c in coms:
        for n in c:
            part[n] = ids
        ids += 1

    mod = community_modularity(part, g)
    return mod

def separability(g, coms):
    ms = coms.number_of_edges()
    cs = len(tuple(nx.edge_boundary(g, coms)))
    try:
        ratio = float(ms) / cs
    except:
        return 1
    return ratio

def clustering_coefficient(coms):
    return nx.average_clustering(coms)

def cohesiveness(G, weights=None):
    '''
    Equation: g(S) = minS′⊂S φ(S′) where φ(S′) is the conductance of S′ measured in the induced subgraph by S.
    To iterate over all possible subgraphs of a community would be too inefficient 2^n, therefore we approximate
    the best subgraph (which would have the lowest conductance) by using Local Spectral communitying to find the best
    cut
    (cite: http://cs.stanford.edu/people/jure/pubs/comscore-icdm12.pdf)
    '''
    import algorithms.min_conductance
    if G.vcount() <= 2:
        val = 1
    else:
        #TODO: Consider using G_i.mincut() instead.
        val, vc = G.min_conductance(weights=weights)
    return val

def pquality_summary(graph, partition, output_file):
    """
    Compute community quality metrics and output to a CSV file.

    Parameters:
    - graph: NetworkX graph object
    - partition: List of community node sets
    - output_file: Path to output CSV file
    """
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([
            "FOMD",
            "TPR",
            "Cut Ratio",
            "Conductance",
            "Flake-ODF",
            "Separability",
            "Density",
            "Cohesiveness",
            "Clustering Coefficient"
        ])

        for community_nodes in tqdm(partition):
            community_graph = graph.subgraph(community_nodes)
            ig_graph = ig.Graph.from_networkx(community_graph)

            ied = internal_edge_density(community_graph)
            fomd = fraction_over_median_degree(community_graph)
            cr = cut_ratio(graph, community_graph)
            cond = conductance(graph, community_graph)
            flake = flake_odf(graph, community_graph)
            tpr = triangle_participation_ratio(community_graph)
            sep = separability(graph, community_graph)
            coh = cohesiveness(ig_graph)
            cc = clustering_coefficient(community_graph)

            writer.writerow([
                fomd,
                tpr,
                cr,
                cond,
                flake,
                sep,
                ied,
                coh,
                cc
            ])

