import algorithms.spectral
from igraph import Graph, VertexCover
import uuid

def __get_weight_attr(G, metric_name, weights):
    '''
    :G graph
    :metric_name the name of the metric calling this function
    :weights

    return (weight_attr, remove) where weight_attr is the weight attribute name, and "remove" is a boolean
        as to wether to get rid of the weight attribute. This happens if we create the weight attribute for
        the purpose of the metric

    '''

    #if the weights parameter is a string then the graph utilizes weights
    if isinstance(weights, str):
      return (weights, False)
    #if the weights is being used for something else, then we
    elif weights is not None:
      attr_name = uuid.uuid5(uuid.NAMESPACE_DNS, '{}.circulo.lab41'.format(metric_name))
      G.es[attr_name] = weights
      return (attr_name, True)
    return (None, False)

def __remove_weight_attr(G, attr_name, remove):
    if remove:
      del G.es[uid]

def __weighted_sum(external_edges, w_attr):
  return len(external_edges) if w_attr is None else sum([ e[w_attr] for e in external_edges ])

def external_edges(cover):
    array_of_sets = [ [] for v in cover ]

    # Go through membership array and find nodes in each community
    community_sets = {}
    membership_arrays = cover.membership
    for vertex_id, vertex_communities in enumerate(membership_arrays):
        for community in vertex_communities:
            if community not in community_sets:
                community_sets[community] = set()
            community_sets[community].add(vertex_id)

    # Move from dictionary of communities to array of communities
    # TODO: Maybe we don't have to do this??
    membership_sets = [community_sets[community] for community in sorted(community_sets.keys())]

    for (edge, crossing) in zip(cover.graph.es, cover.crossing()):
        if crossing:
            src,dst = edge.tuple
            for i, membership_set in enumerate(membership_sets):
                #print(i)
                if src in membership_set and dst not in membership_set:
                    array_of_sets[i].append(edge)
                elif dst in membership_set and src not in membership_set:
                    array_of_sets[i].append(edge)

    return array_of_sets


def conductance(cover, weights=None, allow_nan=False):
    '''
    Conductance is the ratio between the (weighted) number of external (boundary) edges in a cluster and the cluster's total (weighted) number of edges
    '''
    w_attr, remove = __get_weight_attr(cover.graph, 'conductance', weights)

    mode = "nan" if allow_nan else 0
    rv = []
    external_edges = cover.external_edges()
    for i in range(len(cover)):
        int_edges_cnt = __weighted_sum(cover.subgraph(i).es(), w_attr)
        ext_edges_cnt = __weighted_sum(external_edges[i], w_attr)
        denominator = (2.0*int_edges_cnt+ext_edges_cnt)

        rv += [ext_edges_cnt/denominator if denominator > 0 else float(mode)]


    __remove_weight_attr(cover.graph, w_attr, remove)
    return rv

def min_conductance(G, weights=None, tries=3):
    '''
    Returns the minimum conductance of a Graph by using spectral clustering to ``approximate'' the minimum ratio-cut.
    http://www.kyb.mpg.de/fileadmin/user_upload/files/publications/attachments/Luxburg07_tutorial_4488%5b0%5d.pdf
    '''
    (rv_val, rv_vc) = (float("inf"), None)
    for i in range(0,tries):
        #Obtain a cut of G, it should already be a minimum
        curr_vc = G.community_spectral(k=2, weights=weights, which='NCut')
        curr_val = max(curr_vc.as_cover().conductance())
        if curr_val < rv_val :
            (rv_val, rv_vc) = (curr_val, curr_vc)


    return rv_val, rv_vc

def separability(cover, weights=None, allow_nan = False):
    '''
    Separability is the ratio between the (weighted) number of internal edges in a cluster and its (weighted) number of external (boundary) edges.
    '''

    mode = "nan" if allow_nan else 0
    w_attr, remove = __get_weight_attr(cover.graph, 'separability', weights)
    rv = []
    external_edges = cover.external_edges()
    for i in range(len(cover)):
        int_edges_cnt = __weighted_sum(cover.subgraph(i).es(), w_attr)
        ext_edges_cnt = __weighted_sum(external_edges[i], w_attr)
        rv += [1.0*int_edges_cnt/ext_edges_cnt if ext_edges_cnt > 0 else float(mode)]

    __remove_weight_attr(cover.graph, w_attr, remove)
    return rv


VertexCover.conductance = conductance
VertexCover.external_edges = external_edges
Graph.min_conductance = min_conductance
VertexCover.separability = separability