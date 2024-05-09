import networkx as nx
import numpy as np
import operator
from pquality.PartitionQuality import conductance, triangle_participation_ratio
import matplotlib.pyplot as plt

def approx_page_rank(G, seed, beta, epsilon):
	r = np.zeros((max(G.nodes()), ), dtype="float32")
	q = r.copy()
	q[seed - 1] = 1.

	candidates = [seed]
	track = {seed : 1}
	track_r = {}

	while len(candidates) > 0:
		#Choose any vertex u where q_u / d_u >= epsilon
		u = candidates[np.random.randint(len(candidates))]

		#Push u, r, q
		r_prime = r.copy()
		q_prime = q.copy()

		r_prime[u - 1] = r[u - 1] + (1 - beta) * q[u - 1]
		track_r[u] = r_prime[u - 1]
		q_prime[u - 1] = 0.5 * beta * q[u - 1]

		for v in G[u].keys():
			q_prime[v - 1] = q[v - 1] + 0.5 * beta * (q[u - 1]/G.degree(u))
			track[v] = q_prime[v - 1]

		r = r_prime.copy()
		q = q_prime.copy()
		candidates = [k for k in track.keys() if q[k - 1] / G.degree(k) >= epsilon]
	return r, track_r

def detect_community(G, seed, beta, epsilon, alpha):
    r, track_r = approx_page_rank(G, seed, beta, epsilon)
    sorted_r = sorted(track_r.items(), key=operator.itemgetter(1), reverse=True)

    k_star = []
    k = 0

    scores = []
    candidate_k = None
    highest_score = 0

    while k < len(sorted_r):
        k += 1
        c = G.subgraph(map(lambda t: t[0], sorted_r[:k]))
        curr_score = conductance(G, c)
        scores.append(curr_score)
        if curr_score > highest_score:
            highest_score = curr_score
        if candidate_k is None and k > 1 and curr_score > scores[k-2] and curr_score  * alpha < highest_score:
            candidate_k = k - 1
        elif candidate_k is not None and curr_score > alpha * scores[candidate_k - 1]:
            k_star.append(candidate_k)
            candidate_k = None
            highest_score = scores[k-1]
        elif candidate_k is not None and scores[k-1] < scores[candidate_k - 1]:
            candidate_k = None

    plt.plot(range(1, k+1), scores)
    for k in k_star:
        plt.plot(k, scores[k-1], 'ro')
    plt.xlabel('k')
    plt.ylabel('Score')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Score vs. k')
    plt.show()

    return [G.subgraph(map(lambda t: t[0], sorted_r[:k])) for k in k_star]

def detect_community_tpr(G, seed, beta, epsilon, alpha):
    r, track_r = approx_page_rank(G, seed, beta, epsilon)
    sorted_r = sorted(track_r.items(), key=operator.itemgetter(1), reverse=True)

    k_star = []
    k = 0

    scores = []
    candidate_k = None
    highest_score = 0

    while k < len(sorted_r):
        k += 1
        c = G.subgraph(map(lambda t: t[0], sorted_r[:k]))
        curr_score = triangle_participation_ratio(c)
        scores.append(curr_score)
        if curr_score > highest_score:
            highest_score = curr_score
        if candidate_k is None and k > 1 and curr_score > scores[k-2] and curr_score  * alpha < highest_score:
            candidate_k = k - 1
        elif candidate_k is not None and curr_score > alpha * scores[candidate_k - 1]:
            k_star.append(candidate_k)
            candidate_k = None
            highest_score = scores[k-1]
        elif candidate_k is not None and scores[k-1] < scores[candidate_k - 1]:
            candidate_k = None

    plt.plot(range(1, k+1), scores)
    for k in k_star:
        plt.plot(k, scores[k-1], 'ro')
    plt.xlabel('k')
    plt.ylabel('Score')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Score vs. k')
    plt.show()

    return [G.subgraph(map(lambda t: t[0], sorted_r[:k])) for k in k_star]

def detect_community_both(G, seed, beta, epsilon, alpha, plot=True):
    r, track_r = approx_page_rank(G, seed, beta, epsilon)
    sorted_r = sorted(track_r.items(), key=operator.itemgetter(1), reverse=True)

    k_star_conductance = []
    k_star_tpr = []
    k = 0

    scores_conductance = []
    scores_tpr = []
    candidate_k_conductance = None
    candidate_k_tpr = None
    highest_score_conductance = 0
    highest_score_tpr = 0

    while k < len(sorted_r):
        k += 1
        c = G.subgraph(map(lambda t: t[0], sorted_r[:k]))
        curr_score_conductance = conductance(G, c)
        curr_score_tpr = triangle_participation_ratio(c)
        scores_conductance.append(curr_score_conductance)
        if k < 3:
            scores_tpr.append(1)
        else:
            scores_tpr.append(curr_score_tpr)
        if curr_score_conductance > highest_score_conductance:
            highest_score_conductance = curr_score_conductance
        if curr_score_tpr > highest_score_tpr:
            highest_score_tpr = curr_score_tpr
        if candidate_k_conductance is None and k > 1 and curr_score_conductance > scores_conductance[k-2] and curr_score_conductance * alpha < highest_score_conductance:
            candidate_k_conductance = k - 1
        elif candidate_k_conductance is not None and curr_score_conductance > alpha * scores_conductance[candidate_k_conductance - 1]:
            k_star_conductance.append(candidate_k_conductance)
            candidate_k_conductance = None
            highest_score_conductance = scores_conductance[k-1]
        elif candidate_k_conductance is not None and scores_conductance[k-1] < scores_conductance[candidate_k_conductance - 1]:
            candidate_k_conductance = None

        if candidate_k_tpr is None and k > 1 and curr_score_tpr > scores_tpr[k-2] and curr_score_tpr * alpha < highest_score_tpr:
            candidate_k_tpr = k - 1
        elif candidate_k_tpr is not None and curr_score_tpr > alpha * scores_tpr[candidate_k_tpr - 1]:
            k_star_tpr.append(candidate_k_tpr)
            candidate_k_tpr = None
            highest_score_tpr = scores_tpr[k-1]
        elif candidate_k_tpr is not None and scores_tpr[k-1] < scores_tpr[candidate_k_tpr - 1]:
            candidate_k_tpr = None
    if plot:
        plt.plot(range(1, k+1), scores_conductance, label='Conductance')
        plt.plot(range(1, k+1), scores_tpr, label='TPR')
        plt.plot(k_star_conductance, list(scores_conductance[k-1] for k in k_star_conductance), 'ro', label='Conductance k*')
        plt.plot(k_star_tpr, list(scores_tpr[k-1] for k in k_star_tpr), 'go', label='TPR k*')
        plt.xlabel('k')
        plt.ylabel('Score')
        plt.xscale('log')
        plt.yscale('log')
        plt.title('Score vs. k')
        plt.legend()
        plt.show()

    return [G.subgraph(map(lambda t: t[0], sorted_r[:k])) for k in k_star_conductance], [G.subgraph(map(lambda t: t[0], sorted_r[:k])) for k in k_star_tpr]

#c1 ground truth, c2 predicted

def evaluate_f1(c1_nodes, c2_nodes):

	relv = len(set(c1_nodes).intersection(c2_nodes))
	nrelv = len(c1_nodes) - relv
	irelv = len(c2_nodes) - relv
	
	return 2 * relv / (2 * relv + nrelv + irelv)