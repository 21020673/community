#!/usr/bin/env python
#
# Copyright (c) 2014 In-Q-Tel, Inc/Lab41, All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import sys
import os
import argparse
import json
from operator import itemgetter
import glob

import numpy as np
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from scipy.cluster.vq import kmeans2, whiten
import pandas as pd


metric_names = [
    "Fraction Over Median Degree",
    "Triangle Participation Ratio",
    "Cut Ratio",
    "Conductance",
    "Flake Out Degree Fraction",
    "Separability",
    "Density",
    "Cohesiveness",
    "Clustering Coefficient"
    ]



#enables a single legend for all graphs rather
#than a legend for each graph
global_lines = list()


def analyze_metric_file(metrics_path, out_path):
    '''
    Analyzes a given metrics results file

    Args:
        metrics_df: path to the metrics results file
        out_path: directory where to store the results

    Return:
        array of rates of which metric is the most correlated
    '''

    job_name = os.path.splitext(os.path.basename(metrics_path))[0]
    if job_name == "dblp":
        job_name = "DBLP"
    elif job_name == "lj":
        job_name = "LiveJournal"

    #we select nine of the metrics that we wish to test
    metrics_df = pd.read_csv(metrics_path)
    l = metrics_df.values.tolist()
    k = len(l)

    fig = figure(figsize=(16, 10))
    fig.suptitle('Goodness Metrics Indicators: ' + job_name, fontsize=20)

    #iterate over each metric, treating it as the comparator
    #aggregate result is a list of rates
    rates = [run(l, feature_idx + 5, k) for feature_idx in range(4)]

    #plot the results
    global global_lines
    fig.legend(global_lines, ('FOMD', 'TPR', 'CR', 'Cond', 'FODF'), loc='right')
    plt.savefig(os.path.join(out_path, job_name + ".png"), format='png')

    return rates




def running_avg(l):
    '''
    Quick hack of a running average function

    Args:
        l: list of numbers to calculate running average on
    Returns:
        list of the running average at each index
    '''

    r = list()

    total_sum = 0.0

    for idx, item  in enumerate(l):
        total_sum += item
        r.append( total_sum/(idx+1.0))

    return r




def get_rankings(running_avgs, feature_idx):
    '''
    Args:
        running_avgs: a list of running averages, one for each metric
        feature_idx: the primary metric index

    Returns:
        list of top correlating indices for the feature index
    '''

    totals = np.zeros(len(running_avgs))

    for cross_section in list(zip(*running_avgs)):
        #this will always be the feature
        m = max(cross_section)
        diffs = [m-x for x in cross_section]
        for i, diff in enumerate(diffs):
            totals[i]+=diff

    total_max = max(totals)

    if total_max == 0:
        return totals

    #normalize
    totals = totals / float(max(totals))

    matches = [i for i, v in enumerate(totals) if v < .15 and i != feature_idx]

    #we need to get the corresponding total diff for each match so we can sort them
    l = list(zip(matches, [totals[i] for i in matches]))

    #return the sorted list of top correlated metric for primary metric
    return  [i  for i,_ in  sorted(l, key = itemgetter(1))]




def run(metrics_list, feature_idx, k):
    '''
    Creates graph depicting metrics relative to a single metrics specified by
    the feature_idx

    Args:
        metrics_list: List of lists of metrics
        feature_idx: Index of the metric being tested. Refer to the zipped list in main
        k: max size of num of communities we are examining
        metric_name: name of the metric being tested
    '''

    global global_lines

    x_values = range(k)
    num_features = len(metrics_list[0])
    plt.subplot(221 + feature_idx - 5)
    plt.xscale("log")

    vertical = list()

    #basically iterate through features, plotting each one except the primary
    for i in range(5):
        s = sorted(metrics_list, key = itemgetter(i), reverse=True)[:k]
        #now get the running average for the main metric
        running_avgs = running_avg([v[feature_idx] for v in s])
        vertical.append(running_avgs)
        #to keep colors consistent, we need to use a global list of 2D lines
        if len(global_lines) < num_features:
            line, = plt.plot(x_values, running_avgs)
            global_lines.append(line)
        else:
            plt.plot(x_values, running_avgs, color=global_lines[i].get_color())


    plt.ylabel(metric_names[feature_idx])
    plt.xlabel("Rank")

    return get_rankings(vertical, feature_idx)



def main():
    parser = argparse.ArgumentParser(description='Experiment of Correlations in Goodness Metrics')
    parser.add_argument('metrics_path', help="path to metrics results")
    parser.add_argument('--out_path', default='indicators_results',  help="path to save results")
    args = parser.parse_args()

    #create output directory if not exists
    if not os.path.exists(args.out_path):
        os.mkdir(args.out_path)


    rates_agg = list()

    num_files = 0
    for i, f in enumerate(glob.glob(args.metrics_path)):
        results = analyze_metric_file(f, args.out_path)
        rates_agg.append(results)
        num_files+=1


if __name__ == "__main__":
    main()
