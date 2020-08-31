import pandas as pd
import os
import numpy as np
import readBootstraps
import networkx as nx
import itertools
import logging
import logging.handlers
import argparse
import json


import pandas as pd
import os
import numpy as np
import readBootstraps
import networkx as nx
import itertools
import logging
import logging.handlers
import argparse
import json
from scipy.stats.stats import pearsonr
from collections import defaultdict


# Read the correlation values between the transcripts
def bootstrap_reader(base ,experiments):
    logging.basicConfig(level=logging.INFO)
    logging.info('Reading bootstrap files from {} directory'.format(base))
    bootstrap_df_vec = {}
    for exp in experiments:
        txp_names , bootstrapList = readBootstraps.getBootstraps(
            '{}/{}/'.format(base,exp)
        )
        gibbs_df_salmon = pd.DataFrame(bootstrapList).T
        gibbs_df_salmon.index = txp_names
        bootstrap_df_vec[exp] = gibbs_df_salmon
    return bootstrap_df_vec


# Make graph from the equivalence classes
def read_eqfle(base,
               exp,
               bootstrap_df,
               l = 0,
               build_graph = False,
               count_isolated = False
              ):

    exp_dir = os.path.sep.join([base,exp])
    quant_file = os.path.sep.join([exp_dir,'quant.sf'])
    eq_file = os.path.sep.join([exp_dir,'aux_info','eq_classes.txt'])

    logging.info('quant file: {}, eq file {}'.format(quant_file, eq_file))

    quant_dict = pd.read_csv(
        quant_file,
        sep = '\t',
        usecols = ['Name','NumReads']
    ).set_index('Name').to_dict()['NumReads']

    G = nx.Graph()
    eqClasses = {}
    eqClassNormWeights = {}
    edgeToEqClass = defaultdict(list)
    tnames = []
    tnamemap = {}

    with open(eq_file) as ifile:
        numTran = int(ifile.readline().rstrip())
        numEq = int(ifile.readline().rstrip())
        logging.info('eq file: {}; # tran = {}; # eq = {}'.format(eq_file, numTran, numEq))
        firstSamp = True

        if firstSamp:
            for i in range(numTran):
                tnames.append(ifile.readline().rstrip())
                tnamemap[tnames[-1]] = len(tnames) - 1
            diagCounts = np.zeros(len(tnames))
            sumCounts = np.zeros(len(tnames))
            ambigCounts = np.zeros(len(tnames))
        else:
            for i in range(numTran):
                ifile.readline()

        epsilon =  np.finfo(float).eps

        single_nodes = set()
        autocorr_map = {}
        count_sum_map = {}
        sum_weight = 0
        seen_edge = set()

        for i in range(numEq):
            print(i,end='\r')
            toks = ifile.readline().rstrip().split('\t')
            nt = int(toks[0])


            # the transcript
            tids = tuple(list(map(int,(toks[1:nt+1]))))
            # convertedtids = [transcriptNameMap[tnames[tid]] for tid in tids]
            weights = tuple(list(map(float,toks[nt+1:-1])))
            count = int(toks[-1])

            if tids + weights in eqClasses:
                eqClasses[tids + weights] += count
            else:
                eqClasses[tids + weights] = count

            if len(tids) == 1 and not(count_isolated):
                if G.has_node(tids[0]):
                    G.remove_node(tids[0])
                single_nodes.add(tids[0])
                continue

            normWeight = 0
            if build_graph:
                combos = list(
                    itertools.combinations(np.arange(len(tids)).astype(int),2)
                )
                for combo in combos:
                    t1, t2 = tids[combo[0]],tids[combo[1]]
                    if (not((t1 in single_nodes) or
                        (t2 in single_nodes)) or
                        count_isolated
                    ):
                        w = (weights[combo[0]] * quant_dict[tnames[t1]]
                            +
                            weights[combo[1]] * quant_dict[tnames[t2]]
                        )

                        if w != 0:
                            w = np.log(w)
                        edgeToEqClass[(t2,t2)].append(i)
                        if not (t1,t2) in autocorr_map:
                            autocorr_map[(t1,t2)] = pearsonr(
                                bootstrap_df.loc[tnames[t1]].values,
                                bootstrap_df.loc[tnames[t2]].values,
                            )[0]
                            count_sum_map[(t1,t2)] = (
                                quant_dict[tnames[t1]] +
                                quant_dict[tnames[t2]]
                            )

                        if not (t1,t2) in seen_edge:
                            # print(combo)
                            G.add_edge(t1, t2, weight= w, count = count)
                            G[t1][t2]['corr'] = autocorr_map[(t1,t2)]
                            G[t1][t2]['weight'] = count*w
                        else:
                            G[t1][t2]['weight'] +=  count*w
                            G[t1][t2]['count'] += count

                        normWeight += count*w
                eqClassNormWeights[i] = normWeight

    # normalize
    if build_graph:
        for (u,v) in G.edges():
            for i in edgeToEqClass[(u,v)]:
                G[u][v]['weight'] = G[u][v]['weight'] / eqClassNormWeight[i]
    return tnames,tnamemap,eqClasses,G



def main():
    parser = argparse.ArgumentParser(
        description='Construct the graph from equivalence classes and weights'
        )
    parser.add_argument(
        '-b', '--base', required=True, type=str, help='base directory that contains the experiments')
    parser.add_argument(
        '-e', '--exp', required=True, type=str, help='file that contains experiment names')
    parser.add_argument(
        '-o', '--outdir', required=True, type=str, help='directory to store the networkx graphs')
    parser.add_argument(
        '--merge', required=False,  default = False, action='store_true', help='merge the graphs')
    parser.add_argument(
        '--keepIsolated', required=False,  default = False, action='store_true', help='keep the isolated nodes')
    parser.add_argument(
        '-l', required=False, type=float, default=0.5, help='regularizer')

    args = parser.parse_args()

    base = args.base
    expfile = args.exp
    outdir = args.outdir
    merge = args.merge
    keepIsolated = args.keepIsolated
    reg = args.l

    logging.basicConfig(level=logging.INFO)
    logging.info('Reading the experiments ')
    experiments = []
    with open(expfile, 'r') as fp:
        for line in fp:
            experiments += [line.strip()]

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    bootstrap_df_vec = bootstrap_reader(base, experiments)
    tnames_written = False
    for exp in experiments:
        tnames,tnamemap,eqclasses,G = read_eqfle(
            base,
            exp,
            bootstrap_df_vec[exp],
            l = reg,
            build_graph = True,
            count_isolated = keepIsolated
        )
        graph_file = os.path.sep.join([outdir, exp+'.G.pk'])
        nx.write_gpickle(G, graph_file)
        if not tnames_written:
            tname_file = os.path.sep.join([outdir, 'tnames.list'])
            with open(tname_file, 'w') as fp:
                for tname in tnames:
                    fp.write('{}\n'.format(tname))

            tnames_written = True



if __name__ == "__main__":
    main()

