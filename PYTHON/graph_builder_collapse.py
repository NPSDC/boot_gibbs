import pandas as pd
import gzip
import os
import numpy as np
import readBootstraps
import networkx as nx
import itertools
import logging
import logging.handlers
import argparse
import json
import pickle as pi
from collections import defaultdict


def match(l1, l2):
    return[ l2.index(x) if x in l1 else None for x in l1 ]
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
def read_eqfile(base, exp, nEq = None):

    exp_dir = os.path.sep.join([base,exp])
    eq_file = os.path.sep.join([exp_dir,'aux_info','eq_classes.txt.gz'])

    logging.info('eq file {}'.format(eq_file))

    G = nx.Graph()
    eqClasses = {}
    tnames = []
    tnamemap = {}
    single_nodes = set()
    eqclasses_name = dict()
    rep = [] ## Equivalence class transcripts getting repeated
    t_not = set()
    
    with gzip.open(eq_file) as ifile:
        numTran = int(ifile.readline().decode('utf-8').rstrip())
        numEq = int(ifile.readline().decode('utf-8').rstrip())
        logging.info('eq file: {}; # tran = {}; # eq = {}'.format(eq_file, numTran, numEq))
        firstSamp = True ### Meaning

        if firstSamp:
            for i in range(numTran):
                tnames.append(ifile.readline().decode('utf-8').rstrip())
                tnamemap[len(tnames) - 1] = tnames[-1]
        else:
            for i in range(numTran):
                ifile.readline().decode('utf-8').rstrip()

        if(not nEq is None):
            numEq=nEq
        for i in range(numEq):
            print(i,end='\r')
            toks = ifile.readline().decode('utf-8').rstrip().split('\t')
            nt = int(toks[0])

            # the transcript
            tids = tuple(list(map(int,(toks[1:nt+1]))))
            count = int(toks[-1])
            
            flagEqExists = tids in eqClasses # if already exists, then add counts
            if not flagEqExists:
                eqClasses[tids] = count
                eqclasses_name[i] = tids
            else:
                eqClasses[tids] += count
                rep.append(tids)
                
            if len(tids) == 1:
                if G.has_node(tids[0]):
                    G.remove_node(tids[0])
                single_nodes.add(tids[0])
                t_not.add(tids[0])
                continue

            
            combos = list(
                itertools.combinations(np.arange(len(tids)).astype(int),2)
            )
            t_not_flag = 0 ## Flag for setting vertices that do not satisfy the use case
            if(not flagEqExists):
                for t in tids:
                    if t in G.nodes():
                        t_not_flag = 1
                        break
            for combo in combos:
                t1, t2 = sorted([tids[combo[0]],tids[combo[1]]])
                
                if(t1 in t_not or t2 in t_not):
                    t_not_flag = 1

                if not (t1,t2) in G.edges:
                    G.add_edge(t1, t2, count = count, eqClass = [i])
                    #seen_edge.add((t1, t2))
                else:
                    G[t1][t2]['count'] += count
                    if(not flagEqExists): ##Only if a different equivalence class
                        G[t1][t2]['eqClass'].append(i)
                        t_not_flag = 1

            if(t_not_flag == 1):
                for t in tids:
                    t_not.add(t)
            
        return tnames,tnamemap,eqClasses,G,rep,single_nodes,t_not,eqclasses_name

def extract_req_cliques(G, t_not):
    cliques = list(nx.find_cliques(G))
    clique_inds = list()
    c_others = list()
    for i in range(len(cliques)):
        cl = set(cliques[i])
        if(len(cl.difference(t_not)) == len(cl)):
            cl = list(cl)
            c_others.append(i)
            neigh = list(G.neighbors(cl[0]))
            l_neigh = len(neigh)
            req_flag = True
            for n in cl[1:]:
                if(l_neigh != len(list(G.neighbors(n)))):
                    req_flag = False
                    break
            if(req_flag):
                clique_inds.append(i)
    return cliques, clique_inds, c_others

def extract_transcripts(tnamemap, cliques, cl_inds):

    cliques = [cliques[i] for i in cl_inds]
    def get_trans(cl):
        return [tnamemap[i] for i in cl]
    return(map(get_trans, cliques))

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
    
    args = parser.parse_args()

    base = args.base
    expfile = args.exp
    outdir = args.outdir
    
    logging.basicConfig(level=logging.INFO)
    logging.info('Reading the experiments ')
    experiments = []
    with open(expfile, 'r') as fp:
        for line in fp:
            experiments += [line.strip()]

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    tnames_written = False
    for exp in experiments:
        tnames,tnamemap,eqclasses,G,rep,single_nodes,t_not,eqclass_name = read_eqfile(
            base,
            exp)
        cliques, clique_inds, c_others = extract_req_cliques(G, t_not)

        
        ds_dict = {"tnames":tnames, "tnamemap":tnamemap, "eqclasses":eqclasses, "rep":rep, "single_nodes":single_nodes, 
        "t_not":t_not, "eqclass_name":eqclass_name, "cliques":cliques, "clique_inds":clique_inds,"c_others":c_others}

        graph_file = os.path.sep.join([outdir, exp+'.G.pk'])
        others = os.path.sep.join([outdir, exp + '_oth.pi'])
        pi.dump(ds_dict, open(others, "wb"))
        nx.write_gpickle(G, graph_file)


if __name__ == "__main__":
    main()

