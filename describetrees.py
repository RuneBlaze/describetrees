import argparse
import dendropy
from dendropy.calculate.treecompare import false_positives_and_negatives
import argparse
import toml
import logging
import pathlib
import os
from os.path import join, exists
from glob import glob
from fnmatch import fnmatch
import progressbar
import numpy as np
from joblib import Parallel, delayed

parser = argparse.ArgumentParser(description='describe trees')
parser.add_argument('file', type=str)
parser.add_argument('-l', '--level', type=int, default=30)
parser.add_argument('--on', type=str, required=True)
parser.add_argument('--ref', type=str, required=True)
parser.add_argument('--against', type=str, required=True)
parser.add_argument('-c', '--cores', type=int, default=-1)

args = parser.parse_args()

logging.basicConfig(level=args.level)

DATA = toml.load(args.file)
WORKINGDIR = pathlib.Path(args.file).parent.resolve()
CORES = args.cores
datasets = DATA["datasets"]

def compare_trees(tr1, tr2):
    # Find leaf labels that are in both trees
    lb1 = set([l.taxon.label for l in tr1.leaf_nodes()])
    lb2 = set([l.taxon.label for l in tr2.leaf_nodes()])
    com = lb1.intersection(lb2)
    # Restrict trees to shared leaf set
    if com != lb1 or com != lb2:
        com = list(com)
        tns = dendropy.TaxonNamespace(com)
        tr1.retain_taxa_with_labels(com)
        tr1.migrate_taxon_namespace(tns)
        tr2.retain_taxa_with_labels(com)
        tr2.migrate_taxon_namespace(tns)
    com = list(com)
    # Update tree bipartitions
    tr1.update_bipartitions()
    tr2.update_bipartitions()
    # Compute number of leaves and number of internal edges
    nl = len(com)
    ei1 = len(tr1.internal_edges(exclude_seed_edge=True))
    ei2 = len(tr2.internal_edges(exclude_seed_edge=True))
    # Compute number of false positives and false negatives
    [fp, fn] = false_positives_and_negatives(tr1, tr2)
    # Compute symmetric difference rate
    sd = float(fp + fn) / (ei1 + ei2)
    # Compute Robinson-Foulds error rate
    rf = float(fp + fn) / (2 * nl - 6)
    # print((nl, ei1, ei2, fp, ei2, fn, ei1, sd, rf))
    # if ei2 == 0:
    #     print(tr1, tr2)
    return(nl, ei1, ei2, fp, fn, sd, rf)

def comparetreestr(i, tr1, tr2):
    tax = dendropy.TaxonNamespace()
    tr1 = dendropy.Tree.get(data=tr1,
                            schema='newick',
                            rooting='force-unrooted',
                            taxon_namespace=tax)
    # import treeswift as ts
    # tr2ts = ts.read_tree_newick(fh.read())
    # for n in tr2ts.traverse_postorder(True, True):
    #     if n.label and "_" in n.label:
    #         n.label = n.label.split("_")[0]
    tr2 = dendropy.Tree.get(data=tr2,
                            schema='newick',
                            rooting='force-unrooted',
                            taxon_namespace=tax)
    # Unroot trees
    tr1.collapse_basal_bifurcation(set_as_unrooted_tree=True)
    tr2.collapse_basal_bifurcation(set_as_unrooted_tree=True)
    # Compute RF distance
    [nl, ei1, ei2, fp, fn, sd, rf] = compare_trees(tr1, tr2)
    # if ei1 * ei2 == 0:
    #     print(i)
    #     print(tr1)
    #     print(tr2)
    #     exit()
    return [nl, ei1, ei2, fp, fn, sd, rf]

def rf_of_newicks(i, tr1, tr2):
    return comparetreestr(i, tr1, tr2)[-2]

def flatmap(func, *iterable):
    import itertools
    return itertools.chain.from_iterable(map(func, *iterable))

def compare_from_paths(refpath, cmppath):
    logging.debug(f'refpath: {refpath}, cmptree: {cmppath}')
    with open(refpath, "r") as fh: reflines = fh.readlines()
    with open(cmppath, "r") as fh: cmplines = fh.readlines()
    if len(reflines) > 1:
        Parallel(n_jobs=CORES)(delayed(rf_of_newicks)(r, c) for i, (r, c) in enumerate(zip(reflines, cmplines)))
    else:
        return Parallel(n_jobs=CORES)(delayed(rf_of_newicks)(i, reflines[0], c) for i, c in enumerate(cmplines))

for k in datasets:
    if fnmatch(k, args.on):
        ds = datasets[k]
        logging.info(f'dealing with dataset {k}')
        datasetdir = join(WORKINGDIR, ds["basedir"])
        reftrees = []
        cmptrees = []
        for rep in glob(join(datasetdir, ds["replicates"])):
            reftrees.append(join(rep, ds["trees"][args.ref]))
            cmptrees.append(join(rep, ds["trees"][args.against]))
            assert exists(reftrees[-1])
            assert exists(cmptrees[-1])
            logging.debug(f'reftree: {reftrees[-1]}, cmptree: {cmptrees[-1]}')
        cmp_from_path = lambda x: compare_from_paths(x[0], x[1])
        rfs = []
        for i in progressbar.progressbar(range(len(reftrees))):
            l = reftrees[i]
            r = cmptrees[i]
        # for l, r in progressbar.progressbar(zip(reftrees, cmptrees)):
            rfs += compare_from_paths(l, r)
        print(f"for {k}, mean nRF = {np.mean(rfs)}")