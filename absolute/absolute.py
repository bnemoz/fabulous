#!/usr/bin/env python3

import os
import sys
import string
import random
import multiprocessing
import itertools
import fastcluster as fc
from scipy.cluster.hierarchy import fcluster
from collections import Counter
import subprocess as sp
import re
# from tqdm import tqdm
import abstar
import abutils
from abutils import Sequence
from abutils.utils.utilities import nested_dict_lookup
from mnemonic import Mnemonic
import dnachisel
import numpy as np
import pandas as pd
from Bio.Seq import Seq



def infer_input(input, ):
    _input = str(input)
    carret_count = _input.count('>')
    if carret_count == 0:
        return 'unformatted sequence'
    elif carret_count == 1:
        return 'single fasta'
    elif carret_count >= 2:
        return 'multi fasta'


def cleaner(sequence, pure_DNA=False, ):
    _seq = sequence.replace(" ", "")
    _seq = _seq.replace("-", "")
    _seq = _seq.replace("\n","")
    if pure_DNA:
        _seq = _seq.replace("U", "T")
    return _seq.upper()


def generate_random_label(length=8):
    """Generate a random label of a specified length."""
    letters = string.ascii_letters + string.digits
    return ''.join(random.choice(letters) for _ in range(length))


def infer_species(sequence, ):
    _species = 'human'
    return _species


def number(ab):
    ab['numbering'] = {}
    ab['numbering']['kabat'] = anarci_wrap(ab, 'kabat')
    ab['numbering']['IMGT'] = anarci_wrap(ab, 'IMGT')
    ab['numbering']['chothia'] = anarci_wrap(ab, 'chothia')
    ab['numbering']['martin'] = anarci_wrap(ab, 'Martin')
    ab['numbering']['Aho'] = anarci_wrap(ab, 'Aho')
    ab['numbering']['wolfguy'] = anarci_wrap(ab, 'Wolfguy')
    return ab


def anarci_wrap(ab, numbering_scheme='IMGT', raw=False):
    scheme_dict = {'IMGT':'i',
                   'kabat':'k',
                   'chothia':'c',
                   'Martin':'m',
                   'Aho':'a',
                   'Wolfguy':'w',
                   }
    cmd = f"ANARCI -i {ab['vdj_aa']} --scheme {scheme_dict[numbering_scheme]} --hmmerpath /home/serveradmin/antibody/.venv/bin"
    anarci_cmd = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True, universal_newlines=True)
    stdout, stderr = anarci_cmd.communicate()
    if raw:
        return stdout
    else:
        raw = stdout.splitlines()[7:-1]
        numbering = {}
        try:
            for i, e in enumerate(raw):
                elems = e.split(' ')
                numbering[i] = [elems[1], elems[-1]]
        except TypeError:
            numbering = "Sorry :-/ Numbering scheme cannot be assessed"
        return numbering


def longest_substring(string):
    longest = ""
    for i in range(len(string)):
        for j in range(i + 3, len(string) + 1, 3):
            substring = string[i:j]
            if len(substring) > len(longest):
                longest = substring
    return longest


def check_isotype(ab):
    """Based on the isotype aligment score from AbStar, determines if isotype is reliable (returns True) or unreliable (return False)"""
    if len(ab['isotype_alignment']['midline']) >= 8:
        return True
    else:
        return False
        

def clonify_python(
    jsons,
    distance_cutoff=0.32,
    shared_mutation_bonus=0.65,
    length_penalty_multiplier=2,
    preclustering=False,
    preclustering_threshold=0.65,
    preclustering_field="cdr3_nt",
    lineage_field="lineage",
    lineage_size_field="lineage_size",
    return_assignment_dict=False,
):    
    vgene_key = "v_gene.gene"
    jgene_key = "j_gene.gene"
    cdr3_key = "cdr3_aa"
    muts_key = "var_muts_nt.muts"

    # group sequences by V/J genes
    vj_group_dict = {}
    for ab in jsons:
        if ab['vdj_nt'] is None:
            continue
        # build new Sequence objects using just the data we need
        s = Sequence(ab['vdj_nt'], id=ab['seq_id'])
        s["v_gene"] = nested_dict_lookup(ab, vgene_key.split("."))
        s["j_gene"] = nested_dict_lookup(ab, jgene_key.split("."))
        s["cdr3"] = nested_dict_lookup(ab, cdr3_key.split("."))
        
        muts = nested_dict_lookup(ab, muts_key.split("."), [])
        s["mutations"] = [f"{m['position']}:{m['was']}>{m['is']}" for m in muts]

        required_fields = ["v_gene", "j_gene", "cdr3", "mutations"]
        if preclustering:
            s["preclustering"] = nested_dict_lookup(ab, preclustering_field.split("."))
            required_fields.append("preclustering")
        if any([s[v] is None for v in required_fields]):
            continue
            
        # group sequences by VJ gene use
        vj = f"{s['v_gene']}__{s['j_gene']}"
        if vj not in vj_group_dict:
            vj_group_dict[vj] = []
        vj_group_dict[vj].append(s)
        
    # assign lineages
    assignment_dict = {}
    mnemo = Mnemonic("english")
    for vj_group in vj_group_dict.values():
        # preclustering
        if preclustering:
            seq_dict = {s.id: s for s in vj_group}
            cluster_seqs = [Sequence(s[preclustering_field], id=s.id) for s in vj_group]
            clusters = cluster(cluster_seqs, threshold=preclustering_threshold)
            groups = [[seq_dict[i] for i in c.seq_ids] for c in clusters]
        else:
            groups = [
                vj_group,
            ]
        for group in groups:
            if len(group) == 1:
                seq = group[0]
                assignment_dict[seq.id] = "_".join(
                    mnemo.generate(strength=128).split()[:6]
                )
                continue
            # build a distance matrix
            dist_matrix = []
            async_results = []
            with multiprocessing.Pool() as pool:
                for s1, s2 in itertools.combinations(group, 2):
                    async_results.append(pool.apply_async(_clonify_distance, args=(s1, s2, shared_mutation_bonus, length_penalty_multiplier)))
                for a in async_results:
                    dist_matrix.append(a.get())

            # cluster
            linkage_matrix = fc.linkage(
                dist_matrix, method="average", preserve_input=False
            )
            cluster_list = fcluster(
                linkage_matrix, distance_cutoff, criterion="distance"
            )
            
            # rename clusters
            cluster_ids = list(set(cluster_list))
            cluster_names = {
                c: "_".join(mnemo.generate(strength=128).split()[:6])
                for c in cluster_ids
            }
            renamed_clusters = [cluster_names[c] for c in cluster_list]
            
            # assign sequences
            for seq, name in zip(vj_group, renamed_clusters):
                assignment_dict[seq.id] = name
            lineage_size_dict = Counter(assignment_dict.values())
            
    # return assignments
    if return_assignment_dict:
        return assignment_dict
        
    for ab in jsons:
        try:
            ab[lineage_field] = assignment_dict[ab['seq_id']]
            ab[lineage_size_field] = lineage_size_dict[ab[lineage_field]] if lineage_size_dict[ab[lineage_field]] != 0 else 1
        except KeyError:
            ab[lineage_field] = None
            ab[lineage_size_field] = np.nan
        except UnboundLocalError:
            ab[lineage_size_field] = np.nan
    return jsons


def _clonify_distance(s1, s2, shared_mutation_bonus, length_penalty_multiplier):
    if len(s1["cdr3"]) == len(s2["cdr3"]):
        dist = sum([i != j for i, j in zip(s1["cdr3"], s2["cdr3"])])
    else:
        dist = distance(s1["cdr3"], s2["cdr3"])
    length_penalty = abs(len(s1["cdr3"]) - len(s2["cdr3"])) * length_penalty_multiplier
    length = min(len(s1["cdr3"]), len(s2["cdr3"]))
    shared_mutations = list(set(s1["mutations"]) & set(s2["mutations"]))
    mutation_bonus = len(shared_mutations) * shared_mutation_bonus
    score = (dist + length_penalty - mutation_bonus) / length
    return max(score, 0.001)  # distance values can't be negative


def optimize(sequence, ):
    # To be implemented
    return sequence



