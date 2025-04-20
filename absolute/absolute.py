#!/usr/bin/env python3


#############################################################################################################
#                                                                                                           #
#               Imports                                                                                     #
#                                                                                                           #
#############################################################################################################

from dataclasses import dataclass
from enum import Enum
from collections import Counter

import os
import string
import random
import datetime
# import re

import pandas as pd
import numpy as np

# # from Bio.Seq import Seq
import dnachisel as dc

from genbank import create_gb_from_ab

import abstar
import abutils
from abutils import Sequence
from abutils.tools import alignment
from abutils.io import read_fasta
# from abutils.io import read_airr, list_files

# import antpack
from antpack import SingleChainAnnotator
# from antpack import PairedChainAnnotator
# from antpack import HumanizationTool

import subprocess as sp


#############################################################################################################
#                                                                                                           #
#               Key values and dataclasses                                                                  #
#                                                                                                           #
#############################################################################################################


# Seamless cloning sequences:
oh_5 = {'IGH': 'gctgggttttccttgttgctattctcgagggtgtccagtgt',
        'IGK': 'atcctttttctagtagcaactgcaaccggtgtacac',
        'IGL': 'atcctttttctagtagcaactgcaaccggtgtacac'}
oh_3 = {'IGH': 'gctagcaccaagggcccatcggtcttcc',
        'IGK': 'cgtacggtggctgcaccatctgtcttcatc',
        'IGL': 'ggtcagcccaaggctgccccctcggtcactctgttcccgccctcgagtgaggagcttcaagccaacaaggcc'}

# Full length signal/leader sequences:
leader_nt = {'IGH': 'ATGGAACTGGGGCTCCGCTGGGTTTTCCTTGTTGCTATTCTCGAGGGTGTCCAGTGT',
             'IGK': 'ATGGGTTGGTCATGTATCATCCTTTTTCTAGTAGCAACTGCAACCGGTGTACAC',
             'IGL': 'ATGGGTTGGTCATGTATCATCCTTTTTCTAGTAGCAACTGCAACCGGT'}
leader_aa = {'IGH': 'MELGLRWVFLVAILEGVQC',
             'IGK': 'MGWSCIILFLVATATGVH',
             'IGL': 'MGWSCIILFLVATATG'}

# Isotype dictionary:
isotypes = {'IGHG1': 'IgG1',
            'IGHG2': 'IgG2',
            'IGHG3': 'IgG3',
            'IGHG4': 'IgG4',
            'IGHE': 'IgE',
            'IGHA1': 'IgA1',
            'IGHA2': 'IgA2',
            'IGHM': 'IgM',
            'IGHD': 'IgD',
            'IGK': 'kappa',
            'IGL': 'lambda'
            }


class InputType(Enum):
    DNA = 1
    PROTEIN = 2

class Species(Enum):
    HUMAN = 1
    MOUSE = 2
    MONKEY = 3
    HUMANIZED = 4

@dataclass
class FabulousAb:
    """Dataclass for keeping track of an antibody in the Fab'ulous environment."""
    name: str
    raw_input: str
    input_type: InputType
    formatted_input: str
    species: Species
    
    def __init__(self, name, raw_input, input_type, formatted_input, species):
        self.name = name
        self.raw_input = cleaner(raw_input)
        self.input_type = input_type
        self.formatted_input = cleaner(formatted_input, pure_DNA=True)
        self.species = species

    def __str__(self):
        return (self.name+'\n'+self.input_type)
    def __name__(self):
        return (self.name+'\n'+self.input_type)




#############################################################################################################
#                                                                                                           #
#               Python functions to be used by the API                                                      #
#                                                                                                           #
#############################################################################################################


def preprocessing(sequence_id, raw_input, species="human", debug=False):
    """Pre-processes the input to figure out input type (raw unformatted, single-FASTA or multi-FASTA) and the origin (nucleotides/DNA or amino-acids/proteins)
    Formats the input into a suitable format (a list of FabulousAb dataclass instances) for the downstream analysis"""
    
    # input_seq_type = infer_input(raw_input)
    # if debug:
    #     print(input_seq_type)
        
    # _seq = str(raw_input).replace(" ", "").replace("-", "").replace("\n","")

    residue_type = infer_residues(raw_input, )
    if debug:
        print(residue_type)
    
    if residue_type == 'protein':
        formatted_input = reverse_translate(raw_input, )
    elif residue_type == 'DNA':
        formatted_input = raw_input

    ab = FabulousAb(name=sequence_id, raw_input=raw_input, input_type=residue_type, formatted_input=formatted_input, species=species)
    return ab


def infer_input(sequence, ):
    """Infer the input format of a sequence: unformatted, single fasta, or multi fasta."""
    _input = str(sequence)
    carret_count = _input.count('>')
    if carret_count == 0:
        return 'unformatted sequence'
    elif carret_count == 1:
        return 'single fasta'
    elif carret_count >= 2:
        return 'multi fasta'


def generate_random_label(length=8):
    """Generate a random label of a specified length."""
    letters = string.ascii_letters + string.digits
    return ''.join(random.choice(letters) for _ in range(length))


def infer_residues(sequence, ):
    """Infer the rersidue type of a sequence: DNA or protein."""
    _input = str(sequence)
    letters = set([l.lower() for l in _input])
    if len(letters) < 6:
        return "DNA"
    elif len(letters) > 10:
        return "protein"


def reverse_translate(sequence, ):
    """Returns the reverse translated sequence (NT) of an input amino acid sequence."""
    try:
        dna = dc.reverse_translate(sequence, randomize_codons=True, table='Standard')
        return dna
    except:
        return None


def cleaner(sequence, pure_DNA=False, ):
    """Clean a sequence by removing spaces, dashes, and newlines."""
    if sequence is not None:
        _seq = str(sequence)
        _seq = sequence.replace(" ", "")
        _seq = _seq.replace("-", "")
        _seq = _seq.replace("\n","")
        if pure_DNA:
            _seq = _seq.replace("U", "T")
        return _seq.upper()
    

def enforce_modulo3(dna_seq: str) -> str:
    """
    Adjusts a DNA sequence to ensure its length is a multiple of 3 by 
    removing nucleotides from the end if necessary.

    Parameters:
        dna_seq (str): The input DNA sequence.

    Returns:
        str: The adjusted DNA sequence with length as a multiple of 3.
    """
    remainder = len(dna_seq) % 3
    if remainder:
        dna_seq = dna_seq[:-remainder]  # Remove the last 'remainder' nucleotides
    return dna_seq



def antibody_identification(fabulous_ab, debug=False, ):
    """Identify the antibody germline and CDRs using AbStar. Return the results as a Sequence object with annotations in a JSON dictionnary. Also encodes several key/values important for optimization and cloning"""
    
    errors = []

    # Initial annotation with AbStar
    _seq = Sequence(fabulous_ab.formatted_input, id=fabulous_ab.name)
    try:
        ab = abstar.run(_seq, germline_database=fabulous_ab.species, verbose=debug)
        if debug:
            print(datetime.datetime.now(), "AbStar run successful")
    except Exception as e:
        if debug:
            print(datetime.datetime.now(), "AbStar error: ", e)
        errors.append(str(e))
        return None, errors
    
    # Adding Fab'ulous specific annotations to the AbStar output
    ab["input_type"] = fabulous_ab.input_type
    ab["user_input"] = fabulous_ab.raw_input
    ab['v_gaps'] = ab['v_sequence_gapped'].count('.')
    ab['sequence_vdjc_aa_gapped'] = gapper(ab['sequence_vdjc_aa'], ab['sequence_vdjc_gapped'])
    ab['sequence_aa_gapped'] = gapper(ab['sequence_aa'], ab['sequence_gapped'])
    
    try:
        ab['chain'] = 'Heavy' if ab['locus'] == 'IGH' else 'Kappa' if ab['locus'] == 'IGK' else 'Lambda' if ab['locus'] == 'IGL' else None
    except:
        ab['chain'] = "Unknown"
        errors.append("Chain could not be determined")
    try:
        ab['isotype'] = isotypes[ab['c_call'].split('*')[0]] if ab['locus'] == 'IGH' else isotypes[ab['locus']]
    except:
        ab['isotype'] = "Unknwon"
        errors.append("Isotype could not be determined")

    # Calculating SHM 
    ab['SHM_v_nt'] = (1 - ab['v_identity']) * 100
    ab['SHM_v_aa'] = (1 - ab['v_identity_aa']) * 100
    if ab['j_identity'] is not None:
        ab['SHM_vj_nt'] = ((1 - ab['vj_identity']) + (1 - ab['j_identity'])) * 100
        ab['SHM_vj_aa'] = ((1 - ab['vj_identity_aa']) + (1 - ab['j_identity_aa'])) * 100
    else:
        ab['SHM_vj_nt'] = None
        ab['SHM_vj_aa'] = None
        errors.append("SHM VJ could not be calculated")

    # Adding leader information if available
    try:
        vdj_start = ab['sequence_oriented'].find(ab['fwr1'])
        upstream = ab['sequence_oriented'][:vdj_start]
        ab['leader'] = upstream
        if upstream != "":
            ab = assign5prime(ab)
    except:
        errors.append("5' leader sequence not found")

    # Adding trailing sequence annotations if available
    try:
        vdj_stop = ab['sequence_oriented'].find(ab['fwr4']) + len(ab['fwr4'])
        downstream = ab['sequence_oriented'][vdj_stop:]
        ab['trailer'] = downstream
        if downstream != "":
            ab = assign3prime(ab)
    except:
        errors.append("3' trailing sequence not found")

    return ab, errors



def optimize(ab, species='human', debug=False, ):
    sequence = enforce_modulo3(ab.sequence)

    species_dict = {'human':'h_sapiens', 'mouse':'m_musculus', }
    species = species_dict[species.lower()]

    optimize = dc.DnaOptimizationProblem(sequence=sequence, constraints=[dc.EnforceTranslation(), 
                                                        dc.EnforceGCContent(maxi=0.56), 
                                                        dc.EnforceGCContent(maxi=0.64, window=60), 
                                                        dc.UniquifyAllKmers(10), ], 
                                            logger=None,
                                            objectives=[dc.CodonOptimize(species=species)])
    optimize.resolve_constraints(final_check=True)

    try:
        optimize.optimize()
        
        ab['optimized_vdj'] = optimize.sequence
        ab['optimized_species'] = species
        ab['optimizations'] = optimize.sequence_edits_as_array().tolist()
        ab['optimizations_count'] = int(optimize.number_of_edits())
        ab['optimized_gc_content'] = (optimize.sequence.count('G') + optimize.sequence.count('C'))/len(optimize.sequence)*100

        ab['optimization_timestamp'] = str(datetime.datetime.now())
        errors = []

    except Exception as e:
        if debug:
            print("Optimization failed: error", e)
            errors = [f"Optimization failed: error {e}"]

    return ab, errors


def clone(ab, vector, debug=False, ):
    vector_type = vector.get('type')

    if ab['optimized_vdj'] is not None:
        input_seq = ab['optimized_vdj']
    else:
        input_seq = ab['sequence']
        errors = ["Warning: No optimized sequence available. Using the original sequence for cloning"]

    if not vector_type == "custom":
        if ab['locus'] == 'IGH':
            clonable = Sequence((oh_5['IGH'] + input_seq + oh_3['IGH']), id=ab.id)
            ab['clonable'] = clonable.sequence
        elif ab['locus'] == 'IGK':
            clonable = Sequence((oh_5['IGK'] + input_seq + oh_3['IGK']), id=ab.id)
            ab['clonable'] = clonable.sequence
        elif ab['locus'] == 'IGL':
            clonable = Sequence((oh_5['IGL'] + input_seq + oh_3['IGL']), id=ab.id)
            ab['clonable'] = clonable.sequence

    else:
        vector_seq = vector.get('sequence')
        vector_name = vector.get('sequence_id')
        vector_5oh = vector.get('5prime_overhang')
        vector_3oh = vector.get('3prime_overhang')

        clonable = Sequence((vector_5oh + input_seq + vector_3oh), id=str(ab.id+"_"+vector_name))
        plasmid = Sequence(vector_seq+clonable, id=vector_name)
        ab['clonable'] = clonable.sequence
        ab['plasmid'] = plasmid.sequence

    if debug:
        return ab, errors
    else:
        return ab, {}


def numbering(ab, scheme, algo='anarci', debug=False, ):
    if algo == 'anarci':
        ab['numbering'] = {}
        ab['numbering'][scheme] = anarci_wrap(ab, scheme, debug)
        return ab
    elif algo == 'antpack':
        ab['numbering'] = {}
        ab['numbering'][scheme] = antpack_wrap(ab, scheme, debug)
        return ab


def anarci_wrap(ab, numbering_scheme='IMGT', debug=False):
    scheme_dict = {'IMGT':'i',
                   'kabat':'k',
                   'chothia':'c',
                   'Martin':'m',
                   'Aho':'a',
                   'Wolfguy':'w',
                   }
    cmd = f"ANARCI -i {ab['sequence']} --scheme {scheme_dict[numbering_scheme]} --hmmerpath /home/serveradmin/antibody/.venv/bin"
    anarci_cmd = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True, universal_newlines=True)
    stdout, stderr = anarci_cmd.communicate()
    if debug:
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


def antpack_wrap(ab, numbering_scheme='IMGT', debug=False):
    scheme_dict = {'IMGT':'imgt',
                   'kabat':'kabat',
                   'Martin':'martin',
                   'Aho':'aho',
                   }
    
    if numbering_scheme not in scheme_dict.keys():
        return "Sorry :-/ Numbering scheme not supported"
    
    sequence = ab['sequence_aa']
    chain = 'H' if ab['locus'] == 'heavy' else 'K' if ab['locus'] == 'kappa' else 'L'
    sc_annotator = SingleChainAnnotator([chain, ], scheme = scheme_dict[numbering_scheme])
    numbered, percent_identity, chain_type, err_message = sc_annotator.analyze_seq(sequence)
    try:
        output = [(a,z) for a,z in zip(numbered, sequence)]
    except:
        output = err_message

    return output

def longest_substring(string):
    longest = ""
    for i in range(len(string)):
        for j in range(i + 3, len(string) + 1, 3):
            substring = string[i:j]
            if len(substring) > len(longest):
                longest = substring
    return longest    


def assign5prime(ab):
    leaders = read_fasta('./refs/L.fasta')
    alns = alignment.semiglobal_alignment(query=ab['leader'], targets=leaders)
    maxi = max(Counter([a.score for a in alns]))
    aln = [a for a in alns if a.score == maxi][0]
    ab['signal_peptide_id'] = aln.target_id
    ab['signal_peptide_score'] = maxi
    ab['signal_peptide_sequence'] = aln.query.sequence[aln.query_begin:]
    ab['signal_peptide_sequence_aa'] = dc.translate(aln.query.sequence[aln.query_begin:])
    ab['signal_peptide_start'] = aln.query_begin

    ab['5-UTR'] = aln.query.sequence[:aln.query_begin]
    
    part1 = read_fasta('./refs/L1.fasta')
    alns = alignment.semiglobal_alignment(query=ab['signal_peptide_sequence'], targets=part1)
    maxi = max(Counter([a.score for a in alns]))
    aln = [a for a in alns if a.score == maxi][0]
    ab['L_part_1'] = aln.query.sequence[:aln.target_end+1]
    ab['L_part_2'] = aln.query.sequence[aln.target_end+1:]
    
    return ab


def assign3prime(ab):
    ab = assign_ch1(ab)
    if ab['locus'] != 'IGH':
        return ab
    else:
        ab = assign_ch2(ab)
        ab = assign_ch3(ab)
        if any([(ab['c_call'].startswith('IgE')), (ab['c_call'].startswith('IgM'))]):
            ab = assign_ch4(ab)
        ab = assign_h(ab)
        return ab


def assign_ch1(ab):
    query = ab['trailer']
    if ab['locus'] == 'IGH':
        ch1 = read_fasta('./refs/CH1.fasta')
    else:
        ch1 = read_fasta('./refs/CL1.fasta')
    alns = alignment.semiglobal_alignment(query=query, targets=ch1)
    maxi = max(Counter([a.score for a in alns]))
    aln = [a for a in alns if a.score == maxi][0]
    ab['ch1'] = aln.target_id
    ab['ch1_sequence'] = aln.aligned_query[:aln.target_end+1].replace('-','')
    ab['ch1_sequence_aa'] = dc.translate(ab.sequence[-1:]+ab['ch1_sequence'])
    return ab


def assign_ch2(ab):
    query = ab['trailer'].split(ab['ch1_sequence'])[-1]
    if len(query) < 260:
        return ab
    ch2 = read_fasta('./refs/CH2.fasta')
    alns = alignment.semiglobal_alignment(query=query, targets=ch2)
    maxi = max(Counter([a.score for a in alns]))
    aln = [a for a in alns if a.score == maxi][0]
    ab['ch2'] = aln.target_id
    ab['ch2_sequence'] = aln.aligned_query[:aln.target_end+1].replace('-','')
    ab['ch2_sequence_aa'] = dc.translate(ab['ch1_sequence'][-1:]+ab['ch2_sequence'])
    return ab


def assign_ch3(ab):
    try:
        query = ab['trailer'].split(ab['ch2_sequence'])[1]
        if len(query) < 260:
            return ab
        ch3 = read_fasta('./refs/CH3.fasta')
        alns = alignment.semiglobal_alignment(query=query, targets=ch3)
        maxi = max(Counter([a.score for a in alns]))
        aln = [a for a in alns if a.score == maxi][0]
        ab['ch3'] = aln.target_id
        ab['ch3_sequence'] = aln.aligned_query[:aln.target_end+1].replace('-','')
        ab['ch3_sequence_aa'] = dc.translate(ab['ch2_sequence'][-2:]+ab['ch3_sequence'])
        return ab
    except:
        return ab


def assign_ch4(ab):
    query = ab['trailer'].split(ab['ch3_sequence'])[1]
    ch4 = read_fasta('./refs/CH4.fasta')
    alns = alignment.semiglobal_alignment(query=query, targets=ch4)
    maxi = max(Counter([a.score for a in alns]))
    aln = [a for a in alns if a.score == maxi][0]
    ab['ch4'] = aln.target_id
    ab['ch4_sequence'] = aln.aligned_query[:aln.target_end+1].replace('-','')
    ab['ch4_sequence_aa'] = dc.translate(ab['ch4_sequence'])
    return ab


def assign_h(ab):
    query = ab['trailer']
    h = read_fasta('./refs/H.fasta')
    alns = alignment.semiglobal_alignment(query=query, targets=h)
    maxi = max(Counter([a.score for a in alns]))
    aln = [a for a in alns if a.score == maxi][0]
    ab['hinge'] = aln.target_id
    ab['hinge_sequence'] = aln.aligned_query[aln.query_begin:aln.query_end+1].replace('-','')
    ab['hinge_sequence_aa'] = dc.translate(aln.aligned_query[aln.query_begin-1]+ab['hinge_sequence'],)
    return ab


def abnotator(ab, debug=False, ):
    try:
        leader, trailer = ab['sequence_input'].split(ab.sequence)
        ab['leader'] = leader
        ab['trailer'] = trailer
    except:
        ab['leader'] = None
        ab['trailer'] = None
        return ab

    ab = assign5prime(ab)
    ab = assign3prime(ab)
  
    return ab

def gapper(seq_aa_to_gap, seq_nt_gapped):
    """Returns the gapped AA sequence matching the codon positions in a gapped NT sequence."""
    seq_aa = str(seq_aa_to_gap)
    seq_nt = str(seq_nt_gapped)

    nt_clean = seq_nt.replace('-', '').replace('.', '')
    translated = str(dc.translate(nt_clean))
    assert translated == seq_aa, f"Mismatch in AA translation:\nExpected: {seq_aa}\nGot: {translated}"
    assert len(nt_clean) % 3 == 0, "Ungapped NT sequence length should be a multiple of 3"

    gapped_aa = ""
    i = 0  # index in ungapped AA sequence

    for j in range(0, len(seq_nt), 3):
        codon = seq_nt[j:j+3]
        if len(codon) < 3:
            break
        if '.' in codon or '-' in codon:
            gapped_aa += '.'  # gap placeholder
        else:
            gapped_aa += seq_aa[i]
            i += 1

    return gapped_aa


def make_gb_file(ab, debug=False, ):
    """Creates a GenBank file from the antibody sequence and annotations."""

    to_file = f'/tmp/{generate_random_label(16)}.gb'

    gb = create_gb_from_ab(ab, to_file=to_file)

    return gb


def get_clusters(data):
    return



def get_phylogeny(data):
    return



def build_tree(data):
    return