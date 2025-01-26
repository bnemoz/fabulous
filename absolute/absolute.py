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

import abstar
import abutils
from abutils import Sequence, alignment
from abutils.io import read_fasta
# from abutils.io import read_airr, list_files

# import antpack
from antpack import SingleChainAnnotator
# from antpack import PairedChainAnnotator
from antpack import HumanizationTool

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


# Define GraphQL schema
type_defs = """
    type Query {
        hello: String!
        antibody(sequence_id: String!, sequence: String!, species: String = "human"): Antibody!
    }

    type Antibody {
        sequence_id: String!
        raw_input: String!
        input_type: String!
        formatted_input: String!
        species: String!
    }
"""

# Create resolvers
query = QueryType()

@query.field("hello")
def resolve_hello(*_):
    return "Hello, GraphQL!"

@query.field("antibody")
def resolve_antibody(_, info, sequence_id, sequence, species="human"):
    try:
        preprocessed = preprocessing(sequence_id, sequence, species=species)
        result = antibody_identification(preprocessed)
        return {
            "sequence_id": sequence_id,
            "raw_input": preprocessed.raw_input,
            "input_type": str(preprocessed.input_type),
            "formatted_input": preprocessed.formatted_input,
            "species": species,
        }
    except Exception as e:
        raise Exception(f"Error processing antibody: {e}")

# Create executable schema
schema = make_executable_schema(type_defs, query)




#############################################################################################################
#                                                                                                           #
#               Python functions to be used by the API                                                      #
#                                                                                                           #
#############################################################################################################


def billing(user, token, app):
    with open(f'/home/serveradmin/Apps/fabulous/billing/billing_{user}.csv', 'a') as f:
        f.write(f"{user}\t{token}\t{app}\t{datetime.datetime.now()}\n")
    return None