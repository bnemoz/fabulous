#!/usr/bin/env python3


#############################################################################################################
#                                                                                                           #
#               Imports                                                                                     #
#                                                                                                           #
#############################################################################################################


from dataclasses import dataclass
# from enum import Enum
# from collections import Counter

# import os
# import string
# import random
# import datetime
# import re

# import pandas as pd
# import numpy as np

# # from Bio.Seq import Seq
# import dnachisel as dc

# import abstar
# import abutils
# from abutils import Sequence, alignment
# from abutils.io import read_fasta
# from abutils.io import read_airr, list_files

# import antpack
# from antpack import SingleChainAnnotator
# from antpack import PairedChainAnnotator
from antpack import HumanizationTool

# import subprocess as sp


#############################################################################################################
#                                                                                                           #
#               Python functions to be used by the API                                                      #
#                                                                                                           #
#############################################################################################################




def single_humanize(ab, temp:float = 1.25, debug=False, ):

    h_tool = HumanizationTool()
    original = ab['sequence']

    try:
        score, mutations, humanized = h_tool.suggest_mutations(original, 
                                                        excluded_positions = [],
                                                        s_thresh = float(temp),
                                                        )
        # aln = alignment.semiglobal_alignment(original, humanized)
        percent_change = round(len(mutations)/len(humanized) * 100, 2)

        ab['humanized'] = humanized
        ab['humanization_score'] = score
        ab['humanization_mutations'] = mutations
        ab['humanization_percent_change'] = percent_change

        if debug:
            print(f"Score: {score}")
            print(f"Percent change: {percent_change}")
    except Exception as e:
        ab['humanized'] = None
        ab['humanization_score'] = None
        ab['humanization_mutations'] = None
        ab['humanization_percent_change'] = None
        if debug:
            print(e)
            
    return ab



def multi_humanize(sequence, oracles, iterations, seq_per_it, final_output, mutables_fwr, mutables_cdr, debug=False, ):

    return sequence