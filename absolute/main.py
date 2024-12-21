#############################################################################################################
#                                                                                                           #
#               Imports                                                                                     #
#                                                                                                           #
#############################################################################################################


from flask import Flask, request, jsonify
from flask_graphql import GraphQLView
from graphene import ObjectType, String, Schema, Field, List
from flask_cors import CORS

from dataclasses import dataclass
from enum import Enum
from collections import Counter

import string
import random
# import re

# from Bio.Seq import Seq
import dnachisel as dc
from abutils import Sequence
from abutils.io import read_fasta
# from abutils.io import read_airr, list_files
from abutils import alignment
# import antpack
from antpack import SingleChainAnnotator
# from antpack import PairedChainAnnotator
import abutils
import subprocess as sp
import abstar






#############################################################################################################
#                                                                                                           #
#               Defining App                                                                                #
#                                                                                                           #
#############################################################################################################
 

# Create GraphQL schema
schema = Schema(query=Query)

# Create Flask app
app = Flask(__name__)
cors = CORS(app, 
            origins=['*', 'http://localhost:3000'], 
            allow_headers=['Content-Type', 'Access-Control-Allow-Origin', 'Access-Control-Allow-Headers', 'Access-Control-Allow-Methods'],
            methods=['GET', 'POST', 'OPTIONS'],
            )

# Add /graphql endpoint
app.add_url_rule(
    '/graphql',
    view_func=GraphQLView.as_view(
        'graphql',
        schema=schema,
        graphiql=True  # Enables GraphiQL interface for testing
    )
)






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


# Define GraphQL types
class AntibodyType(ObjectType):
    sequence_id = String()
    raw_input = String()
    input_type = String()
    formatted_input = String()
    species = String()


# Create GraphQL Query class
class Query(ObjectType):
    antibody = Field(
        AntibodyType,
        sequence_id=String(required=True),
        sequence=String(required=True),
        species=String(default_value="human")
    )

    # Corrected resolver for 'antibody'
    def resolve_antibody(self, info, sequence_id, sequence, species):
        try:
            # Preprocess the input
            preprocessed = preprocessing(sequence_id, sequence, species=species)
            
            # Identify antibody
            result = antibody_identification(preprocessed)
            
            # Build the response object
            return AntibodyType(
                sequence_id=sequence_id,
                raw_input=preprocessed.raw_input,
                input_type=preprocessed.input_type.name if isinstance(preprocessed.input_type, Enum) else str(preprocessed.input_type),
                formatted_input=preprocessed.formatted_input,
                species=species
            )
        except Exception as e:
            # Return a GraphQL-compatible error
            raise Exception(f"Error processing antibody: {e}")





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



def antibody_identification(fabulous_ab, debug=False, ):
    """Identify the antibody germline and CDRs using AbStar. Return the results as a Sequence object with annotations in a JSON dictionnary. Also encodes several key/values important for optimization and cloning"""
    
    # Initial annotation with AbStar
    _seq = Sequence(fabulous_ab.formatted_input, id=fabulous_ab.name)
    ab = abstar.run(_seq, germline_database=fabulous_ab.species, verbose=debug)

    # Adding Fab'ulous specific annotations to the AbStar output
    ab["input_type"] = fabulous_ab.input_type
    ab["fabulous_input"] = fabulous_ab.raw_input

    return ab


def optimize(ab, species='homo sapiens', debug=False, ):
    return ab


def clone(ab, debug=False, ):
    # To-Do
    return ab


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








#############################################################################################################
#                                                                                                           #
#               API routes                                                                                  #
#                                                                                                           #
#############################################################################################################


@app.route('/', methods=['GET'])
def index():
    return "Welcome to the Antibody Identification API!"


@app.route('/id', methods=['GET', 'POST'])
def id():
    if request.method == 'POST':
        # For POST requests, fetch JSON data or form-encoded data
        data = request.get_json() or request.form
        sequence_id = data.get('sequence_id')
        sequence = data.get('sequence')
        species = data.get('species', 'human')  # Default to "human" if not provided
    else:
        # For GET requests, use query parameters
        sequence_id = request.args.get('sequence_id')
        sequence = request.args.get('sequence')
        species = request.args.get('species', 'human')

    # Process the sequence
    try:
        preprocessed = preprocessing(sequence_id, sequence, species=species)
    except Exception as e:
        return jsonify({"error": str(e)}), 400
    try:
        result = antibody_identification(preprocessed, debug=False)
    except Exception as e:
        return jsonify({"error": str(e)}), 400

    # Return the results in a consistent format
    if result:
        return jsonify(dict(result))  # Convert the dictionary to JSON for the response
    else:
        return jsonify({"error": "No results found"}), 400
    

@app.route('/ids', methods=['POST'])
def ids():
    data = request.get_json() or request.form

    preprocessed = []

    for obj in data:
        sequence_id = obj.get('sequence_id')
        sequence = obj.get('sequence')
        species = obj.get('species', 'human')
        ab = preprocessing(sequence_id, sequence, species=species)
        preprocessed.append(ab)

    results = []

    if preprocessed is not None:
        for _item in preprocessed:
            ab = antibody_identification(_item, debug=False)
            results.append(ab)

    return jsonify(results)


@app.route('/optimize', methods=['GET', 'POST'])
def optimize():
    ab = request.get_json() or request.form
    ab = optimize(ab)
    return ab


@app.route('/clone', methods=['GET', 'POST'])
def clone():
    ab = request.get_json() or request.form
    ab = clone(ab)
    return ab


@app.route('/number', methods=['POST'])
def number():
    ab = request.get_json() or request.form
    ab = numbering(ab)
    return ab


@app.route('/annotate', methods=['POST'])
def annotate():
    ab = request.get_json() or request.form
    ab = abnotator(ab)
    return ab


@app.route('/clonotype', methods=['POST'])
def clonotype():
    data = request.get_json() or request.form
    
    clusters = abutils.tools.cluster.cluster(data, )
    return clusters


@app.route('/phylogeny', methods=['POST'])
def phylogeny():
    data = request.get_json() or request.form
    
    return data





####### App caller 

if __name__ == '__main__':
    app.run(debug=True, port=8282)
