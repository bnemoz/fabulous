#############################################################################################################
#                                                                                                           #
#               Imports                                                                                     #
#                                                                                                           #
#############################################################################################################


from flask import Flask, request, jsonify
# from flask_graphql import GraphQLView
# from graphene import ObjectType, String, Schema, Field, List
from flask_cors import CORS
from ariadne import QueryType, make_executable_schema, graphql_sync


from dataclasses import dataclass
from enum import Enum
from collections import Counter


import string
import random
import datetime
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
from antpack import HumanizationTool
import abutils
import subprocess as sp
import abstar






#############################################################################################################
#                                                                                                           #
#               Defining App                                                                                #
#                                                                                                           #
#############################################################################################################
 

# Create Flask app
app = Flask(__name__)
cors = CORS(app, 
            origins=['*', 'http://localhost:3000'], 
            allow_headers=['Content-Type', 'Access-Control-Allow-Origin', 'Access-Control-Allow-Headers', 'Access-Control-Allow-Methods'],
            methods=['GET', 'POST', 'OPTIONS'],
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
    ab["fabulous_input"] = fabulous_ab.raw_input

    try:
        ab['chain'] = 'Heavy' if ab['locus'] == 'IGH' else 'Kappa' if ab['locus'] == 'IGK' else 'Lambda' if ab['locus'] == 'IGL' else None
    except:
        ab['chain'] = "Unknown"
        errors.append("Chain could not be determined")
    try:
        ab['isotype'] = isotypes[ab['c_call']] if ab['locus'] == 'IGH' else isotypes[ab['locus']]
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

    return ab, errors



def optimize(ab, species='h_sapiens', debug=False, ):
    sequence = ab.sequence

    optimize = dc.DnaOptimizationProblem(sequence=sequence, constraints=[dc.EnforceTranslation(), 
                                                        dc.EnforceGCContent(maxi=0.56), 
                                                        dc.EnforceGCContent(maxi=0.64, window=60), 
                                                        dc.UniquifyAllKmers(10), ], 
                                            logger=None,
                                            objectives=[dc.CodonOptimize(species=species)])
    optimize.resolve_constraints(final_check=True)
    optimize.optimize()

    ab['optimized_vdj'] = optimize.sequence
    ab['optimized_species'] = species
    ab['optimizations'] = optimize.sequence_edits_as_array()
    ab['optimizations_count'] = optimize.number_of_edits()
    ab['optimized_gc_content'] = (optimize.sequence.count('G') + optimize.sequence.count('C'))/len(optimize.sequence)*100

    ab['optimization_timestamp'] = str(datetime.datetime.now())

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


def single_humanize(ab, temp:float = 1.25, debug=False, ):

    h_tool = HumanizationTool()
    original = ab['sequence']

    score, mutations, humanized = h_tool.suggest_mutations(original, 
                                                     excluded_positions = [],
                                                     s_thresh = knob,
                                                    )
    # aln = alignment.semiglobal_alignment(original, humanized)
    percent_change = round(len(mutations)/len(humanized) * 100, 2)

    ab['humanized'] = humanized
    ab['humanization_score'] = score
    ab['humanization_mutations'] = mutations
    ab['humanization_percent_change'] = percent_change

    return ab


def multi_humanize(sequence, oracles, iterations, seq_per_it, final_output, mutables_fwr, mutables_cdr, debug=False, ):

    return sequence



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
        userid = data.get('userid')
        authtoken = data.get('authtoken')
        debug = data.get('debug', False)
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
        result, errors = antibody_identification(preprocessed, debug=debug)
        if errors and debug:
            for err in errors:
                print(f"Ecountered this error:{err}")
    except Exception as e:
        return jsonify({"error": str(e)}), 400

    # Return the results in a consistent format
    if result:
        billing(user=userid, token=authtoken, app='id')
        return jsonify({"results":dict(result)})
    else:
        return jsonify({"error": "No results found (Fab'ulous Id App)"}), 400
    

# @app.route('/ids', methods=['POST'])
# def ids():
#     data = request.get_json() or request.form

#     preprocessed = []

#     for n in data:
#         obj = data[n]
#         sequence_id = obj.get('sequence_id')
#         sequence = obj.get('sequence')
#         species = obj.get('species', 'human')
#         userid = data.get('userid')
#         authtoken = data.get('authtoken')
#         debug = data.get('debug', False)

#         ab = preprocessing(sequence_id, sequence, species=species)
#         billing(user=userid, token=authtoken, app='ids')
#         preprocessed.append(ab)

#     results = {}

#     if preprocessed is not None:
#         for n, _item in enumerate(preprocessed):
#             ab, errors = antibody_identification(_item, debug=debug)
#             if errors and debug:
#                 for err in errors:
#                     print(f"Ecountered this error:{err}")
#             results[n] = ab

#     if results != {}:            
#         return jsonify({"results":dict(results)})
#     else:
#         return jsonify({"error": "No results found (Fab'ulous Ids App)"}), 400


@app.route('/ids', methods=['POST'])
def ids():
    data = request.get_json()
    if not data:
        return jsonify({"error": "Invalid or missing JSON payload"}), 400

    userid = data.get('userid')
    authtoken = data.get('authtoken')
    debug = data.get('debug', False)
    sequences = data.get('sequences')

    if not sequences or not isinstance(sequences, dict):
        return jsonify({"error": "Invalid or missing 'sequences' field"}), 400

    results = {}
    errors = []
    for seq_id, seq_data in sequences.items():
        sequence_id = seq_data.get('sequence_id')
        sequence = seq_data.get('sequence')
        species = seq_data.get('species', 'human')  # Default to 'human'

        # Process the sequence
        try:
            preprocessed = preprocessing(sequence_id, sequence, species=species)
        except Exception as e:
            errors.append({"sequence_id": sequence_id, "error": str(e)})
            continue

        try:
            result, seq_errors = antibody_identification(preprocessed, debug=debug)
            if seq_errors and debug:
                errors.extend([{"sequence_id": sequence_id, "error": err} for err in seq_errors])

            # Extract the dictionary from result.annotations
            if hasattr(result, "annotations") and isinstance(result.annotations, dict):
                results[seq_id] = result.annotations
                billing(user=userid, token=authtoken, app='ids')
            else:
                errors.append({"sequence_id": sequence_id, "error": "Result annotations are missing or invalid"})
        except Exception as e:
            errors.append({"sequence_id": sequence_id, "error": str(e)})
        

    response = {"results": results}
    if errors and debug:
        response["errors"] = errors

    return jsonify(response), 200 if results else 400




@app.route('/optimize', methods=['POST'])
def optimize():
    raw = request.get_json() or request.form
    header, data = raw
    ab, species = data
    userid = header.get('userid')
    authtoken = header.get('authtoken')
    ab = optimize(ab, species=species, )
    billing(user=userid, token=authtoken, app='optimize')
    return ab


@app.route('/clone', methods=['POST'])
def clone():
    raw = request.get_json() or request.form
    header, data = raw
    ab, vector = data
    userid = header.get('userid')
    authtoken = header.get('authtoken')
    ab = clone(ab, vector)
    billing(user=userid, token=authtoken, app='clone')
    return ab


@app.route('/number', methods=['POST'])
def number():
    raw = request.get_json() or request.form
    header, data = raw
    ab, numbering_scheme = data
    userid = header.get('userid')
    authtoken = header.get('authtoken')
    ab = numbering(ab, numbering_scheme)
    billing(user=userid, token=authtoken, app='numbering')
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


# @app.route('/humanize', methods=['POST'])
# def _humanize():
#     raw = request.get_json() or request.form
#     header, data = raw
#     ab = data.get('antibody')
#     knob = data.get('knob', 1.25)
#     userid = header.get('userid')
#     authtoken = header.get('authtoken')
#     ab = humanize(ab, knob=knob)
#     billing(user=userid, token=authtoken, app='humanize')
#     return ab


@app.route('/humanize', methods=['POST'])
def humanize():
    # Parse the input payload
    payload = request.get_json() or request.form
    if not payload:
        return jsonify({"error": "Invalid or missing payload"}), 400

    userid = payload.get('userid')
    authtoken = payload.get('authtoken')
    model = payload.get('model')
    debug = payload.get('debug', False)

    if model == "single":
        # Extract relevant data for single model
        temp = payload.get('temp')
        sequence_data = payload.get('sequence')
        
        if not sequence_data or not isinstance(sequence_data, dict):
            return jsonify({"error": "Invalid or missing 'sequence' for single model"}), 400

        sequence_id = sequence_data.get('sequence_id')
        sequence = sequence_data.get('sequence')
        species = sequence_data.get('species', 'Mouse')  # Default species

        ab = preprocessing(sequence_id, sequence, species=species)

        # Perform single humanization
        try:
            ab = single_humanize(ab=ab, temp=temp, debug=debug)
            billing(user=userid, token=authtoken, app='humanize_single')
        except Exception as e:
            return jsonify({"error": str(e)}), 500

    elif model == "multi":
        # Extract relevant data for multi model
        oracles = payload.get('oracles')
        iterations = payload.get('iterations')
        seq_per_it = payload.get('seqs_per_it')
        final_output = payload.get('final_output')
        mutables_fwr = payload.get('mutables_fwr')
        mutables_cdr = payload.get('mutables_cdr')
        initiation_data = payload.get('initiation')

        if not initiation_data or not isinstance(initiation_data, dict):
            return jsonify({"error": "Invalid or missing 'initiation' for multi model"}), 400

        sequence_id = initiation_data.get('sequence_id')
        sequence = initiation_data.get('sequence')
        species = initiation_data.get('species', 'Chimpanzee')  # Default species

        ab = preprocessing(sequence_id, sequence, species=species)

        # Perform multi humanization
        try:
            ab = multi_humanize(
                sequence=ab,
                oracles=oracles,
                iterations=iterations,
                seq_per_it=seq_per_it,
                final_output=final_output,
                mutables_fwr=mutables_fwr,
                mutables_cdr=mutables_cdr,
                debug=debug
            )
            billing(user=userid, token=authtoken, app='humanize_multi')
        except Exception as e:
            return jsonify({"error": str(e)}), 500

    else:
        return jsonify({"error": "Invalid model type. Must be 'single' or 'multi'"}), 400

    # Return the results
    return jsonify({"result": ab}), 200



@app.route("/graphql", methods=["POST"])
def graphql_server():
    data = request.get_json()
    success, result = graphql_sync(schema, data, context_value=request, debug=True)
    status_code = 200 if success else 400
    return jsonify(result), status_code


####### App caller 

if __name__ == '__main__':
    app.run(debug=True, port=8282)
