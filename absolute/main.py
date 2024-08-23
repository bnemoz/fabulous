from flask import Flask, request, jsonify
from flask_cors import CORS
from dataclasses import dataclass
from enum import Enum
import dnachisel as dc
from abutils import Sequence
from abutils.io import read_fasta
import subprocess as sp
import abstar
import string
import random
import re

app = Flask(__name__)
cors = CORS(app, origins=['*'], allow_headers=['Content-Type', 'Access-Control-Allow-Origin', 'Access-Control-Allow-Headers', 'Access-Control-Allow-Methods'])


####### Declare some key values and dataclasses for future use

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



####### Python functions to be used by the API 

def preprocessing(raw_input, species="human", debug=False):
    """Pre-processes the input to figure out input type (raw unformatted, single-FASTA or multi-FASTA) and the origin (nucleotides/DNA or amino-acids/proteins)
    Formats the input into a suitable format (a list of FabulousAb dataclass instances) for the downstream analysis"""
    
    input_seq_type = infer_input(raw_input)
    if debug:
        print(input_seq_type)
        
    clean_inputs = []
    
    if input_seq_type == "unformatted sequence":
        if raw_input == "":
            return None
        fasta_original_input = Sequence(raw_input, id=generate_random_label(8))
        residue_type = infer_residues(fasta_original_input.sequence, )
        if debug:
            print(residue_type)
        if residue_type == "protein":
            fasta_nt_input = reverse_translate(fasta_original_input, )
        elif residue_type == "DNA":
            fasta_nt_input = fasta_original_input
        ab = FabulousAb(name=fasta_original_input.id, raw_input=fasta_original_input.sequence, input_type=residue_type, formatted_input=fasta_nt_input.sequence, species=species)
        clean_inputs.append(ab)
        return clean_inputs
        
    elif input_seq_type == "single fasta":
        fasta_original_input = read_fasta(raw_input)[0]
        residue_type = infer_residues(fasta_original_input.sequence, )
        if debug:
            print(residue_type)
        if residue_type == "protein":
            fasta_nt_input = reverse_translate(fasta_original_input, )
        elif residue_type == "DNA":
            fasta_nt_input = fasta_original_input
        ab = FabulousAb(name=fasta_original_input.id, raw_input=fasta_original_input.sequence, input_type=residue_type, formatted_input=fasta_nt_input.sequence, species=species)
        clean_inputs.append(ab)
        return clean_inputs
    
    elif input_seq_type == "multi fasta" :
        _sequences = read_fasta(raw_input)
        for _seq in _sequences:
            residue_type = infer_residues(_seq.sequence)
            if residue_type == "DNA":
                ab = FabulousAb(name=_seq.id, raw_input=_seq.sequence, input_type=residue_type, formatted_input=_seq.sequence, species=species)
                clean_inputs.append(ab)
            elif residue_type == "protein":
                fasta_nt_input = reverse_translate(_seq, )
                ab = FabulousAb(name=_seq.id, raw_input=_seq.sequence, input_type=residue_type, formatted_input=fasta_nt_input.sequence, species=species)
                clean_inputs.append(ab)
        return clean_inputs


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
    if len(letters) < 5:
        return "DNA"
    elif len(letters) > 10:
        return "protein"


def reverse_translate(sequence, ):
    """Returns the reverse translated sequence (NT) of an input amino acid sequence."""
    try:
        dna = dc.reverse_translate(sequence.sequence, randomize_codons=True, table='Standard')
        ab_nt = Sequence(dna, id=sequence.id)
        return ab_nt
    except:
        return None

def cleaner(sequence, pure_DNA=False, ):
    """Clean a sequence by removing spaces, dashes, and newlines."""
    _seq = str(sequence)
    _seq = sequence.replace(" ", "")
    _seq = _seq.replace("-", "")
    _seq = _seq.replace("\n","")
    if pure_DNA:
        _seq = _seq.replace("U", "T")
    return _seq.upper()


def antibody_identification(fabulous_ab, debug=False, ):
    """Identify the antibody germline and CDRs using AbStar. Return the results as a Sequence object with annotations in a JSON dictionnary. Also encodes several key/values important for optimization and cloning"""
    _seq = Sequence(fabulous_ab.formatted_input, id=fabulous_ab.name)
    ab = abstar.run(_seq, germ_db=fabulous_ab.species, output_type="json", verbose=debug)
    ab["input_type"] = fabulous_ab.input_type
    ab["fabulous_input"] = fabulous_ab.raw_input
    ab['fr4_nt_mod3'] = longest_substring(ab['fr4_nt'])
    ab['germ_alignments_nt']['var']['html_mid'] = re.sub(' ', '.', ab['germ_alignments_nt']['var']['midline'])
    ab['germ_alignments_nt']['join']['html_mid'] = re.sub(' ', '.', ab['germ_alignments_nt']['join']['midline'])
    ab['leader_nt'] = leader_nt['IGH' if ab['chain'] == "heavy" else 'IGL' if ab['chain'] == 'lambda' else 'IGK'].upper()
    ab['leader_aa'] = leader_aa['IGH' if ab['chain'] == "heavy" else 'IGL' if ab['chain'] == 'lambda' else 'IGK'].upper()
    return ab


def optimize(ab, debug=False, ):
    # To-Do
    pass


def clone(ab, debug=False, ):
    # To-Do
    pass


def numbering(ab, debug=False, ):
    ab['numbering'] = {}
    ab['numbering']['kabat'] = anarci_wrap(ab, 'kabat', debug)
    ab['numbering']['IMGT'] = anarci_wrap(ab, 'IMGT', debug)
    ab['numbering']['chothia'] = anarci_wrap(ab, 'chothia', debug)
    ab['numbering']['martin'] = anarci_wrap(ab, 'Martin', debug)
    ab['numbering']['Aho'] = anarci_wrap(ab, 'Aho', debug)
    ab['numbering']['wolfguy'] = anarci_wrap(ab, 'Wolfguy', debug)
    return ab


def anarci_wrap(ab, numbering_scheme='IMGT', debug=False):
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


def longest_substring(string):
    longest = ""
    for i in range(len(string)):
        for j in range(i + 3, len(string) + 1, 3):
            substring = string[i:j]
            if len(substring) > len(longest):
                longest = substring
    return longest    





####### API routes 

@app.route('/', methods=['GET'])
def index():
    return "Welcome to the Antibody Identification API!"


@app.route('/id', methods=['GET'])
def id():
    sequence = request.args.get('sequence')
    sequence = sequence.replace("%3E", ">")
    
    try:
        species = request.args.get('species')
    except:
        species = "human"
    
    preprocessed = preprocessing(sequence, species=species) 

    results = []
    if preprocessed is not None:
        for _item in preprocessed:
            ab = antibody_identification(_item, debug=False)
            results.append(ab)

    return (dict(results[0]))


@app.route('/optimize', methods=['GET', 'POST'])
def optimize():
    # To-Do
    pass


@app.route('/clone', methods=['GET', 'POST'])
def clone():
    # To-Do
    pass

@app.route('/number', methods=['GET', 'POST'])
def number():
    #To-do
    pass


####### App caller 

if __name__ == '__main__':
    app.run(debug=True, port=8282)
