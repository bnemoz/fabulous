from flask import Flask, request, jsonify
from flask_cors import CORS
import dnachisel as dc
from abutils import Sequence
from abutils.io import read_fasta
import abstar
import string
import random

app = Flask(__name__)
cors = CORS(app, origins=['*'], allow_headers=['Content-Type', 'Access-Control-Allow-Origin', 'Access-Control-Allow-Headers', 'Access-Control-Allow-Methods'])


####### Declare some key values for future use

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



####### Python functions to be used by the API 

def infer_input(input, ):
    """Infer the input format of a sequence: unformatted, single fasta, or multi fasta."""
    _input = str(input)
    carret_count = _input.count('>')
    if carret_count == 0:
        return 'unformatted sequence'
    elif carret_count == 1:
        return 'single fasta'
    elif carret_count >= 2:
        return 'multi fasta'


def infer_residues(input, ):
    """Infer the rersidue type of a sequence: DNA or protein."""
    _input = str(input)
    letters = set([l.lower for l in _input])
    if len(letters) < 5:
        return "DNA"
    elif len(letters) > 10:
        return "protein"
    

def cleaner(sequence, pure_DNA=False, ):
    """Clean a sequence by removing spaces, dashes, and newlines."""
    _seq = str(sequence)
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


def antibody_identification(sequence, species):
    """Main function to provide antibody identification. Returns a JSON object with nested dictionaries."""
    ab_json = abstar.run(sequence, germ_db=species, verbose=False)
    return ab_json


def reverse_translate(sequence, species):
    """Returns the reverse translated sequence (NT) of an input amino acid sequence."""
    ab_aa = Sequence(sequence, id=generate_random_label())
    dna = dc.reverse_translate(ab_aa.sequence, randomize_codons=True, table='Standard')
    ab_nt = Sequence(dna, id=ab_aa.id)
    return ab_nt




####### API routes 

@app.route('/', methods=['GET'])
def index():
    return "Welcome to the Antibody Identification API!"


@app.route('/id', methods=['GET'])
def id():
    sequence = request.args.get('sequence')
    seq_type = infer_input(sequence)
    residues = infer_residues(sequence)

    if seq_type == 'single fasta':
        sequence = read_fasta(sequence)
    elif seq_type == 'unformatted sequence':
        sequence = cleaner(sequence)
        sequence = Sequence(sequence, id=generate_random_label())
        seq_type = 'single fasta'
    elif seq_type == 'multi fasta':
        sequence = read_fasta(sequence)

    if sequence in ['', None]:
        ab_json = {}
        return ab_json

    if seq_type == 'single fasta':
        if residues == 'protein':
            ab_nt = reverse_translate(sequence, species='human')
            ab_json = antibody_identification(ab_nt, species='human')
        elif residues == 'DNA':
            ab_json = antibody_identification(sequence, species='human')
        return dict(ab_json)
    
    elif seq_type == 'multi fasta':
        antibodies = {}
        for ab in sequence:
            ab_json = antibody_identification(ab, species='human')
            antibodies[ab.id] = dict(ab_json)
        return antibodies


@app.route('/clone', methods=['GET', 'POST'])
def clone():
    sequence = request.args.get('sequence')
    clone = {}
    clone['original'] = sequence
    # clone['optimized'] = 
    return clone



####### App caller 

if __name__ == '__main__':
    app.run(debug=True, port=8282)
