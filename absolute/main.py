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
# from collections import Counter

# import string
# import random
# import datetime
# import re

# import abutils

from absolute import preprocessing, antibody_identification, optimize, clone, numbering, abnotator, get_clusters, get_phylogeny, build_tree
from humanize import single_humanize, multi_humanize
from billing import billing, get_bill, authenticate



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



#############################################################################################################
#                                                                                                           #
#               API routes                                                                                  #
#                                                                                                           #
#############################################################################################################


@app.route('/', methods=['GET'])
def index():
    return "Welcome to the Fab'ulous Antibody API!"



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
        billing(user=userid, token=authtoken, app='Identification (single)')
        if debug:
            return jsonify({"result": dict(result), "errors": errors})
        return jsonify({"results":dict(result)})
    else:
        return jsonify({"Single Identification App error": "No results found (Fab'ulous Id App)"}), 400
    


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

    if not authenticate(userid, authtoken):
        return jsonify({"Authentification App error": "Invalid or missing authentication"}), 400

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
                billing(user=userid, token=authtoken, app='Identification (Batch)')
            else:
                errors.append({"sequence_id": sequence_id, "error": "Result annotations are missing or invalid"})
        except Exception as e:
            errors.append({"sequence_id": sequence_id, "Batch Identification App error": str(e)})
        

    response = {"results": results}
    if errors and debug:
        response["errors"] = errors

    return jsonify(response), 200 if results else 400



@app.route('/optimize', methods=['POST'])
def optimizing():
    data = request.get_json()
    if not data:
        return jsonify({"error": "Invalid or missing JSON payload"}), 400

    userid = data.get('userid')
    authtoken = data.get('authtoken')
    debug = data.get('debug', False)
    sequence_id = data.get('sequence_id')
    sequence = data.get('sequence')
    species = data.get('species', 'human')  # Default species

    if not sequence or not isinstance(sequence, str):
        return jsonify({"error": "Invalid or missing 'sequence' field"}), 400

    if not authenticate(userid, authtoken):
        return jsonify({"Authentification App error": "Invalid or missing authentication"}), 400

    preprocessed = preprocessing(sequence_id, sequence, species=species, debug=debug)
    ab, e1 = antibody_identification(preprocessed, debug=debug)
    ab, e2 = optimize(ab, species=species, debug=debug)

    errors = e1 + e2

    keys_to_keep = ['sequence_id', 'user_input', 'optimized_vdj', 'optimizations', 'optimizations_count', 'optimized_gc_content', 'optimization_timestamp']
    response = {k: ab[k] for k in keys_to_keep}

    billing(user=userid, token=authtoken, app='Optimize')

    if debug:
        return jsonify({"result": response, "errors": errors}), 200
    else:
        return jsonify({"result": response}), 200



@app.route('/clone', methods=['POST'])
def clone():
    payload = request.get_json() or request.form

    userid = payload.get('userid')
    authtoken = payload.get('authtoken')
    debug = payload.get('debug', False)

    if not authenticate(userid, authtoken):
        return jsonify({"Authentification App error": "Invalid or missing authentication"}), 400

    ab, vector = payload.get('sequence'), payload.get('vector')
    ab, error = clone(ab=ab, vector=vector, debug=debug)

    billing(user=userid, token=authtoken, app='Clone')

    if debug:
        return jsonify({"result": ab, "errors": error}), 200
    else:
        return ab



@app.route('/number', methods=['POST'])
def number():
    raw = request.get_json() or request.form
    header, data = raw
    ab, numbering_scheme = data
    userid = header.get('userid')
    authtoken = header.get('authtoken')
    ab = numbering(ab, numbering_scheme)
    billing(user=userid, token=authtoken, app='Numbering')
    return ab



@app.route('/annotate', methods=['POST'])
def annotate():
    ab = request.get_json() or request.form
    ab = abnotator(ab)
    return ab



@app.route('/clonotype', methods=['POST'])
def clonotype():
    data = request.get_json() or request.form
    clusters = get_clusters(data, )
    return clusters



@app.route('/phylogeny', methods=['POST'])
def phylogeny():
    data = request.get_json() or request.form
    phylo = get_phylogeny(data, )
    return phylo



@app.route('/tree', methods=['POST'])
def tree():
    data = request.get_json() or request.form
    tree = build_tree(data, )
    return tree



# @app.route('/stabilize', methods=['POST'])
# def stabilize():
    # data = request.get_json() or request.form
    # stabilized = stabilize(data, )
    # return stabilized



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

        ab = preprocessing(sequence_id, sequence, species=species, debug=debug)
        ab, errors = antibody_identification(ab, debug=debug)

        keys_to_keep = ['sequence_id', 'user_input', 'humanized', 'humanization_score', 'humanization_mutations', 'humanization_percent_change']

        # Perform single humanization
        try:
            ab = single_humanize(ab=ab, temp=temp, debug=debug)
            response = {k: ab[k] for k in keys_to_keep}
            billing(user=userid, token=authtoken, app='Humanize (single)')
        except Exception as e:
            return jsonify({"Single Humanization App error": str(e)}), 500
        
        if debug == True:
            return jsonify({"result": response, "errors": errors}), 200
        else:
            return jsonify({"result": response}), 200

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

        ab = preprocessing(sequence_id, sequence, species=species, debug=debug)
        ab = antibody_identification(ab, debug=debug)

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
            billing(user=userid, token=authtoken, app='Humanize (multi)')
        except Exception as e:
            return jsonify({"Multi Humanize App error": str(e)}), 500

    else:
        return jsonify({"error": "Invalid model type. Must be 'single' or 'multi'"}), 400

    # Return the results
    return jsonify({"result": ab}), 200



@app.route('/bill', methods=['POST'])
def bill():
    data = request.get_json() or request.form
    userid = data.get('userid')
    authtoken = data.get('authtoken')
    
    if not authenticate(userid, authtoken):
        return jsonify({"Authentification App error": "Invalid or missing authentication"}), 400

    try:
        bill_df = get_bill(userid)

        # Convert the DataFrame to JSON-compatible dictionary
        bill_json = bill_df.to_dict(orient="records")  # Convert rows to a list of dictionaries
        return jsonify({"bill": bill_json}), 200
    except Exception as e:
        return jsonify({"Authentification App error": str(e)}), 500




@app.route("/graphql", methods=["POST"])
def graphql_server():
    data = request.get_json()
    success, result = graphql_sync(schema, data, context_value=request, debug=True)
    status_code = 200 if success else 400
    return jsonify(result), status_code



#############################################################################################################
#                                                                                                           #
#               API caller                                                                                  #
#                                                                                                           #
#############################################################################################################

if __name__ == '__main__':
    app.run(debug=True, port=8282)
