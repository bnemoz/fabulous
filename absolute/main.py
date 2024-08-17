from flask import Flask, request, jsonify
from flask_cors import CORS
import abstar

app = Flask(__name__)
cors = CORS(app, origins=['*'], allow_headers=['Content-Type', 'Access-Control-Allow-Origin', 'Access-Control-Allow-Headers', 'Access-Control-Allow-Methods'])


def antibody_identification(sequence, species):
    ab_json = abstar.run(sequence, germ_db=species, verbose=False)
    return ab_json


@app.route('/', methods=['GET'])
def index():
    return "Welcome to the Antibody Identification API!"


@app.route('/id', methods=['GET'])
def id():
    sequence = request.args.get('sequence')

    if sequence == None:
        ab_json = {}
        return ab_json


    ab_json = antibody_identification(sequence, species='human')
    return dict(ab_json)


@app.route('/clone', methods=['GET', 'POST'])
def clone():
    sequence = request.args.get('sequence')
    clone = {}
    clone['original'] = sequence
    # clone['optimized'] = 
    return clone


if __name__ == '__main__':
    app.run(debug=True, port=8282)
