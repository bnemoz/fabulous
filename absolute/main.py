from flask import Flask, request, jsonify
from flask_cors import CORS
import abstar

app = Flask(__name__)
cors = CORS(app, origins=['*'], )

@app.route('/api/antibody', methods=['GET', 'POST'])
def main():
    sequence = request.args.get('sequence')
    return jsonify(
        {
        abstar.run(sequence, verbose=False)
        }
    )

if __name__ == '__main__':
    app.run(debug=True, port=8282)