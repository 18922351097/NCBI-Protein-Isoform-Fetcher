import os
from flask import Flask, render_template, request, jsonify
from Bio import Entrez, SeqIO
from io import StringIO

app = Flask(__name__)
app.secret_key = os.environ.get("FLASK_SECRET_KEY") or "a secret key"

# Set your email for NCBI
Entrez.email = "your_email@example.com"

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/fetch_sequence', methods=['POST'])
def fetch_sequence():
    sequence_id = request.form.get('sequence_id')
    
    try:
        # Fetch the sequence from NCBI
        handle = Entrez.efetch(db="nucleotide", id=sequence_id, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        
        sequence_data = {
            'id': record.id,
            'description': record.description,
            'sequence': str(record.seq)
        }
        
        return jsonify(success=True, data=sequence_data)
    except Exception as e:
        return jsonify(success=False, error=str(e))

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)
