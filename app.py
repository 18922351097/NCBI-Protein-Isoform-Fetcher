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
    sequence_ids = request.form.get('sequence_ids', '').split(',')
    sequence_ids = [seq_id.strip() for seq_id in sequence_ids if seq_id.strip()]
    
    if not sequence_ids:
        return jsonify(success=False, error="Please enter at least one valid sequence ID.")
    
    results = []
    errors = []
    
    for sequence_id in sequence_ids:
        try:
            # Fetch the sequence from NCBI
            handle = Entrez.efetch(db="nucleotide", id=sequence_id, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            
            sequence_data = {
                'id': record.id,
                'description': record.description,
                'sequence': str(record.seq)
            }
            results.append(sequence_data)
        except Exception as e:
            errors.append(f"Error fetching sequence {sequence_id}: {str(e)}")
    
    if results:
        return jsonify(success=True, data=results, errors=errors)
    else:
        return jsonify(success=False, error="Unable to fetch any sequences. Please check the IDs and try again.", errors=errors)

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)
