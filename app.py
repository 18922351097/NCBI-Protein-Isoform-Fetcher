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

def get_sequence_id(query):
    try:
        handle = Entrez.esearch(db="gene", term=query + "[Gene Name]", retmax=1)
        record = Entrez.read(handle)
        if record["Count"] == "0":
            return None
        gene_id = record["IdList"][0]
        
        handle = Entrez.elink(dbfrom="gene", db="nucleotide", id=gene_id)
        record = Entrez.read(handle)
        if not record[0]["LinkSetDb"]:
            return None
        return record[0]["LinkSetDb"][0]["Link"][0]["Id"]
    except Exception as e:
        print(f"Error searching for gene: {str(e)}")
        return None

@app.route('/fetch_sequence', methods=['POST'])
def fetch_sequence():
    queries = request.form.get('queries', '').split(',')
    queries = [query.strip() for query in queries if query.strip()]
    
    if not queries:
        return jsonify(success=False, error="Please enter at least one valid sequence ID or gene name.")
    
    results = []
    errors = []
    
    for query in queries:
        try:
            sequence_id = get_sequence_id(query)
            if sequence_id is None:
                sequence_id = query
            
            # Fetch DNA sequence
            dna_handle = Entrez.efetch(db="nucleotide", id=sequence_id, rettype="fasta", retmode="text")
            dna_record = SeqIO.read(dna_handle, "fasta")
            
            # Fetch RNA sequence
            rna_handle = Entrez.efetch(db="nucleotide", id=sequence_id, rettype="fasta_cds_na", retmode="text")
            rna_record = SeqIO.read(rna_handle, "fasta")
            
            # Fetch protein sequence
            protein_handle = Entrez.efetch(db="protein", id=sequence_id, rettype="fasta", retmode="text")
            protein_record = SeqIO.read(protein_handle, "fasta")
            
            sequence_data = {
                'id': dna_record.id,
                'description': dna_record.description,
                'dna_sequence': str(dna_record.seq),
                'rna_sequence': str(rna_record.seq),
                'protein_sequence': str(protein_record.seq),
                'ncbi_link': f"https://www.ncbi.nlm.nih.gov/nuccore/{sequence_id}"
            }
            results.append(sequence_data)
        except Exception as e:
            errors.append(f"Error fetching sequence for {query}: {str(e)}")
    
    if results:
        return jsonify(success=True, data=results, errors=errors)
    else:
        return jsonify(success=False, error="Unable to fetch any sequences. Please check the IDs or gene names and try again.", errors=errors)

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)
