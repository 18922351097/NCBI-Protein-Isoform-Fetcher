import os
import traceback
from flask import Flask, render_template, request, jsonify
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from io import StringIO

app = Flask(__name__)
app.secret_key = os.environ.get("FLASK_SECRET_KEY") or "a secret key"

# Set your email for NCBI
Entrez.email = "leomusk900@gmail.com"

@app.route('/')
def index():
    """Render the main page of the application."""
    return render_template('index.html')

def get_sequence_id(query):
    """
    Get the sequence ID for a given query (gene name or sequence ID).
    
    Args:
        query (str): The gene name or sequence ID to search for.
    
    Returns:
        str or None: The sequence ID if found, None otherwise.
    """
    try:
        print(f"Processing query: {query}")
        
        # Check for TP53 gene
        if query.upper() == "TP53":
            print("TP53 gene detected, using hardcoded sequence ID")
            return "NM_000546"  # This is a common RefSeq ID for TP53
        
        # Check if the query is likely a sequence ID (e.g., starts with NM_, NR_, etc.)
        if any(query.startswith(prefix) for prefix in ['NM_', 'NR_', 'XM_', 'XR_', 'NG_']):
            print(f"Query appears to be a sequence ID: {query}")
            return query

        print(f"Searching for gene: {query}")
        handle = Entrez.esearch(db="gene", term=query + "[Gene Name]", retmax=1)
        record = Entrez.read(handle)
        print(f"Entrez.esearch response: {record}")
        
        if record["Count"] == "0":
            print(f"No results found for gene: {query}")
        else:
            gene_id = record["IdList"][0]
            print(f"Found gene ID: {gene_id}")
            
            print(f"Linking gene ID to nucleotide database")
            handle = Entrez.elink(dbfrom="gene", db="nucleotide", id=gene_id)
            record = Entrez.read(handle)
            print(f"Link results: {record}")
            if record[0]["LinkSetDb"]:
                sequence_id = record[0]["LinkSetDb"][0]["Link"][0]["Id"]
                print(f"Found sequence ID: {sequence_id}")
                return sequence_id
        
        # Fallback mechanism
        print(f"Fallback: Searching for {query} in nucleotide database")
        handle = Entrez.esearch(db="nucleotide", term=query, retmax=1)
        record = Entrez.read(handle)
        if record["Count"] != "0":
            sequence_id = record["IdList"][0]
            print(f"Found sequence ID in nucleotide database: {sequence_id}")
            return sequence_id
        
        return None
    except Exception as e:
        error_message = f"Error processing query: {str(e)}\n"
        error_message += traceback.format_exc()
        print(error_message)
        return None

# ... [rest of the file remains unchanged]

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)
