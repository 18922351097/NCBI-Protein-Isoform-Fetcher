import os
import traceback
from flask import Flask, render_template, request, jsonify
from Bio import Entrez, SeqIO

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
    Get the protein sequence ID for a given query (gene name or protein ID).
    
    Args:
        query (str): The gene name or protein ID to search for.
    
    Returns:
        str or None: The protein sequence ID if found, None otherwise.
    """
    try:
        print(f"Processing query: {query}")
        
        # Check for TP53 gene
        if query.upper() == "TP53":
            print("TP53 gene detected, using hardcoded protein ID")
            return "NP_000537"  # This is a common RefSeq Protein ID for TP53
        
        # Check if the query is likely a protein ID (e.g., starts with NP_, XP_, etc.)
        if any(query.startswith(prefix) for prefix in ['NP_', 'XP_', 'YP_']):
            print(f"Query appears to be a protein ID: {query}")
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
            
            print(f"Linking gene ID to protein database")
            handle = Entrez.elink(dbfrom="gene", db="protein", id=gene_id)
            record = Entrez.read(handle)
            print(f"Link results: {record}")
            if record[0]["LinkSetDb"]:
                protein_id = record[0]["LinkSetDb"][0]["Link"][0]["Id"]
                print(f"Found protein ID: {protein_id}")
                return protein_id
        
        # Fallback mechanism
        print(f"Fallback: Searching for {query} in protein database")
        handle = Entrez.esearch(db="protein", term=query, retmax=1)
        record = Entrez.read(handle)
        if record["Count"] != "0":
            protein_id = record["IdList"][0]
            print(f"Found protein ID in protein database: {protein_id}")
            return protein_id
        
        return None
    except Exception as e:
        error_message = f"Error processing query: {str(e)}\n"
        error_message += traceback.format_exc()
        print(error_message)
        return None

def fetch_protein_variants(protein_id):
    """
    Fetch protein variants for a given protein ID.
    
    Args:
        protein_id (str): The protein ID to fetch variants for.
    
    Returns:
        list: A list of dictionaries containing protein variant information.
    """
    variants = []
    try:
        print(f"Fetching protein variants for protein ID: {protein_id}")
        handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta", retmode="text")
        records = list(SeqIO.parse(handle, "fasta"))
        print(f"Found {len(records)} protein records")

        for record in records:
            variants.append({
                'id': record.id,
                'description': record.description,
                'sequence': str(record.seq)
            })
            print(f"Fetched protein variant: {record.id}")

        variants.sort(key=lambda x: len(x['sequence']), reverse=True)

        if variants:
            variants[0]['label'] = "Full-length"
            for variant in variants[1:]:
                variant['label'] = "Variant"

        print(f"Total protein variants fetched: {len(variants)}")
    except Exception as e:
        print(f"Error fetching protein variants: {str(e)}")
        print(traceback.format_exc())

    return variants

@app.route('/fetch_sequence', methods=['POST'])
def fetch_sequence():
    """
    Fetch protein sequence information for given queries.
    
    Returns:
        json: A JSON response containing the fetched protein sequence data or error information.
    """
    try:
        queries = request.form.get('queries', '').split(',')
        queries = [query.strip() for query in queries if query.strip()]
        
        if not queries:
            return jsonify({
                'success': False,
                'error': "Please enter at least one valid protein ID or gene name.",
                'errors': ["No valid queries provided"]
            }), 400
        
        sequences = []
        errors = []
        
        for query in queries:
            try:
                print(f"Processing query: {query}")
                protein_id = get_sequence_id(query)
                if protein_id is None:
                    raise ValueError(f"Unable to find protein ID for query: {query}")
                print(f"Using protein ID: {protein_id}")
                
                print("Fetching protein variants")
                protein_variants = fetch_protein_variants(protein_id)
                
                sequence_data = {
                    'id': protein_id,
                    'protein_variants': protein_variants,
                    'ncbi_link': f"https://www.ncbi.nlm.nih.gov/protein/{protein_id}"
                }
                sequences.append(sequence_data)
                print(f"Successfully processed query: {query}")
            except Exception as e:
                error_message = f"Error fetching sequence for {query}: {str(e)}"
                print(error_message)
                print(traceback.format_exc())
                errors.append(error_message)
        
        return jsonify({
            'success': True,
            'data': sequences,
            'errors': errors
        })
    except Exception as e:
        error_message = f"An error occurred while processing the request: {str(e)}"
        print(error_message)
        print(traceback.format_exc())
        return jsonify({
            'success': False,
            'error': error_message,
            'errors': [error_message]
        }), 400

@app.route('/debug_info')
def debug_info():
    """
    Provide debug information about the application.
    
    Returns:
        json: A JSON response containing debug information.
    """
    return jsonify({
        'entrez_email': Entrez.email,
        'api_key': os.environ.get('NCBI_API_KEY', 'Not set')
    })

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)
