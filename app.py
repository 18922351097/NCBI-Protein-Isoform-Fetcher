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

def get_gene_id(query):
    """
    Get the gene ID for a given query (gene name or protein ID).
    
    Args:
        query (str): The gene name or protein ID to search for.
    
    Returns:
        str or None: The gene ID if found, None otherwise.
    """
    try:
        print(f"Processing query: {query}")
        
        # Check for TP53 gene
        if query.upper() == "TP53":
            print("TP53 gene detected, using hardcoded gene ID")
            return "7157"  # This is the gene ID for TP53
        
        print(f"Searching for gene: {query}")
        handle = Entrez.esearch(db="gene", term=f"{query}[Gene Name] AND human[ORGN]", retmax=1)
        record = Entrez.read(handle)
        print(f"Entrez.esearch response: {record}")
        
        if record["Count"] == "0":
            print(f"No results found for gene: {query}")
            return None
        else:
            gene_id = record["IdList"][0]
            print(f"Found gene ID: {gene_id}")
            return gene_id
        
    except Exception as e:
        error_message = f"Error processing query: {str(e)}\n"
        error_message += traceback.format_exc()
        print(error_message)
        return None

def fetch_protein_variants(gene_id):
    """
    Fetch protein isoforms for a given gene ID.
    
    Args:
        gene_id (str): The gene ID to fetch protein isoforms for.
    
    Returns:
        list: A list of dictionaries containing protein isoform information.
    """
    variants = []
    try:
        print(f"Fetching protein isoforms for gene ID: {gene_id}")
        handle = Entrez.esearch(db="protein", term=f"{gene_id}[Gene ID] AND refseq[filter]", retmax=100)
        record = Entrez.read(handle)
        protein_ids = record["IdList"]
        print(f"Found {len(protein_ids)} protein isoforms")

        for protein_id in protein_ids:
            handle = Entrez.efetch(db="protein", id=protein_id, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            
            # Check for indicators of functionality
            functional = False
            if record.features:
                for feature in record.features:
                    if feature.type in ["CDS", "mat_peptide", "domain"]:
                        functional = True
                        break
            
            if "isoform" in record.description.lower() or "variant" in record.description.lower() or "functional" in record.description.lower():
                functional = True

            variants.append({
                'id': record.id,
                'description': record.description,
                'sequence': str(record.seq),
                'size': len(record.seq),
                'ncbi_link': f"https://www.ncbi.nlm.nih.gov/protein/{record.id}",
                'functional': functional
            })
            print(f"Fetched protein isoform: {record.id}, Functional: {functional}")

        variants.sort(key=lambda x: len(x['sequence']), reverse=True)
        if variants:
            variants[0]['label'] = "Longest Isoform"
            for variant in variants[1:]:
                variant['label'] = "Isoform"

        print(f"Total protein isoforms fetched: {len(variants)}")
    except Exception as e:
        print(f"Error fetching protein isoforms: {str(e)}")
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
                'error': "Please enter at least one valid gene name.",
                'errors': ["No valid queries provided"]
            }), 400
        
        sequences = []
        errors = []
        
        for query in queries:
            try:
                print(f"Processing query: {query}")
                gene_id = get_gene_id(query)
                if gene_id is None:
                    raise ValueError(f"Unable to find gene ID for query: {query}")
                print(f"Using gene ID: {gene_id}")
                
                print("Fetching protein isoforms")
                protein_variants = fetch_protein_variants(gene_id)
                
                sequence_data = {
                    'gene_name': query,
                    'gene_id': gene_id,
                    'protein_variants': protein_variants,
                    'ncbi_link': f"https://www.ncbi.nlm.nih.gov/gene/{gene_id}"
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
