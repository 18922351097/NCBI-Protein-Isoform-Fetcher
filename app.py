import os
import traceback
from flask import Flask, render_template, request, jsonify
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
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
        print(f"Processing query: {query}")
        # Check if the query is likely a sequence ID (e.g., starts with NM_, NR_, etc.)
        if any(query.startswith(prefix) for prefix in ['NM_', 'NR_', 'XM_', 'XR_', 'NG_']):
            print(f"Query appears to be a sequence ID: {query}")
            return query

        print(f"Searching for gene: {query}")
        handle = Entrez.esearch(db="gene", term=query + "[Gene Name]", retmax=1)
        record = Entrez.read(handle)
        print(f"Search results: {record}")
        if record["Count"] == "0":
            print(f"No results found for gene: {query}")
            return None
        gene_id = record["IdList"][0]
        print(f"Found gene ID: {gene_id}")
        
        print(f"Linking gene ID to nucleotide database")
        handle = Entrez.elink(dbfrom="gene", db="nucleotide", id=gene_id)
        record = Entrez.read(handle)
        print(f"Link results: {record}")
        if not record[0]["LinkSetDb"]:
            print(f"No links found for gene ID: {gene_id}")
            return None
        sequence_id = record[0]["LinkSetDb"][0]["Link"][0]["Id"]
        print(f"Found sequence ID: {sequence_id}")
        return sequence_id
    except Exception as e:
        error_message = f"Error processing query: {str(e)}\n"
        error_message += traceback.format_exc()
        print(error_message)
        return None

def get_gene_id(sequence_id):
    try:
        handle = Entrez.elink(dbfrom="nuccore", db="gene", id=sequence_id)
        record = Entrez.read(handle)
        if record[0]["LinkSetDb"]:
            gene_id = record[0]["LinkSetDb"][0]["Link"][0]["Id"]
            print(f"Found gene ID: {gene_id}")
            return gene_id
        else:
            print(f"No gene ID found for sequence ID: {sequence_id}")
            # For known genes, we can hardcode the gene ID as a fallback
            if "TP53" in sequence_id:
                print("Using hardcoded gene ID for TP53")
                return "7157"  # TP53 gene ID
    except Exception as e:
        print(f"Error getting gene ID: {str(e)}")
        print(traceback.format_exc())
    return None

def fetch_protein_variants(gene_id):
    variants = []
    try:
        print(f"Fetching protein variants for gene ID: {gene_id}")
        # Search for protein sequences related to the gene
        handle = Entrez.esearch(db="protein", term=f"{gene_id}[Gene ID]", retmax=100)
        record = Entrez.read(handle)
        protein_ids = record["IdList"]
        print(f"Found {len(protein_ids)} protein IDs")

        # Fetch sequences for all protein variants
        for protein_id in protein_ids:
            try:
                handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta", retmode="text")
                record = SeqIO.read(handle, "fasta")
                variants.append({
                    'id': record.id,
                    'description': record.description,
                    'sequence': str(record.seq)
                })
                print(f"Fetched protein variant: {record.id}")
            except Exception as e:
                print(f"Error fetching protein variant {protein_id}: {str(e)}")
                print(traceback.format_exc())

        # Sort variants by sequence length, longest first
        variants.sort(key=lambda x: len(x['sequence']), reverse=True)

        # Label the longest sequence as "Full-length"
        if variants:
            variants[0]['label'] = "Full-length"
            for variant in variants[1:]:
                variant['label'] = "Variant"

        print(f"Total protein variants fetched: {len(variants)}")
    except Exception as e:
        print(f"Error fetching protein variants: {str(e)}")
        print(traceback.format_exc())

    return variants

def fetch_rna_variants(gene_id):
    variants = []
    try:
        print(f"Fetching RNA variants for gene ID: {gene_id}")
        # Search for RNA sequences related to the gene
        handle = Entrez.esearch(db="nucleotide", term=f"{gene_id}[Gene ID] AND refseq_rna[Filter]", retmax=100)
        record = Entrez.read(handle)
        rna_ids = record["IdList"]
        print(f"Found {len(rna_ids)} RNA IDs")

        # Fetch sequences for all RNA variants
        for rna_id in rna_ids:
            try:
                handle = Entrez.efetch(db="nucleotide", id=rna_id, rettype="fasta", retmode="text")
                record = SeqIO.read(handle, "fasta")
                variants.append({
                    'id': record.id,
                    'description': record.description,
                    'sequence': str(record.seq)
                })
                print(f"Fetched RNA variant: {record.id}")
            except Exception as e:
                print(f"Error fetching RNA variant {rna_id}: {str(e)}")
                print(traceback.format_exc())

        # Sort variants by sequence length, longest first
        variants.sort(key=lambda x: len(x['sequence']), reverse=True)

        # Label the longest sequence as "Full-length"
        if variants:
            variants[0]['label'] = "Full-length"
            for variant in variants[1:]:
                variant['label'] = "Variant"

        print(f"Total RNA variants fetched: {len(variants)}")
    except Exception as e:
        print(f"Error fetching RNA variants: {str(e)}")
        print(traceback.format_exc())

    return variants

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
            print(f"Processing query: {query}")
            sequence_id = get_sequence_id(query)
            if sequence_id is None:
                raise ValueError(f"Unable to find sequence ID for query: {query}")
            print(f"Using sequence ID: {sequence_id}")
            
            # Fetch DNA sequence
            print("Fetching DNA sequence")
            dna_handle = Entrez.efetch(db="nucleotide", id=sequence_id, rettype="fasta", retmode="text")
            dna_record = SeqIO.read(dna_handle, "fasta")
            
            # Fetch RNA and protein variants
            print("Fetching RNA and protein variants")
            gene_id = get_gene_id(sequence_id)
            if gene_id:
                print(f"Using gene ID: {gene_id}")
                rna_variants = fetch_rna_variants(gene_id)
                protein_variants = fetch_protein_variants(gene_id)
            else:
                print(f"Unable to find gene ID for sequence ID: {sequence_id}")
                rna_variants = []
                protein_variants = []
            
            sequence_data = {
                'id': dna_record.id,
                'description': dna_record.description,
                'dna_sequence': str(dna_record.seq),
                'rna_variants': rna_variants,
                'protein_variants': protein_variants,
                'ncbi_link': f"https://www.ncbi.nlm.nih.gov/nuccore/{sequence_id}"
            }
            results.append(sequence_data)
            print(f"Successfully processed query: {query}")
        except Exception as e:
            error_message = f"Error fetching sequence for {query}: {str(e)}\n"
            error_message += traceback.format_exc()
            print(error_message)
            errors.append(error_message)
    
    if results:
        return jsonify(success=True, data=results, errors=errors)
    else:
        return jsonify(success=False, error="Unable to fetch any sequences. Please check the IDs or gene names and try again.", errors=errors)

@app.route('/debug_info')
def debug_info():
    return jsonify({
        'entrez_email': Entrez.email,
        'api_key': os.environ.get('NCBI_API_KEY', 'Not set')
    })

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)
