# NCBI Protein Isoform Fetcher

## Introduction
The NCBI Protein Isoform Fetcher is a web application that allows users to fetch protein isoforms for given gene names or IDs from the National Center for Biotechnology Information (NCBI) database. This tool is particularly useful for researchers and bioinformaticians who need quick access to protein sequence data for various genes.

## Features
- Fetch protein isoforms for multiple genes simultaneously
- Display detailed information about each isoform, including:
  - Protein ID
  - Description
  - Sequence length
  - Full amino acid sequence
- Provide direct links to NCBI entries for further information
- Allow downloading of protein sequences in plain text format
- Responsive web interface with a dark theme for comfortable viewing

## Installation

### Prerequisites
- Python 3.x
- Flask
- Biopython

### Steps
1. Clone the repository:
   ```
   git clone https://github.com/your-username/ncbi-protein-isoform-fetcher.git
   cd ncbi-protein-isoform-fetcher
   ```

2. Install the required packages:
   ```
   pip install flask biopython
   ```

3. Set up the environment variables:
   - Create a `.env` file in the project root and add the following:
     ```
     FLASK_SECRET_KEY=your_secret_key_here
     NCBI_API_KEY=your_ncbi_api_key_here
     ```

## Usage

1. Start the Flask server:
   ```
   python main.py
   ```

2. Open a web browser and navigate to `http://localhost:5000`

3. Enter one or more gene names or IDs in the input field, separated by commas

4. Click the "Fetch Protein Isoforms" button

5. View the results, including protein isoform details and sequences

6. (Optional) Download individual protein sequences using the provided buttons

## API Endpoints

- `/`: Main page (GET)
- `/fetch_sequence`: Fetch protein sequences (POST)
- `/debug_info`: Get debug information (GET)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is open source and available under the [MIT License](LICENSE).

## Acknowledgements

- [NCBI](https://www.ncbi.nlm.nih.gov/) for providing the biological data
- [Biopython](https://biopython.org/) for the excellent bioinformatics tools
- [Flask](https://flask.palletsprojects.com/) for the web framework
- [Bootstrap](https://getbootstrap.com/) for the frontend design

## Contact

For any questions or concerns, please open an issue on the GitHub repository.
