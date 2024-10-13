document.addEventListener('DOMContentLoaded', function() {
    const form = document.getElementById('sequenceForm');
    const resultDiv = document.getElementById('result');
    const sequenceDataDiv = document.getElementById('sequenceData');
    const errorDiv = document.getElementById('error');
    const loadingDiv = document.getElementById('loading');

    form.addEventListener('submit', function(e) {
        e.preventDefault();
        const queries = document.getElementById('queries').value.trim();
        if (validateInput(queries)) {
            fetchSequences(queries);
        }
    });

    function validateInput(input) {
        if (input === '') {
            displayError('Please enter at least one gene name.');
            return false;
        }
        return true;
    }

    function fetchSequences(queries) {
        showLoading(true);
        fetch('/fetch_sequence', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/x-www-form-urlencoded',
            },
            body: `queries=${encodeURIComponent(queries)}`
        })
        .then(response => response.json())
        .then(data => {
            showLoading(false);
            if (data.success) {
                displayResult(data.data, data.errors);
            } else {
                displayError(data.error, data.errors);
            }
        })
        .catch(error => {
            showLoading(false);
            displayError('An error occurred while fetching the sequences.', [error.message]);
        });
    }

    function displayResult(sequences, errors) {
        sequenceDataDiv.innerHTML = '';
        sequences.forEach((sequence, index) => {
            const sequenceElement = document.createElement('div');
            sequenceElement.className = 'mb-4';
            sequenceElement.innerHTML = `
                <h4>Gene: ${sequence.gene_name} (ID: ${sequence.gene_id})</h4>
                <p>NCBI Link: <a href="${sequence.ncbi_link}" target="_blank" class="text-info">${sequence.ncbi_link}</a></p>
                <p>Total Isoforms: ${sequence.total_isoforms}</p>
                <div id="proteinIsoforms-${index}"></div>
            `;
            sequenceDataDiv.appendChild(sequenceElement);

            displayProteinIsoforms(sequence.protein_variants, index);
        });

        if (errors && errors.length > 0) {
            displayErrorList(errors);
        }

        resultDiv.classList.remove('d-none');
        errorDiv.classList.add('d-none');
    }

    function displayProteinIsoforms(variants, sequenceIndex) {
        const isoformsContainer = document.getElementById(`proteinIsoforms-${sequenceIndex}`);
        if (!variants || variants.length === 0) {
            isoformsContainer.innerHTML = '<p>No protein isoforms found.</p>';
            return;
        }

        let isoformHtml = '';

        variants.forEach((variant, index) => {
            const functionalLabel = variant.functional ? '<span class="badge bg-success ms-2">Functional</span>' : '';
            isoformHtml += `
                <div id="protein-isoform-${sequenceIndex}-${index}" class="protein-isoform mb-4">
                    <h5>Isoform #${variant.index}: ${variant.id}${functionalLabel}</h5>
                    <p>Description: ${variant.description}</p>
                    <p>Size: ${variant.size} aa</p>
                    <p>NCBI Link: <a href="${variant.ncbi_link}" target="_blank" class="text-info">${variant.ncbi_link}</a></p>
                    <pre class="bg-dark text-light p-3 rounded mt-2">${variant.sequence}</pre>
                    <button class="btn btn-secondary mt-2" onclick="downloadSequence('${variant.sequence}', '${variant.id}_protein.txt')">Download Protein Sequence</button>
                </div>
            `;
        });

        isoformsContainer.innerHTML = isoformHtml;
    }

    function displayError(message, errors = []) {
        errorDiv.innerHTML = `<p>${message}</p>`;
        if (errors.length > 0) {
            displayErrorList(errors);
        }
        errorDiv.classList.remove('d-none');
        resultDiv.classList.add('d-none');
    }

    function displayErrorList(errors) {
        const errorList = document.createElement('ul');
        errorList.className = 'list-group mb-3';
        errors.forEach(error => {
            const errorItem = document.createElement('li');
            errorItem.className = 'list-group-item list-group-item-danger';
            errorItem.textContent = error;
            errorList.appendChild(errorItem);
        });
        errorDiv.appendChild(errorList);
    }

    function showLoading(show) {
        if (show) {
            loadingDiv.classList.remove('d-none');
            resultDiv.classList.add('d-none');
            errorDiv.classList.add('d-none');
        } else {
            loadingDiv.classList.add('d-none');
        }
    }
});

function downloadSequence(sequenceData, fileName) {
    const blob = new Blob([sequenceData], { type: 'text/plain' });
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.style.display = 'none';
    a.href = url;
    a.download = fileName;
    document.body.appendChild(a);
    a.click();
    window.URL.revokeObjectURL(url);
}
