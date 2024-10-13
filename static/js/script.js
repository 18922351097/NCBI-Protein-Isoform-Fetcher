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
            displayError('Please enter at least one sequence ID or gene name.');
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
                displayError(data.error);
            }
        })
        .catch(error => {
            showLoading(false);
            displayError('An error occurred while fetching the sequences.');
        });
    }

    function displayResult(sequences, errors) {
        sequenceDataDiv.innerHTML = '';
        sequences.forEach((sequence, index) => {
            const sequenceElement = document.createElement('div');
            sequenceElement.className = 'mb-4';
            sequenceElement.innerHTML = `
                <h4>ID: ${sequence.id}</h4>
                <p>Description: ${sequence.description}</p>
                <p>NCBI Link: <a href="${sequence.ncbi_link}" target="_blank" class="text-info">${sequence.ncbi_link}</a></p>
                <ul class="nav nav-tabs" id="sequenceTabs-${index}" role="tablist">
                    <li class="nav-item" role="presentation">
                        <button class="nav-link active" id="dna-tab-${index}" data-bs-toggle="tab" data-bs-target="#dna-${index}" type="button" role="tab">DNA</button>
                    </li>
                    <li class="nav-item" role="presentation">
                        <button class="nav-link" id="rna-tab-${index}" data-bs-toggle="tab" data-bs-target="#rna-${index}" type="button" role="tab">RNA</button>
                    </li>
                    <li class="nav-item" role="presentation">
                        <button class="nav-link" id="protein-tab-${index}" data-bs-toggle="tab" data-bs-target="#protein-${index}" type="button" role="tab">Protein</button>
                    </li>
                </ul>
                <div class="tab-content" id="sequenceTabContent-${index}">
                    <div class="tab-pane fade show active" id="dna-${index}" role="tabpanel">
                        <pre class="bg-dark text-light p-3 rounded mt-2">${sequence.dna_sequence}</pre>
                    </div>
                    <div class="tab-pane fade" id="rna-${index}" role="tabpanel">
                        <pre class="bg-dark text-light p-3 rounded mt-2">${sequence.rna_sequence}</pre>
                    </div>
                    <div class="tab-pane fade" id="protein-${index}" role="tabpanel">
                        <pre class="bg-dark text-light p-3 rounded mt-2">${sequence.protein_sequence}</pre>
                    </div>
                </div>
            `;
            sequenceDataDiv.appendChild(sequenceElement);
        });

        if (errors && errors.length > 0) {
            const errorList = document.createElement('ul');
            errorList.className = 'list-group mb-3';
            errors.forEach(error => {
                const errorItem = document.createElement('li');
                errorItem.className = 'list-group-item list-group-item-danger';
                errorItem.textContent = error;
                errorList.appendChild(errorItem);
            });
            sequenceDataDiv.appendChild(errorList);
        }

        resultDiv.classList.remove('d-none');
        errorDiv.classList.add('d-none');
    }

    function displayError(message) {
        errorDiv.textContent = message;
        errorDiv.classList.remove('d-none');
        resultDiv.classList.add('d-none');
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
