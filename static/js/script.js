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
            displayError('Please enter at least one protein ID or gene name.');
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
                <h4>ID: ${sequence.id}</h4>
                <p>NCBI Link: <a href="${sequence.ncbi_link}" target="_blank" class="text-info">${sequence.ncbi_link}</a></p>
                <div id="proteinVariants-${index}"></div>
            `;
            sequenceDataDiv.appendChild(sequenceElement);

            displayProteinVariants(sequence.protein_variants, index);
        });

        if (errors && errors.length > 0) {
            displayErrorList(errors);
        }

        resultDiv.classList.remove('d-none');
        errorDiv.classList.add('d-none');
    }

    function displayProteinVariants(variants, sequenceIndex) {
        const variantsContainer = document.getElementById(`proteinVariants-${sequenceIndex}`);
        if (!variants || variants.length === 0) {
            variantsContainer.innerHTML = '<p>No protein variants found.</p>';
            return;
        }

        let variantHtml = `
            <select class="form-select mb-2" id="protein-variant-select-${sequenceIndex}" onchange="showSelectedProteinVariant(${sequenceIndex})">
                ${variants.map((variant, index) => `
                    <option value="${index}">${variant.label}: ${variant.id}</option>
                `).join('')}
            </select>
        `;

        variants.forEach((variant, index) => {
            variantHtml += `
                <div id="protein-variant-${sequenceIndex}-${index}" class="protein-variant ${index === 0 ? '' : 'd-none'}">
                    <h5>${variant.label}: ${variant.id}</h5>
                    <p>Description: ${variant.description}</p>
                    <pre class="bg-dark text-light p-3 rounded mt-2">${variant.sequence}</pre>
                    <button class="btn btn-secondary mt-2" onclick="downloadSequence('${variant.sequence}', '${variant.id}_protein.txt')">Download Protein Sequence</button>
                </div>
            `;
        });

        variantsContainer.innerHTML = variantHtml;
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

function showSelectedProteinVariant(sequenceIndex) {
    const select = document.getElementById(`protein-variant-select-${sequenceIndex}`);
    const selectedIndex = select.value;
    const variants = document.querySelectorAll(`#proteinVariants-${sequenceIndex} .protein-variant`);
    variants.forEach((variant, index) => {
        if (index.toString() === selectedIndex) {
            variant.classList.remove('d-none');
        } else {
            variant.classList.add('d-none');
        }
    });
}
