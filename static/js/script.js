document.addEventListener('DOMContentLoaded', function() {
    const form = document.getElementById('sequenceForm');
    const resultDiv = document.getElementById('result');
    const sequenceDataDiv = document.getElementById('sequenceData');
    const errorDiv = document.getElementById('error');
    const downloadBtn = document.getElementById('downloadBtn');
    const loadingDiv = document.getElementById('loading');

    form.addEventListener('submit', function(e) {
        e.preventDefault();
        const queries = document.getElementById('queries').value.trim();
        if (validateInput(queries)) {
            fetchSequences(queries);
        }
    });

    downloadBtn.addEventListener('click', function() {
        downloadSequences();
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
            const sequenceElement = document.createElement('pre');
            sequenceElement.className = 'bg-dark text-light p-3 rounded mb-3';
            sequenceElement.textContent = `ID: ${sequence.id}\nDescription: ${sequence.description}\n\nSequence:\n${sequence.sequence}`;
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

    function downloadSequences() {
        const sequenceData = sequenceDataDiv.innerText;
        const blob = new Blob([sequenceData], { type: 'text/plain' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = 'sequences.txt';
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
    }
});
