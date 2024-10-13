document.addEventListener('DOMContentLoaded', function() {
    const form = document.getElementById('sequenceForm');
    const resultDiv = document.getElementById('result');
    const sequenceDataPre = document.getElementById('sequenceData');
    const errorDiv = document.getElementById('error');
    const downloadBtn = document.getElementById('downloadBtn');

    form.addEventListener('submit', function(e) {
        e.preventDefault();
        const sequenceId = document.getElementById('sequenceId').value;
        fetchSequence(sequenceId);
    });

    downloadBtn.addEventListener('click', function() {
        downloadSequence();
    });

    function fetchSequence(sequenceId) {
        fetch('/fetch_sequence', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/x-www-form-urlencoded',
            },
            body: `sequence_id=${encodeURIComponent(sequenceId)}`
        })
        .then(response => response.json())
        .then(data => {
            if (data.success) {
                displayResult(data.data);
            } else {
                displayError(data.error);
            }
        })
        .catch(error => {
            displayError('An error occurred while fetching the sequence.');
        });
    }

    function displayResult(data) {
        sequenceDataPre.textContent = `ID: ${data.id}\nDescription: ${data.description}\n\nSequence:\n${data.sequence}`;
        resultDiv.classList.remove('d-none');
        errorDiv.classList.add('d-none');
    }

    function displayError(message) {
        errorDiv.textContent = message;
        errorDiv.classList.remove('d-none');
        resultDiv.classList.add('d-none');
    }

    function downloadSequence() {
        const sequenceData = sequenceDataPre.textContent;
        const blob = new Blob([sequenceData], { type: 'text/plain' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = 'sequence.txt';
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
    }
});
