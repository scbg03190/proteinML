let index = 0;
function addSequence(sequence) {
    index += 1;
    let element = document.createElement("tr")
    element.innerHTML = `
    <tr>
    <td>${index}</td>
    <td>
        <input type="text" name="mrna" id="mrna" placeholder="mRNA Sequence" value="${sequence}" >
    </td>
    <td>
        <input type="number" name="start-codon" id="startcodon${index}" placeholder="0" required >
    </td>
    <td>
        <input type="number" name="stop-codon" id="stopcodon${index}" placeholder="0" required >
    </td>
    <td class="text-center">
    </td>
</tr>`
    document.getElementById("sequences-table").appendChild(element);
}

document.getElementById("sequence-adder").addEventListener("click", _ => {
    addSequence("")
})

document.getElementById("file").addEventListener("change", _ => {
    let file = document.getElementById("file").files[0];
    let reader = new FileReader();
    reader.readAsText(file, "UTF-8");
    reader.onload = function (e) {
        let data = e.target.result.trim();
        let sequences = [];
        let lines = data.split("\n");
        for (let i = 0; i < lines.length; i++) {
            if (lines[i].startsWith(">")) {
                sequences.push("");
            } else {
                sequences[sequences.length - 1] += lines[i];
            }
        }
        for (let j = 0; j < sequences.length; j++) {
            addSequence(sequences[j]);
        }
    }
})

function is_valid_sequence(gene_sequence, start_codon_index, stop_codon_index) {
    const start_codon = "AUG";
    const stop_codons = ["UAA", "UAG", "UGA"];

    if (
        gene_sequence.slice(start_codon_index - 1, start_codon_index + 2) === start_codon &&
        stop_codons.includes(gene_sequence.slice(stop_codon_index - 1, stop_codon_index + 2))
    ) {
        const cds_sequence = gene_sequence.slice(start_codon_index + 2, stop_codon_index - 1);

        if (cds_sequence.length % 3 !== 0) {
            return "Please make sure coding sequence is in triplets";
        }

        return "✅";
    } else {
        return "Invalid start or stop codon index.";
    }
}

// ${is_valid_sequence(
//     sequence,
//     Number(document.getElementById("startcodon").value),
//     Number(document.getElementById("stopcodon").value)
// )}
document.getElementById("submit-btn").addEventListener("click", _ => {
    let sequencesData = [...document.getElementById("sequences-table").querySelectorAll("tr")].map(tr => {
        let cells = tr.cells;
        console.log(cells)
        cells[4].innerText = is_valid_sequence(
                cells[1].firstElementChild.value,
                Number(cells[2].firstElementChild.value),
                Number(cells[3].firstElementChild.value)
            )
        return [cells[1].firstElementChild.value, Number(cells[2].firstElementChild.value), Number(cells[3].firstElementChild.value)]
    })
    fetch("/initiation_rate_prediction", {
        method: "POST",
        headers: {
            "Content-Type": "application/json"
        },
        body: JSON.stringify(sequencesData)
    }).then(response => response.json()).then(data => {
        let outputBody = document.getElementById("output-body")
        outputBody.innerHTML = ""
        for (let i = 0; i < data.length; i++) {
            let k = data[i];
            let childCode = `
            <tr>
            <td>${i + 1}</td>
            <td>${k.N1}</td>
            <td>${k.N4}</td>
            <td>${k.folding_energy_70}</td>
            <td>${k.folding_energy_80}</td>
            <td>${k.gene_length}</td>
            <td>${k['in_frame AUG']}</td>
            <td>${k.kozak_score}</td>
            <td>${k.length_of_5prime_utr}</td>
            <td>${k.initiation_rate}</td>
        </tr>`
            let child = document.createElement("tr");
            child.innerHTML = childCode;
            outputBody.appendChild(child)
        }
    })
})


document.getElementById("optimize").addEventListener("click", e => {
    let targetI = Number(document.getElementById("target-i").value)
    let iterations = Number(document.getElementById("iterations").value)
    let method = Number(document.getElementById("method").value)
    console.log(method)
    let listOfSequences = [...document.getElementById("sequences-table").querySelectorAll("tr")].map(tr => {
        let cells = tr.cells;
        return {
            value: cells[1].firstElementChild.value,
            start_codon_index: Number(cells[2].firstElementChild.value),
            stop_codon_index: Number(cells[3].firstElementChild.value),
            five_prime_utr: null,
        };
    });
    
    [...document.getElementById("output-body").querySelectorAll("tr")].map((tr, index) => {
        let cells = tr.cells;
        listOfSequences[index].five_prime_utr = Number(cells[8].textContent); // Update the five_prime_utr value
    });
    document.getElementById("loading").style.display = "";
    fetch("/optimize", {
        method: "POST",
        headers: {
            "Content-Type": "application/json"
        },
        body: JSON.stringify({
            targetI: targetI,
            iterations: iterations,
            method: method,
            sequences: listOfSequences
        })
    }).then(response => {
        return response.json()
    }).then(data => {
        for(i = 0; i < data.length; i++) {

            let outputBody = document.getElementById("optimized-output")
            outputBody.innerHTML = ""
            for (let i = 0; i < data.length; i++) {
                let k = data[i];
                let childCode = `
                <tr>
                <td>${k.tir}</td>
                <td>${k.I}</td>
                <td><input type="text" value="${k.gene}"></td>
            </tr>`
                let child = document.createElement("tr");
                child.innerHTML = childCode;
                outputBody.appendChild(child)
            }
        }
        document.getElementById("loading").style.display = "none"
    })
})