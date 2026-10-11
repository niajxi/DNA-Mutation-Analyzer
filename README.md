# DNA Mutation Analyzer
#### A Python-based computational biology tool for detecting DNA sequence variants and analyzing their potential effects on the corresponding protein sequences 
## Overview
The DNA Mutation Analyzer compares a reference DNA sequence with a mutated sequence to identify nucleotide substitutions, insertions and deletions.
It uses pairwise sequence alignment to locate variants and evaluates their effects on codons, amino acid sequences and reading frames.

## Features
* DNA sequence validation
* Detection of nucleotide substitutions, insertions and deletions
* Pairwise DNA sequence alignment using Biopython
* Classification of synonymous, missense, nonsense and stop-loss substitutions
* Identification of in-frame insertions and deletions
* Detection of frameshifts and reading-frame restoration
* DNA-to-protein translation 
* Identification of the first amino acid position where protein sequences diverge
* Alignment-based codon validation for sequences containing multiple variants
* Automated testing using pytest

## Installation
To run the DNA Mutation Analyzer locally, Python 3.13 or a compatible version must be installed.

### 1. Clone the repository

Clone the repository and navigate to the project directory:
```bash
git clone https://github.com/niajxi/DNA-Mutation-Analyzer.git
cd DNA-Mutation-Analyzer
```

### 2. Create a virtual environment
```bash
python3 -m venv .venv
```

Activate the environment:

**macOS/Linux:**
```bash
source .venv/bin/activate
```

**Windows (PowerShell):**
```powershell
.venv\Scripts\Activate.ps1
```

### 3. Install dependencies
```bash
python -m pip install -r requirements.txt
```

### 4. Run the program
```bash
python main.py
```

The program will prompt you to enter a reference DNA sequence and a mutated DNA sequence.

### 5. Run automated tests (optional)

```bash
python -m pytest testing_analyzer.py -v
```
## Usage
The user enters a reference DNA sequence and a mutated DNA sequence. The program generates a report identifying 
detected variants, their classifications, and their potential effects on protein sequences.


### Example Input and Output

**Reference DNA:**
```text
ATGGCTAAACCG
```

**Mutated DNA:**
```text
ATGACTAAACCG
```

**Example Output:**

```text
Combined Variant Report
========================================

Substitution at Position: 4
Nucleotide change: G -> A
Codon change: GCT -> ACT
Amino acid position: 2
Amino acid change: A -> T
Classification: MISSENSE

Protein Analysis
========================================
Reference Protein: MAKP
Mutated Protein: MTKP
Protein divergence begins at amino acid: 2
```

In this example, a nucleotide substitution at position 4 changes the codon from GCT to ACT, resulting in a missense
mutation that replaces alanine (A) with threonine (T) at amino acid position 2.

## Limitations
* The analyzer assumes that DNA sequences are provided in the correct reading frame, beginning at the first nucleotide
* Variant positions may be ambiguous in repetitive DNA sequences because multiple alignments can be equally valid
* The program uses the standard genetic code and does not account for alternative splicing or different genetic codes
* Codon-level classifications are skipped when insertions or deletions disrupt the alignment of the corresponding codons
* The predicted protein effects are based on sequence changes and do not determine whether a variant is pathogenic or not

## Future Improvements
* **FASTA file support:** allow users to upload DNA sequences from FASTA files instead of entering them manually.
* **CSV export:** Enable users to save variant reports as CSV files for further analysis
* **Improved alignment handling:** Support alternative alignments when variant positions are ambiguous
* **Genomic annotations:** Incorporate reference gene annotations to provide more detailed information about variant
locations and their potential consequences
* **Larger sequence support:** Optimize the analyzer for longer DNA sequences
