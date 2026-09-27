#variable creation
reference_dna = input("Enter reference DNA sequence: ").upper()
mutated_dna = input("Enter mutated DNA sequence: ").upper()

#valid dna sequence containing A, T, C and G only
def is_valid_dna(sequence):
    valid_nucleotides = {"A", "T", "G", "C"}
    return set(sequence).issubset(valid_nucleotides)

#finding mutations
def find_mutations(reference, mutated):
    mutations = []
    for index, nucleotide in enumerate(reference):
        mutated_nucleotide = mutated[index]

        if nucleotide != mutated_nucleotide:
            position = index + 1
            mutations.append((position, nucleotide, mutated_nucleotide))
    return mutations

#standard codons
codon_table = {
    # Phenylalanine / Leucine
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    # Leucine
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    # Isoleucine / Methionine
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    # Valine
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    # Serine
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    # Proline
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    # Threonine
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    # Alanine
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    # Tyrosine / Stop
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    # Histidine / Glutamine
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    # Asparagine / Lysine
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    # Aspartic acid / Glutamic acid
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    # Cysteine / Stop / Tryptophan
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    # Arginine
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    # Serine / Arginine
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    # Glycine
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}

#dna to codons to proteins
def translate_dna(dna):
    protein = ""
    for index in range(0, len(dna), 3):
        codon = dna[index:index + 3]
        amino_acid = codon_table[codon]
        protein = protein + amino_acid
    return protein

#classify types of mutations
def compare_proteins(reference, mutated):
    for index, amino_acid in enumerate(reference):
        mutated_amino_acid = mutated[index]

        if amino_acid != mutated_amino_acid:
            position = index + 1
            print(f"Protein change at position {position}: {amino_acid} -> {mutated_amino_acid}")


#length change
def classify_length_change(reference, mutated):
    length_difference = len(mutated) - len(reference)

    if length_difference == 0:
        return "same length"
    elif length_difference % 3 == 0:
        return "in-frame indel"
    else:
        return "frameshift"


#finding the first mismatch to see where indel begins
def find_first_difference(reference, mutated):
    shortest_length = min(len(reference), len(mutated))

    for index in range(shortest_length):
        if reference[index] != mutated[index]:
            return index

    return shortest_length

#distinguishing insertion vs deletion
def classify_indel(reference, mutated):
    if len(mutated) > len(reference):
        return "insertion"
    elif len(mutated) < len(reference):
        return "deletion"
    else:
        return "no indel"

#analyzing indel
def analyze_indel(reference, mutated):
    indel_type = classify_indel(reference, mutated)
    indel_index = find_first_difference(reference, mutated)
    position = indel_index + 1

    length_difference = abs(len(mutated) - len(reference))
    frame_effect = classify_length_change(reference, mutated)

    if indel_type == "insertion":
        changed_sequence = mutated[indel_index:indel_index + length_difference]
    else:
        changed_sequence = reference[indel_index:indel_index + length_difference]

    print("\nINDEL REPORT")
    print("=" * 40)
    print(f"Type: {indel_type.upper()}")
    print(f"Position: {position}")
    print(f"Length change: {length_difference} nucleotide(s)")
    print(f"Reading frame effect: {frame_effect.upper()}")
    if indel_type == "insertion":
        print(f"Inserted nucleotide(s): {changed_sequence}")
    else:
        print(f"Deleted nucleotide(s): {changed_sequence}")


# analyze mutations
def analyze_mutations(reference, mutated):
    mutations = find_mutations(reference, mutated)

    # if the sequences are identical, there are no mutations to analyze
    if len(mutations) == 0:
        print("\nNo mutations detected.")
        return

    print("\nMUTATION REPORT")
    print("=" * 40)

    reported_codons = set()

    for position, reference_nucleotide, mutated_nucleotide in mutations:
        index = position - 1
        codon_number = (index // 3) + 1
        codon_start = (index // 3) * 3

        # skip this codon if we already analyzed it
        if codon_number in reported_codons:
            continue

        reported_codons.add(codon_number)

        # get the reference and mutated codons
        reference_codon = reference[codon_start:codon_start + 3]
        mutated_codon = mutated[codon_start:codon_start + 3]

        # translate the codons
        reference_amino_acid = codon_table[reference_codon]
        mutated_amino_acid = codon_table[mutated_codon]

        # classify the protein consequence
        if mutated_amino_acid == "*":
            classification = "nonsense"
        elif reference_amino_acid == mutated_amino_acid:
            classification = "synonymous"
        else:
            classification = "missense"

        # find ALL nucleotide mutations belonging to this codon
        codon_mutations = []

        for mutation_position, ref_nucleotide, mut_nucleotide in mutations:
            mutation_index = mutation_position - 1
            mutation_codon_number = (mutation_index // 3) + 1

            if mutation_codon_number == codon_number:
                codon_mutations.append(
                    (mutation_position, ref_nucleotide, mut_nucleotide)
                )

        # make readable lists for the report
        positions = []
        nucleotide_changes = []

        for mutation_position, ref_nucleotide, mut_nucleotide in codon_mutations:
            positions.append(str(mutation_position))
            nucleotide_changes.append(
                f"{ref_nucleotide} -> {mut_nucleotide}"
            )

        # print the mutation report
        print(f"\nDNA position(s): {', '.join(positions)}")
        print(f"Nucleotide change(s): {', '.join(nucleotide_changes)}")
        print(f"Codon change: {reference_codon} -> {mutated_codon}")
        print(f"Amino acid change: {reference_amino_acid} -> {mutated_amino_acid}")
        print(f"Classification: {classification.upper()}")
        print("-" * 40)

if not is_valid_dna(reference_dna):
    print("Error: Reference sequence contains invalid nucleotides.")
elif not is_valid_dna(mutated_dna):
    print("Error: Mutated sequence contains invalid nucleotides.")
elif len(reference_dna) != len(mutated_dna):
    analyze_indel(reference_dna, mutated_dna)
elif len(reference_dna) % 3 != 0:
    print("Error: DNA sequence length must be divisible by 3.")
else:
    analyze_mutations(reference_dna, mutated_dna)