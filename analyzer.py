from genetic_code import codon_table

#valid dna sequence containing A, T, C and G only
def is_valid_dna(sequence):
    valid_nucleotides = {"A", "T", "G", "C"}
    if len(sequence) == 0:
        return False
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

#dna to codons to proteins
def translate_dna(dna):
    protein = ""
    for index in range(0, len(dna) - 2, 3):
        codon = dna[index:index + 3]
        amino_acid = codon_table[codon]
        if amino_acid == "*":
            protein = protein + amino_acid
            break
        protein = protein + amino_acid
    return protein

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

    reference_protein = translate_dna(reference)
    mutated_protein = translate_dna(mutated)

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
    print(f"Reference protein: {reference_protein}")
    print(f"Mutated protein: {mutated_protein}")

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
        if reference_amino_acid == mutated_amino_acid:
            classification = "synonymous"
        elif mutated_amino_acid == "*":
            classification = "nonsense"
        elif reference_amino_acid == "*":
            classification = "loss of stop codon"
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
        print(f"Amino acid position: {codon_number}")
        print(f"Amino acid change: {reference_amino_acid} -> {mutated_amino_acid}")
        print(f"Classification: {classification.upper()}")
        print("-" * 40)
