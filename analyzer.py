from Bio import Align
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

#alignment based indel function
def find_indel_with_alignment(reference, mutated):
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -3
    aligner.extend_gap_score = -0.5
    alignments = aligner.align(reference, mutated)

    alignment = alignments[0]
    coordinates = alignment.coordinates
    indels = []

    for index in range(len(coordinates[0]) - 1):
        reference_change = coordinates[0][index +1] - coordinates[0][index]
        mutated_change = coordinates[1][index +1] - coordinates[1][index]

        if reference_change > 0 and mutated_change == 0:
            position = int(coordinates[0][index] + 1)

            deleted_sequence = reference[
                coordinates[0][index] : coordinates[0][index + 1]
            ]
            indels.append(("deletion", position, deleted_sequence))

        elif reference_change == 0 and mutated_change > 0:
            position = int(coordinates[0][index] + 1)

            inserted_sequence = mutated[
                coordinates[1][index] : coordinates[1][index + 1]
            ]
            indels.append(("insertion", position, inserted_sequence))
    return indels

#creates an alignment and prepares an empty list to store substitutions
def find_substitutions_with_alignment(reference, mutated):
    aligner = Align.PairwiseAligner()

    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -3
    aligner.extend_gap_score = -0.5

    alignment = aligner.align(reference, mutated) [0]
    coordinates = alignment.coordinates
    substitutions = []

    for index in range(len(coordinates[0]) - 1):
        ref_start = int(coordinates[0][index])
        ref_end = int(coordinates[0][index + 1])

        mut_start = int(coordinates[1][index])
        mut_end = int(coordinates[1][index + 1])

        if ref_end - ref_start == mut_end - mut_start:
            for offset in range(ref_end - ref_start):
                ref_position = ref_start + offset
                mut_position = mut_start + offset

                if reference[ref_position] != mutated[mut_position]:
                    substitutions.append((
                        ref_position + 1,
                        reference[ref_position],
                        mutated[mut_position]
                    ))
    return substitutions

#combined variant detection
def find_all_variants(reference, mutated):
    aligner = Align.PairwiseAligner()

    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -3
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(reference, mutated)
    alignment = alignments[0]
    coordinates = alignment.coordinates

    substitutions = []
    indels = []

    for index in range(len(coordinates[0]) - 1):
        ref_start = int(coordinates[0][index])
        ref_end = int(coordinates[0][index + 1])

        mut_start = int(coordinates[1][index])
        mut_end = int(coordinates[1][index + 1])

        if ref_end > ref_start and mut_end == mut_start:
            deleted_sequence = reference[ref_start:ref_end]
            indels.append(("deletion", ref_start + 1, deleted_sequence))
        elif ref_end== ref_start and mut_end > mut_start:
            inserted_sequence = mutated[mut_start:mut_end]
            indels.append(("insertion", mut_start + 1, inserted_sequence))
        else:
            for offset in range(ref_end - ref_start):
                ref_position = ref_start + offset
                mut_position = mut_start + offset
                if reference[ref_position] != mutated[mut_position]:
                    substitutions.append((
                        ref_position + 1,
                        reference[ref_position],
                        mutated[mut_position]
                    ))
    return substitutions, indels

#substitution classification
def classify_substitution(reference_amino_acid, mutated_amino_acid):
    if reference_amino_acid == mutated_amino_acid:
        return "synonymous"
    elif mutated_amino_acid == ("*"):
        return "nonsense"
    elif reference_amino_acid == ("*"):
        return "loss of stop codon"
    else:
        return "missense"

#combined variant report
def analyze_all_variants(reference, mutated):
    substitutions, indels = find_all_variants(reference, mutated)

    reference_protein = translate_dna(reference)
    mutated_protein = translate_dna(mutated)
    protein_difference_index = find_first_difference(
        reference_protein,
        mutated_protein
    )

    print("\nCombined Variant Report")
    print("=" * 40)

    for position, ref_nucleotide, mut_nucleotide in substitutions:
        print(f"\nSubstitution at Position: {position}")
        print(f"Nucleotide change: {ref_nucleotide} -> {mut_nucleotide}")

        if not indels:
            codon_index = (position - 1) // 3
            codon_start = codon_index * 3
            reference_codon = reference[codon_start:codon_start + 3]
            mutated_codon = mutated[codon_start:codon_start + 3]

            if len(reference_codon) == 3 and len(mutated_codon) == 3:
                reference_amino_acid = codon_table[reference_codon]
                mutated_amino_acid = codon_table[mutated_codon]

                classification = classify_substitution(reference_amino_acid, mutated_amino_acid)

                print(f"Codon change: {reference_codon} -> {mutated_codon}")
                print(f"Amino acid position: {codon_index + 1}")
                print(f"Amino acid change: {reference_amino_acid} -> {mutated_amino_acid}")
                print(f"Classification: {classification.upper()}")

    for indel_type, position, sequence in indels:
        print(f"\n{indel_type.upper()} at position {position}")
        print(f"Nucleotide(s): {sequence}")

    print("\n Protein Analysis")
    print("=" * 40)
    print(f"Reference Protein: {reference_protein}")
    print(f"Mutated Protein: {mutated_protein}")

    if reference_protein == mutated_protein:
        print("Protein effect: No change detected")
    else:
        print(f"Protein divergence begins at amino acid: {protein_difference_index + 1}")

#analyzing indel
def analyze_indel(reference, mutated):
    indels = find_indel_with_alignment(reference, mutated)

    length_difference = abs(len(mutated) - len(reference))
    frame_effect = classify_length_change(reference, mutated)
    if len(indels) > 1 and len(reference) == len(mutated):
        frame_effect = "net reading frame restored"

    reference_protein = translate_dna(reference)
    mutated_protein = translate_dna(mutated)

    protein_difference_index = find_first_difference(reference_protein, mutated_protein)
    protein_difference_position = protein_difference_index + 1

    print("\nINDEL REPORT")
    print("=" * 40)
    for indel_type, position, changed_sequence in indels:
        print(f"\nType: {indel_type.upper()}")
        print(f"Position: {position}")
        print(f"Length: {len(changed_sequence)} nucleotide(s)")
        if len(changed_sequence) % 3 == 0:
            individual_effect = "in-frame"
        else:
            individual_effect = "frameshift"
        print(f"Individual reading frame effect: {individual_effect.upper()}")

        if indel_type == "insertion":
            print(f"Inserted nucleotide(s): {changed_sequence}")
        else:
            print(f"Deleted nucleotide(s): {changed_sequence}")

        print("-" * 40)
    print(f"Overall reading frame effect: {frame_effect.upper()}")
    print(f"Reference protein: {reference_protein}")
    print(f"Mutated protein: {mutated_protein}")
    print(f"Protein divergence begins at amino acid: {protein_difference_position}")

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
