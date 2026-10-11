from Bio import Align
from genetic_code import codon_table

#valid dna sequence containing A, T, C and G only
def is_valid_dna(sequence):
    valid_nucleotides = {"A", "T", "G", "C"}
    if len(sequence) == 0:
        return False
    return set(sequence).issubset(valid_nucleotides)

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

#finding the first mismatch to see where indel begins
def find_first_difference(reference, mutated):
    shortest_length = min(len(reference), len(mutated))

    for index in range(shortest_length):
        if reference[index] != mutated[index]:
            return index

    return shortest_length

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
            indels.append(("insertion", ref_start + 1, inserted_sequence))
        else:
            for offset in range(ref_end - ref_start):
                ref_position = ref_start + offset
                mut_position = mut_start + offset
                if reference[ref_position] != mutated[mut_position]:
                    substitutions.append((
                        ref_position + 1,
                        mut_position + 1,
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

#locating codons
def get_codon_at_position(dna, position):
    codon_index = (position - 1) // 3
    codon_start = codon_index * 3
    codon = dna[codon_start:codon_start + 3]

    if len(codon) != 3:
        return None
    return codon

#alignment-based codon validation
def codons_align_correctly(reference, mutated, position, mutated_position):
    reference_start = ((position - 1) // 3) * 3
    mutated_start = ((mutated_position - 1) // 3) * 3

    if reference_start + 3 > len(reference):
        return False

    if mutated_start + 3 > len(mutated):
        return False

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -3
    aligner.extend_gap_score = -0.5

    alignment = aligner.align(reference, mutated)[0]
    aligned_positions = {}

    for (ref_start, ref_end), (mut_start, mut_end) in zip(
            alignment.aligned[0],
            alignment.aligned[1]
    ):
        for offset in range(ref_end - ref_start):
            aligned_positions[ref_start + offset] = mut_start + offset

    for offset in range(3):
        if aligned_positions.get(reference_start + offset) != mutated_start + offset:
            return False
    return True

#combined variant report
def analyze_all_variants(reference, mutated):
    substitutions, indels = find_all_variants(reference, mutated)

    if not substitutions and not indels:
        print("\nNo DNA variants detected.")
        return

    reference_protein = translate_dna(reference)
    mutated_protein = translate_dna(mutated)
    protein_difference_index = find_first_difference(
        reference_protein,
        mutated_protein
    )

    print("\nCombined Variant Report")
    print("=" * 40)

    reported_codons = set()

    for position, mutated_position, ref_nucleotide, mut_nucleotide in substitutions:
        print(f"\nSubstitution at Position: {position}")
        print(f"Nucleotide change: {ref_nucleotide} -> {mut_nucleotide}")

        if codons_align_correctly(reference, mutated, position, mutated_position):
            codon_index = (position - 1) // 3
            if codon_index in reported_codons:
                continue
            reported_codons.add(codon_index)
            reference_codon = get_codon_at_position(reference, position)
            mutated_codon = get_codon_at_position(mutated, mutated_position)

            if reference_codon is not None and mutated_codon is not None:
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

        if len(sequence) % 3 == 0:
            classification = "IN FRAME"
        else:
            classification = "FRAMESHIFT"
        print(f"Classification: {classification}")
    if indels:
        length_difference = len(mutated) - len(reference)
        print("\nOverall Reading Frame Analysis:")
        print("=" * 40)

        if length_difference % 3 == 0:
            if any(len(sequence) % 3 !=0 for _, _, sequence in indels):
                print("Reading frame restored after multiple indels")
            else:
                print("Reading frame preserved")
        else:
            print("Net frameshift detected")

    print("\nProtein Analysis")
    print("=" * 40)
    print(f"Reference Protein: {reference_protein}")
    print(f"Mutated Protein: {mutated_protein}")

    if reference_protein == mutated_protein:
        print("Protein effect: No change detected")
    else:
        print(f"Protein divergence begins at amino acid: {protein_difference_index + 1}")

