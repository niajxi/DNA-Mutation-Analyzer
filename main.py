from analyzer import is_valid_dna, analyze_indel, analyze_mutations, find_indel_with_alignment

if __name__ == "__main__":
    reference_dna = input("Enter reference DNA sequence: ").upper()
    mutated_dna = input("Enter mutated DNA sequence: ").upper()

    if len(reference_dna) == 0:
        print("Error: Reference sequence cannot be empty.")
    elif len(mutated_dna) == 0:
        print("Error: Mutated sequence cannot be empty.")
    elif not is_valid_dna(reference_dna):
        print("Error: Reference sequence contains invalid nucleotides.")
    elif not is_valid_dna(mutated_dna):
        print("Error: Mutated sequence contains invalid nucleotides.")


    else:
        indels = find_indel_with_alignment(reference_dna, mutated_dna)
        if indels:
            analyze_indel(reference_dna, mutated_dna)
        elif len(reference_dna) != len(mutated_dna):
            print("Error: Sequence lengths differ, but no indels were detected.")
        elif len(reference_dna) % 3 != 0:
            print("Error: DNA sequence length must be divisible by 3.")
        else:
            analyze_mutations(reference_dna, mutated_dna)