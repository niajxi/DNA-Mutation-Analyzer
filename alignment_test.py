from Bio import Align

aligner = Align.PairwiseAligner()

reference = "ATGAAAAGCT"
mutated = "ATGAAAGCT"

alignments = aligner.align(reference, mutated)

print(alignments[0])

alignments = alignments[0]
print(alignments.coordinates)

coordinates = alignments.coordinates

for index in range(len(coordinates[0]) - 1):
    reference_change = coordinates[0][index + 1] - coordinates[0][index]
    mutated_change = coordinates[1][index + 1] - coordinates[1][index]

    if reference_change > 0 and mutated_change == 0:
        position = coordinates[0][index] + 1

        deleted_sequence = reference[
            coordinates[0][index]:coordinates[0][index + 1]
        ]
        print("Deletion Detected")
        print(f"Position: {position}")
        print(f"Deleted nucleotide(s): {deleted_sequence}")

    elif reference_change == 0 and mutated_change > 0:
        position = coordinates[0][index] + 1

        inserted_sequence = mutated[
            coordinates[1][index]:coordinates[1][index + 1]
        ]
        print("Insertion Detected")
        print(f"Position: {position}")
        print(f"Inserted nucleotide(s): {inserted_sequence}")


