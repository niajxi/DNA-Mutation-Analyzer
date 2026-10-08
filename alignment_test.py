from Bio import Align

from analyzer import find_indel_with_alignment

reference = "ATGGCTAAACCG"
mutated = "ATGACTAAACCG"

indels = find_indel_with_alignment(reference, mutated)

for indel_type, position, sequence in indels:
    print(f"{indel_type.upper()} at position {position}: {sequence}")
