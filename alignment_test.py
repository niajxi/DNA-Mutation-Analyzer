from Bio import Align

from analyzer import find_indel_with_alignment


reference = "ATGAAAGCT"
mutated = "ATGAAAAGCT"

result = find_indel_with_alignment(reference, mutated)

print(result)
