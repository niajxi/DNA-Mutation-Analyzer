from Bio import Align

from analyzer import analyze_all_variants

reference = "ATGAAAGCTTTTCCCGGG"
mutated = "ATGAAATGCTTTTCCGGA"

analyze_all_variants(reference, mutated)