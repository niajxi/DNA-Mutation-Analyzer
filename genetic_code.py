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