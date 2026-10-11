from analyzer import (is_valid_dna, translate_dna, classify_substitution, find_all_variants, codons_align_correctly,
                      analyze_all_variants)

def test_dna_validation():
    assert is_valid_dna("ATGC") is True
    assert is_valid_dna("ATGX") is False
    assert is_valid_dna("") is False

def test_dna_translation():
    assert translate_dna("ATGGCTAAACCG") == "MAKP"
    assert translate_dna("ATGTAA") == "M*"

def test_substitution_classification():
    assert classify_substitution("A", "A") == "synonymous"
    assert classify_substitution("A", "T") == "missense"
    assert classify_substitution("Q", "*") == "nonsense"
    assert classify_substitution("*", "Q") == "loss of stop codon"

def test_substitution_detection():
    reference = "ATGGCTAAACCG"
    mutated = "ATGACTAAACCG"
    substitutions, indels = find_all_variants(reference, mutated)
    assert substitutions == [(4, 4, "G", "A")]
    assert indels == []

def test_insertion_detection():
    reference = "ATGAAAGCTCCG"
    mutated = "ATGAAAGGGGCTCCG"
    substitutions, indels = find_all_variants(reference, mutated)
    assert substitutions == []
    assert indels == [("insertion", 7, "GGG")]

def test_deletion_detection():
    reference = "ATGAAATTTGCTCCG"
    mutated = "ATGAAAGCTCCG"
    substitutions, indels = find_all_variants(reference, mutated)
    assert substitutions == []
    assert indels == [("deletion", 7, "TTT")]

def test_codons_align_after_inframe_insertion():
    reference = "ATGAAAGCTCCG"
    mutated = "ATGAAATTTGCTCAG"
    assert codons_align_correctly(reference, mutated, 11, 14) is True

def test_codons_disrupted_by_frameshift():
    reference = "ATGGCTAAACCG"
    mutated = "ATGACTTAAACCG"
    assert codons_align_correctly(reference, mutated, 4, 4) is False

def test_reading_frame_restoration():
    reference = "ATGAAAGCTTTTCCCGGG"
    mutated = "ATGAAATGCTTTTCCGGG"
    substitutions, indels = find_all_variants(reference, mutated)
    assert substitutions == []
    assert indels == [
        ("insertion", 7, "T"),
        ("deletion", 13, "C")
    ]
    assert (len(mutated) - len(reference)) % 3 == 0

def test_complete_variant_report(capsys):
    reference = "ATGGCTAAACCG"
    mutated = "ATGACTAAACCG"
    analyze_all_variants(reference, mutated)
    output = capsys.readouterr().out
    assert "Combined Variant Report" in output
    assert "Substitution at Position: 4" in output
    assert "Codon change: GCT -> ACT" in output
    assert "Classification: MISSENSE" in output
    assert "Reference Protein: MAKP" in output
    assert "Mutated Protein: MTKP" in output

def test_report_skips_disrupted_codon(capsys):
    reference = "ATGGCTAAACCG"
    mutated = "ATGACTTAAACCG"

    analyze_all_variants(reference, mutated)

    output = capsys.readouterr().out

    assert "Substitution at Position: 4" in output
    assert "Classification: FRAMESHIFT" in output
    assert "Net frameshift detected" in output

    assert "Codon change:" not in output
    assert "Classification: MISSENSE" not in output

    assert "Reference Protein: MAKP" in output
    assert "Mutated Protein: MT*" in output

def test_multiple_substitutions_same_codon(capsys):
    reference = "ATGGCTAAACCG"
    mutated = "ATGACAAAACCG"
    analyze_all_variants(reference, mutated)
    output = capsys.readouterr().out

    assert "Substitution at Position: 4" in output
    assert "Substitution at Position: 6" in output
    assert "Codon change: GCT -> ACA" in output
    assert "Amino acid change: A -> T" in output
    assert output.count("Codon change: GCT -> ACA") == 1
    assert output.count("Substitution at Position:") == 2