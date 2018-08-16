class ProteinProcessor:
    def __init__(self):
        self.codonTable = {"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
                           "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
                           "TAT": "Y", "TAC": "Y", "TAA": ".", "TAG": ".",
                           "TGT": "C", "TGC": "C", "TGA": ".", "TGG": "W",
                           "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
                           "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
                           "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
                           "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
                           "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
                           "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
                           "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
                           "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
                           "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
                           "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
                           "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
                           "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"}

        self.base_complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

    def translate_orf(self, orf):
        protein_sequence = ""
        for codon in range(0, (len(orf) - 1), 3):
            try:
                amino_acid = self.codonTable[orf[codon:codon + 3]]
                protein_sequence += amino_acid
            except KeyError:
                protein_sequence += 'X'
        return protein_sequence

    def _complement(self, s):
        letters = list(s)
        letters = [self.base_complement[base] for base in letters]
        return ''.join(letters)

    def reverse_complement(self, s):
        return self._complement(s[::-1])
