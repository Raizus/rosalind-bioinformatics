

from collections import Counter
import random
import re
from typing import Iterator, Literal, Union

from BioInfoToolkit.Sequences.StringUtils import findMotif
from .SequenceUtils import gc_content, generate_random_sequence
from .structures import DNA_CODONS, NUCLEOTIDE_BASE, RNA_CODONS


class BioSequence:
    """DNA Sequence class"""
    label: str
    seq: str
    seq_type: str
    is_valid: bool

    def __init__(self, seq: str, 
                 seq_type: str = "DNA", 
                 label: str = "no label") -> None:
        self.seq = seq.upper()
        self.label = label
        self.seq_type = seq_type
        self.is_valid = self.__validate()
        assert self.is_valid, f"Provided data is not a valid sequence {self.seq}"

    def __validate(self) -> bool:
        """Check if the sequence is a valid DNA sequence string

        Args:
            dna_seq (str): _description_

        Returns:
            bool | str: Uppercase version of string or false
        """

        return set(NUCLEOTIDE_BASE[self.seq_type]).issuperset(self.seq)

    def get_seq_biotype(self):
        return self.seq_type

    def length(self):
        return len(self.seq)

    def get_seq_info(self) -> str:
        return f"[Label]: {self.label}\n[Sequence]: {self.seq}\n[Byotype]: {self.seq_type}\n[Length]: {len(self.seq)}"

    @classmethod
    def generate_random_seq(cls, length: int, seq_type="DNA", label: str = "no label") -> "BioSequence":
        seq = generate_random_sequence(length, NUCLEOTIDE_BASE[seq_type])
        return BioSequence(seq, seq_type, label)

    def nucleotide_frequency(self) -> dict[str, int]:
        """Counts the frequency of each nucleotide

        Returns:
            dict[str, int]: counts
        """
        c = dict(Counter(self.seq))
        # for nuc in DNA_NUCLEOTIDES:
        #     c.setdefault(nuc, 0)
        return c

    def transcription(self) -> str:
        """Translates DNA to RNA by substituting 'T' by 'U'

        Returns:
            str: rna_seq
        """
        if self.seq_type == "DNA":
            return self.seq.replace("T", "U")
        raise TypeError(
            "Sequence is not a DNA sequence, cannot be transcribed")

    def complement(self) -> str:
        if self.seq_type == "DNA":
            mapping = str.maketrans('ATCG', 'TAGC')
        else:
            mapping = str.maketrans('AUCG', 'UAGC')
        return self.seq.translate(mapping)

    def reverseComplement(self) -> str:
        return self.complement()[::-1]

    def gc_content(self) -> float:
        """Returns GC constent in a DNA / RNA sequence in % value

        Args:
            seq (str): _description_

        Returns:
            _type_: _description_
        """
        return gc_content(self.seq)

    def gc_content_subseq(self, k: int = 20) -> list[float]:
        """Returns GC constent in a DNA / RNA sequence in % value

        Args:
            seq (str): _description_

        Returns:
            _type_: _description_
        """
        res = []
        for i in range(0, len(self.seq) - k + 1, k):
            subseq = self.seq[i:i+k]
            res.append(gc_content(subseq))
        return res

    def translate_seq(self, init_pos: int = 0) -> list[str]:
        """Translates a DNA or RNA string into a aminoacid sequence string

        Args:
            init_pos (int, optional): Position to start translation. Defaults to 0.

        Returns:
            list[str]: aminoacid sequence
        """
        if self.seq_type == "DNA":
            return [DNA_CODONS[self.seq[p:p+3]] for p in range(init_pos, len(self.seq) - 2, 3)]
        else:
            return [RNA_CODONS[self.seq[p:p+3]] for p in range(init_pos, len(self.seq) - 2, 3)]

    def codon_usage(self, aminoacid: str) -> dict[str, float]:
        if self.seq_type == "DNA":
            tmpList = [self.seq[i:i+3]
                       for i in range(0, len(self.seq) - 2, 3)
                       if DNA_CODONS[self.seq[i:i+3]] == aminoacid]
        else:
            tmpList = [self.seq[i:i+3]
                       for i in range(0, len(self.seq) - 2, 3)
                       if RNA_CODONS[self.seq[i:i+3]] == aminoacid]

        freqDict: dict[str, float] = dict(Counter(tmpList))
        totalWeight = sum(freqDict.values())
        for seq in freqDict:
            freqDict[seq] = round(
                freqDict[seq] / totalWeight, 2)  # type: ignore
        return freqDict

    def reading_frames_iterator(self, complement: bool = False):
        if not complement:
            for i in range(3):
                yield [self.seq[p:p+3] for p in range(i, len(self.seq) - 2, 3)]
        else:
            rComp = BioSequence(self.reverseComplement(), self.seq_type)
            for i in range(3):
                yield [rComp.seq[p:p+3] for p in range(i, len(rComp.seq) - 2, 3)]

    def reading_frames_iterator_aminoacids(self) -> Iterator[list[str]]:
        """Generate the six reading frames of a DNA sequence, including the reverse complement

        Args:
            seq (str): _description_
        """
        yield self.translate_seq(0)
        yield self.translate_seq(1)
        yield self.translate_seq(2)
        reverseComplement = BioSequence(self.reverseComplement(), self.seq_type)
        yield reverseComplement.translate_seq(0)
        yield reverseComplement.translate_seq(1)
        yield reverseComplement.translate_seq(2)
        del reverseComplement

    def reading_frames_list_aminoacids(self) -> list[list[str]]:
        """Generate the six reading frames of a DNA sequence, including the reverse complement

        Args:
            seq (str): _description_
        """
        return [frame for frame in self.reading_frames_iterator_aminoacids()]

    @classmethod
    def proteins_from_rf(cls, aa_seq: list[str]) -> list[str]:
        """Finds proteins in a reading frame (aminoacid sequence), if any

        Args:
            aa_seq (list[str]): _description_

        Returns:
            list[str]: protein list
        """
        current_prot: list[str] = []
        proteins: list[str] = []
        for aa in aa_seq:
            if aa == '_':
                # stop codon
                if not current_prot:
                    continue
                for p in current_prot:
                    proteins.append(p)
                current_prot = []
            else:
                if aa == 'M':
                    current_prot.append("")
                for i in range(len(current_prot)):
                    current_prot[i] += aa
        return proteins

    def all_proteins_from_rfs(self, startReadPos: int = 0, endReadPos: int = 0, ordered: bool = False) -> list[str]:
        """Returns all possible DISTINCT proteins from all reading frames

        Args:
            startReadPos (int, optional): _description_. Defaults to 0.
            endReadPos (int, optional): _description_. Defaults to 0.
            ordered (bool, optional): _description_. Defaults to False.

        Returns:
            list[str]: _description_
        """
        if endReadPos > startReadPos:
            tmp_seq = BioSequence(self.seq[startReadPos:endReadPos], self.seq_type)
            rfs_iter = tmp_seq.reading_frames_iterator_aminoacids()
        else:
            rfs_iter = self.reading_frames_iterator_aminoacids()

        res: list[str] = []
        for rf in rfs_iter:
            prots = BioSequence.proteins_from_rf(rf)
            for p in prots:
                if p not in res:
                    res.append(p)

        if ordered:
            return sorted(res, key=len, reverse=True)
        return res

    def find_substring(self, substr: str) -> list[int]:
        """Returns the starting indexes of the substring

        Args:
            substr (str): _description_

        Returns:
            list[int]: _description_
        """
        idxs: list[int] = []
        idx = 0
        while True:
            idx = self.seq.find(substr, idx)
            if idx == -1:
                return idxs
            idxs.append(idx)
            idx += 1

    def find_reverse_palindrome(self, l: int) -> list[int]:
        """Finds reverse palindromes with lenght l in the sequence

        Args:
            l (int): palindrome length

        Returns:
            list[int]: indexes (0-indexed) of the start of each palindromes
        """
        compl = self.complement()
        positions = []
        for i in range(0, self.length()-l + 1):
            # print(f"{i}: {self.seq[i:i+l]} == {compl[i:i+l][::-1]}")
            if self.seq[i:i+l] == compl[i:i+l][::-1]:
                positions.append(i)
        return positions

    def RNASplicing(self, introns: list[str]) -> str:
        """Splices the given introns out of the current sequence,
        and translates the resulting sequence into a protein

        Args:
            introns (list[str]): list of introns

        Returns:
            str: The resulting protein
        """
        intron_pattern = '|'.join(introns)
        intron_pattern_re = re.compile(f"({intron_pattern})")

        splice = intron_pattern_re.sub("", self.seq)
        splice_seq = BioSequence(splice, self.seq_type)

        aas = splice_seq.translate_seq()
        proteins = sorted(BioSequence.proteins_from_rf(
            aas), key=lambda p: len(p), reverse=True)
        return proteins[0]

    def peptideDecoding(self, peptide: str) -> list[str]:
        """Find substrings of a sequence encoding a given amino acid sequence.

        We say that a DNA string Pattern encodes an amino acid string Peptide if the RNA string transcribed from either Pattern or its reverse complement Pattern translates into Peptide.

        Args:
            peptide (str): _description_

        Returns:
            list[str]: list of seubstrings of sequence encoding the peptide
        """
        substrings: list[str] = []
        l = len(peptide)
        for i in range(3):
            aa_seq = ''.join(self.translate_seq(i))
            matches = findMotif(peptide, aa_seq)
            for pos, aa_match in matches.items():
                ss = self.seq[i+pos*3:i+(pos+l)*3]
                substrings.append(ss)

        rComp = BioSequence(self.reverseComplement(), self.seq_type)
        for i in range(3):
            aa_seq = ''.join(rComp.translate_seq(i))
            matches = findMotif(peptide, aa_seq)
            for pos, aa_match in matches.items():
                ss = rComp.seq[i+pos*3:i+(pos+l)*3]
                ss = BioSequence(ss, self.seq_type).reverseComplement()
                substrings.append(ss)
                
        return substrings


