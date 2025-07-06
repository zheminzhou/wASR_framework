import sys
from os.path import isfile
from collections import defaultdict
import numpy as np
from Bio import AlignIO, SeqIO
import config as ttconf
from seq_utils import seq2array, guess_alphabet, alphabets

# Compatibility with Python 2/3
string_types = [str] if sys.version_info[0] == 3 else [str, unicode]

def simple_logger(*args, **kwargs):
    print(args)

class SequenceData(object):
    def __init__(self, aln, ref=None, logger=None, convert_upper=True,
                 sequence_length=None, compress=True, word_length=1, sequence_type=None,
                 fill_overhangs=True, seq_multiplicity=None, ambiguous=None, **kwargs):
        self.logger = logger or simple_logger
        self._aln = None
        self._ref = None
        self.likely_alphabet = None
        self.compressed_to_full_sequence_map = None
        self._multiplicity = None
        self.is_sparse = None
        self.convert_upper = convert_upper
        self.compress = compress
        self.seq_multiplicity = seq_multiplicity or {}
        self.additional_constant_sites = kwargs.get('additional_constant_sites', 0)

        self._full_length = None
        self.full_length = sequence_length
        self._compressed_length = None
        self.word_length = word_length
        self.fill_overhangs = fill_overhangs
        self.ambiguous = ambiguous
        self.sequence_type = sequence_type

        self.ref = ref
        self.aln = aln

    @property
    def aln(self):
        return self._aln

    @aln.setter
    def aln(self, in_aln):
        from Bio.Align import MultipleSeqAlignment
        self._aln, self.is_sparse = None, None
        if in_aln is None:
            return

        if isinstance(in_aln, (dict, defaultdict)):
            self._aln = in_aln
            self.is_sparse = True
        elif isinstance(in_aln, string_types) and isfile(in_aln):
            for fmt in ['fasta', 'phylip-relaxed', 'nexus']:
                try:
                    in_aln = AlignIO.read(in_aln, fmt)
                    break
                except Exception:
                    continue

        if isinstance(in_aln, MultipleSeqAlignment):
            self._aln = {}
            for s in in_aln:
                tmp_name = (
                    s.id if s.id == s.name else
                    s.name if '<unknown' in s.id else
                    s.id if '<unknown' in s.name else
                    s.name
                )
                self._aln[tmp_name] = seq2array(
                    s,
                    convert_upper=self.convert_upper,
                    fill_overhangs=self.fill_overhangs,
                    ambiguous=self.ambiguous
                )
            self.check_alphabet(list(self._aln.values()))
            self.is_sparse = False
            self.logger("SequenceData: loaded alignment.", 1)

        elif isinstance(in_aln, (dict, defaultdict)):
            self.logger("SequenceData: loaded sparse/vcf alignment.", 1)
            self.check_alphabet([self.ref])
            self.is_sparse = True
            self._aln = in_aln
        else:
            raise TypeError(f"SequenceData: loading alignment failed... {in_aln}")

        if self.full_length:
            if self.is_sparse:
                if self.full_length != len(self.ref):
                    self.logger("SequenceData.aln: specified sequence length doesn't match reference length, ignoring sequence length.", 1, warn=True)
                    self._full_length = len(self.ref)
            else:
                if self.full_length < in_aln.get_alignment_length():
                    raise AttributeError("SequenceData.aln: specified sequence length is smaller than alignment length!")
                elif self.full_length > in_aln.get_alignment_length():
                    self.logger("SequenceData.aln: specified sequence length doesn't match alignment length. Treating difference as constant sites.", 2, warn=True)
                    self.additional_constant_sites = max(0, self.full_length - in_aln.get_alignment_length())
        else:
            self.full_length = len(self.ref) if self.is_sparse else in_aln.get_alignment_length()

        self.sequence_names = list(self.aln.keys())
        self.make_compressed_alignment()

    @property
    def full_length(self):
        return self._full_length

    @full_length.setter
    def full_length(self, L):
        if getattr(self, '_full_length', None) is None and L:
            self._full_length = int(L)
        elif L:
            self.logger("Alignment: one_mutation and sequence length can only be specified once!", 1)

    @property
    def compressed_length(self):
        return self._compressed_length

    @property
    def ref(self):
        return self._ref

    @ref.setter
    def ref(self, in_ref):
        if in_ref and isfile(in_ref):
            for fmt in ['fasta', 'genbank']:
                try:
                    in_ref = SeqIO.read(in_ref, fmt)
                    self.logger(f"SequenceData: loaded reference sequence as {fmt} format", 1)
                    break
                except Exception:
                    continue
            else:
                raise TypeError("SequenceData.ref: reference sequence file could not be parsed. Only FASTA and GenBank formats are supported.")

        if in_ref:
            self._ref = seq2array(in_ref, fill_overhangs=False, word_length=self.word_length)
            self.full_length = self._ref.shape[0]
            self.compressed_to_full_sequence_map = None
            self._multiplicity = None

    def multiplicity(self, weight=None):
        return self._multiplicity if weight is None else self._multiplicity * weight

    def check_alphabet(self, seqs):
        self.likely_alphabet = guess_alphabet(seqs)

        if self.sequence_type and self.likely_alphabet != self.sequence_type:
            warning = {
                'nuc': "POSSIBLE ERROR: This does not look like a nucleotide alignment!",
                'aa': "POSSIBLE ERROR: This looks like a nucleotide alignment, you indicated amino acids!"
            }
            self.logger(warning.get(self.sequence_type, ""), 0, warn=True)

        if self.ambiguous is None:
            self.ambiguous = 'N' if self.likely_alphabet == 'nuc' else 'X'

    def make_compressed_alignment(self):
        if not self.compress:
            self._multiplicity = np.ones(self.full_length, dtype=float)
            self.full_to_compressed_sequence_map = np.arange(self.full_length)
            self.compressed_to_full_sequence_map = {p: np.array([p]) for p in np.arange(self.full_length)}
            self._compressed_length = self._full_length
            self.compressed_alignment = self._aln
            return ttconf.SUCCESS

        self.logger("SeqData: making compressed alignment...", 1)
        self.full_to_compressed_sequence_map = np.zeros(self.full_length, dtype=int)
        self.compressed_to_full_sequence_map = {}

        if self.is_sparse:
            from seq_utils import process_sparse_alignment
            tmp = process_sparse_alignment(self.aln, self.ref, self.ambiguous)
            compressed_aln_transpose = tmp["constant_columns"]
            alignment_patterns = tmp["constant_patterns"]
            variable_positions = tmp["variable_positions"]
            self.inferred_const_sites = tmp["constant_up_to_ambiguous"]
            self.nonref_positions = tmp["nonref_positions"]
        else:
            alignment_patterns = {}
            compressed_aln_transpose = []
            aln_transpose = np.array([self.aln[k] for k in self.sequence_names]).T
            variable_positions = np.arange(aln_transpose.shape[0])

        for pi in variable_positions:
            pattern = (
                np.array([self.aln[k][pi] if pi in self.aln[k] else self.ref[pi] for k in self.sequence_names])
                if self.is_sparse else
                np.copy(aln_transpose[pi])
            )

            unique_letters = list(np.unique(pattern))
            if len(unique_letters) == 2 and self.ambiguous in unique_letters:
                other = [c for c in unique_letters if c != self.ambiguous][0]
                pattern[pattern == self.ambiguous] = other
                unique_letters = [other]

            str_pattern = "".join(pattern.astype('U'))
            if len(unique_letters) > 1:
                str_pattern += f"_{pi}"

            if str_pattern not in alignment_patterns:
                alignment_patterns[str_pattern] = (len(compressed_aln_transpose), [pi])
                compressed_aln_transpose.append(pattern)
            else:
                alignment_patterns[str_pattern][1].append(pi)

        for pattern_key, (new_idx, original_positions) in alignment_patterns.items():
            for pos in original_positions:
                self.full_to_compressed_sequence_map[pos] = new_idx
            self.compressed_to_full_sequence_map[new_idx] = np.array(original_positions, dtype=int)

        self._compressed_length = len(alignment_patterns)
        self._multiplicity = np.zeros(self._compressed_length, dtype=float)
        self.compressed_alignment = {}

        for i, name in enumerate(self.sequence_names):
            self.compressed_alignment[name] = np.array(
                [compressed_aln_transpose[p][i] for p in range(self._compressed_length)]
            )

        for i in range(self._compressed_length):
            self._multiplicity[i] = len(self.compressed_to_full_sequence_map[i])

        self.logger("SeqData: compression complete.", 1)
        return ttconf.SUCCESS
