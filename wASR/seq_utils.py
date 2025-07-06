import numpy as np
from Bio import Seq, SeqRecord
from collections import defaultdict

alphabet_synonyms = {
    'nuc': 'nuc', 'nucleotide': 'nuc', 'aa': 'aa', 'aminoacid': 'aa',
    'nuc_nogap': 'nuc_nogap', 'nucleotide_nogap': 'nuc_nogap',
    'aa_nogap': 'aa_nogap', 'aminoacid_nogap': 'aa_nogap',
    'DNA': 'nuc', 'DNA_nogap': 'nuc_nogap'
}

alphabets = {
    "nuc": np.array(['A', 'C', 'G', 'T', '-']),
    "nuc_nogap": np.array(['A', 'C', 'G', 'T']),
    "aa": np.array(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '*', '-']),
    "aa_nogap": np.array(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'])
}

def create_profile_maps():
    nuc_profiles = {
        'A': [1, 0, 0, 0, 0], 'C': [0, 1, 0, 0, 0], 'G': [0, 0, 1, 0, 0], 'T': [0, 0, 0, 1, 0], '-': [0, 0, 0, 0, 1],
        'N': [1, 1, 1, 1, 1], 'X': [1, 1, 1, 1, 1], 'R': [1, 0, 1, 0, 0], 'Y': [0, 1, 0, 1, 0], 'S': [0, 1, 1, 0, 0],
        'W': [1, 0, 0, 1, 0], 'K': [0, 0, 1, 1, 0], 'M': [1, 1, 0, 0, 0], 'D': [1, 0, 1, 1, 0], 'H': [1, 1, 0, 1, 0],
        'B': [0, 1, 1, 1, 0], 'V': [1, 1, 1, 0, 0]
    }
    
    nuc_nogap_profiles = {
        'A': [1, 0, 0, 0], 'C': [0, 1, 0, 0], 'G': [0, 0, 1, 0], 'T': [0, 0, 0, 1], '-': [1, 1, 1, 1],
        'N': [1, 1, 1, 1], 'X': [1, 1, 1, 1], 'R': [1, 0, 1, 0], 'Y': [0, 1, 0, 1], 'S': [0, 1, 1, 0],
        'W': [1, 0, 0, 1], 'K': [0, 0, 1, 1], 'M': [1, 1, 0, 0], 'D': [1, 0, 1, 1], 'H': [1, 1, 0, 1],
        'B': [0, 1, 1, 1], 'V': [1, 1, 1, 0]
    }
    
    aa_order = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '*', '-']
    aa_profiles = {}
    for i, aa in enumerate(aa_order):
        profile = [0] * 22
        profile[i] = 1
        aa_profiles[aa] = profile
    
    aa_profiles['X'] = [1] * 22
    aa_profiles['B'] = [0] * 22
    aa_profiles['B'][2] = 1  # D
    aa_profiles['B'][11] = 1  # N
    aa_profiles['Z'] = [0] * 22
    aa_profiles['Z'][3] = 1  # E
    aa_profiles['Z'][13] = 1  # Q
    
    aa_nogap_order = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    aa_nogap_profiles = {}
    for i, aa in enumerate(aa_nogap_order):
        profile = [0] * 20
        profile[i] = 1
        aa_nogap_profiles[aa] = profile
    
    aa_nogap_profiles['X'] = [1] * 20
    aa_nogap_profiles['B'] = [0] * 20
    aa_nogap_profiles['B'][2] = 1  # D
    aa_nogap_profiles['B'][11] = 1  # N
    aa_nogap_profiles['Z'] = [0] * 20
    aa_nogap_profiles['Z'][3] = 1  # E
    aa_nogap_profiles['Z'][13] = 1  # Q
    
    return {
        'nuc': {k: np.array(v, dtype=float) for k, v in nuc_profiles.items()},
        'nuc_nogap': {k: np.array(v, dtype=float) for k, v in nuc_nogap_profiles.items()},
        'aa': {k: np.array(v, dtype=float) for k, v in aa_profiles.items()},
        'aa_nogap': {k: np.array(v, dtype=float) for k, v in aa_nogap_profiles.items()}
    }

profile_maps = create_profile_maps()

def extend_profile(gtr, aln, logger=None):
    unique_chars = set()
    for seq in aln:
        unique_chars.update(np.unique(seq))
    
    for c in unique_chars:
        if c not in gtr.profile_map:
            gtr.profile_map[c] = np.ones(gtr.n_states)
            if logger:
                logger(f"WARNING: character {c} is unknown. Treating it as missing information", 1, warn=True)

def guess_alphabet(aln):
    total = sum(len(seq) for seq in aln)
    nuc_count = sum(np.sum(np.isin(seq, list('acgtACGT-N'))) for seq in aln)
    return 'nuc' if nuc_count > 0.9 * total else 'aa'

def seq2array(seq, word_length=1, convert_upper=False, fill_overhangs=False, ambiguous='N'):
    if isinstance(seq, str):
        seq_str = seq
    elif isinstance(seq, Seq.Seq):
        seq_str = str(seq)
    elif isinstance(seq, SeqRecord.SeqRecord):
        seq_str = str(seq.seq)
    else:
        raise TypeError(f"seq2array: sequence must be Bio.Seq, Bio.SeqRecord, or string. Got {type(seq)}")

    if convert_upper:
        seq_str = seq_str.upper()

    if word_length == 1:
        seq_array = np.array(list(seq_str))
    else:
        if len(seq_str) % word_length:
            raise ValueError("sequence length has to be multiple of word length")
        seq_array = np.array([seq_str[i*word_length:(i+1)*word_length]
                              for i in range(len(seq_str)//word_length)])

    if fill_overhangs:
        gaps = np.where(seq_array != '-')[0]
        if len(gaps):
            seq_array[:gaps[0]] = ambiguous
            seq_array[gaps[-1]+1:] = ambiguous
        else:
            seq_array[:] = ambiguous

    return seq_array

def seq2prof(seq, profile_map):
    return np.array([profile_map[k] for k in seq])

def prof2seq(profile, gtr, sample_from_prof=False, normalize=True, rng=None):
    if rng is None:
        rng = np.random.default_rng()
    
    if normalize:
        tmp_profile, _ = normalize_profile(profile, return_offset=False)
    else:
        tmp_profile = profile

    if sample_from_prof:
        cumdis = tmp_profile.cumsum(axis=1).T
        randnum = rng.random(size=cumdis.shape[1])
        idx = np.argmax(cumdis >= randnum, axis=0)
    else:
        idx = tmp_profile.argmax(axis=1)
    
    seq = gtr.alphabet[idx]
    prof_values = tmp_profile[np.arange(tmp_profile.shape[0]), idx]
    
    return seq, prof_values, idx

def normalize_profile(in_profile, log=False, return_offset=True):
    if log:
        tmp_prefactor = in_profile.max(axis=1)
        tmp_prof = np.exp(in_profile.T - tmp_prefactor).T
    else:
        tmp_prefactor = 0.0
        tmp_prof = in_profile

    norm_vector = tmp_prof.sum(axis=1)
    normalized = np.einsum('ai,a->ai', tmp_prof, 1.0/norm_vector)
    
    if return_offset:
        offset = np.log(norm_vector) + tmp_prefactor
        return normalized, offset
    else:
        return normalized, None

def process_sparse_alignment(aln, ref, ambiguous_char):
    nseq = len(aln)
    inv_map = defaultdict(list)
    
    for k, v in aln.items():
        for pos, bs in v.items():
            inv_map[pos].append(bs)

    nonref_positions = np.sort(list(inv_map.keys()))
    constant_up_to_ambiguous = []
    nonref_const = []
    nonref_alleles = []
    ambiguous_const = []
    variable_pos = []
    
    for pos, bs in inv_map.items():
        bases = list(np.unique(bs))
        if len(bs) == nseq:
            if (len(bases) <= 2 and ambiguous_char in bases) or len(bases) == 1:
                nonref_const.append(pos)
                if len(bases) == 1:
                    nonref_alleles.append(bases[0])
                else:
                    nonref_alleles.append([x for x in bases if x != ambiguous_char][0])
                
                if ambiguous_char in bases:
                    constant_up_to_ambiguous.append(pos)
            else:
                variable_pos.append(pos)
        else:
            if len(bases) == 1 and bases[0] == ambiguous_char:
                ambiguous_const.append(pos)
                constant_up_to_ambiguous.append(pos)
            else:
                variable_pos.append(pos)

    refMod = np.copy(ref)
    refMod[nonref_const] = nonref_alleles
    states = np.unique(refMod)
    refMod[variable_pos] = '.'

    constant_columns = []
    constant_patterns = {}
    
    for base in states:
        if base == ambiguous_char:
            continue
        p = np.repeat(base, nseq)
        pos = list(np.where(refMod == base)[0])
        if len(pos):
            constant_patterns["".join(p.astype('U'))] = [len(constant_columns), pos]
            constant_columns.append(p)

    return {
        "constant_columns": constant_columns,
        "constant_patterns": constant_patterns,
        "variable_positions": variable_pos,
        "nonref_positions": nonref_positions,
        "constant_up_to_ambiguous": constant_up_to_ambiguous
    }