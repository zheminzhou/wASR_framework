from collections import defaultdict
import numpy as np
import config as ttconf
from seq_utils import alphabets, profile_maps, alphabet_synonyms

def avg_transition(W, pi, gap_index=None):
    if gap_index is None:
        return np.einsum('i,ij,j', pi, W, pi)
    else:
        return (np.einsum('i,ij,j', pi, W, pi) - np.sum(pi * W[:, gap_index]) * pi[gap_index]) / (1 - pi[gap_index])

class GTR(object):
    def __init__(self, alphabet='nuc', prof_map=None, logger=None):
        self.debug = False
        self.is_site_specific = False
        
        if isinstance(alphabet, str):
            if alphabet not in alphabet_synonyms:
                raise AttributeError("Unknown alphabet type specified")
            tmp_alphabet = alphabet_synonyms[alphabet]
            self.alphabet = alphabets[tmp_alphabet]
            self.profile_map = profile_maps[tmp_alphabet]
        else:
            self.alphabet = np.array(alphabet)
            if prof_map is None:
                self.profile_map = {s: x for s, x in zip(self.alphabet, np.eye(len(self.alphabet)))}
            else:
                self.profile_map = {x if type(x) is str else x: k for x, k in prof_map.items()}

        self.state_index = {s: si for si, s in enumerate(self.alphabet)}
        self.state_index.update({s: si for si, s in enumerate(self.alphabet)})
        
        if logger is None:
            def logger_default(*args, **kwargs):
                if self.debug:
                    print(*args)
            self.logger = logger_default
        else:
            self.logger = logger

        self.ambiguous = None
        self.gap_index = None
        self.n_states = len(self.alphabet)
        self.assign_gap_and_ambiguous()

        self.logger("GTR: init with dummy values!", 3)
        self.v = None
        self.v_inv = None
        self.eigenvals = None
        self.assign_rates()

    def assign_gap_and_ambiguous(self):
        n_states = len(self.alphabet)
        self.logger("GTR: with alphabet: " + str([x for x in self.alphabet]), 1)
        
        if any([x.sum() == n_states for x in self.profile_map.values()]):
            amb_states = [c for c, x in self.profile_map.items() if x.sum() == n_states]
            self.ambiguous = 'N' if 'N' in amb_states else amb_states[0]
            self.logger("GTR: ambiguous character: " + self.ambiguous, 2)
        else:
            self.ambiguous = None

        try:
            self.gap_index = self.state_index['-']
        except:
            self.logger("GTR: no gap symbol!", 4, warn=True)
            self.gap_index = None

    @property
    def mu(self):
        return self._mu

    @property
    def Pi(self):
        return self._Pi

    @property
    def W(self):
        return self._W

    @W.setter
    def W(self, value):
        self.assign_rates(mu=self.mu, pi=self.Pi, W=value)

    @Pi.setter
    def Pi(self, value):
        self.assign_rates(mu=self.mu, pi=value, W=self.W)

    @mu.setter
    def mu(self, value):
        self.assign_rates(mu=value, pi=self.Pi, W=self.W)

    @property
    def Q(self):
        Q_tmp = (self.W * self.Pi).T
        Q_diag = -np.sum(Q_tmp, axis=0)
        np.fill_diagonal(Q_tmp, Q_diag)
        return Q_tmp

    def __str__(self):
        multi_site = len(self.Pi.shape) == 2
        
        if multi_site:
            eq_freq_str = "Average substitution rate (mu): " + str(np.round(self.average_rate, 6)) + '\n'
        else:
            eq_freq_str = "Substitution rate (mu): " + str(np.round(self.mu, 6)) + '\n'

        if not multi_site:
            eq_freq_str += "\nEquilibrium frequencies (pi_i):\n"
            for a, p in zip(self.alphabet, self.Pi):
                eq_freq_str += '  ' + a + ': ' + str(np.round(p, 4)) + '\n'

        W_str = "\nSymmetrized rates from j->i (W_ij):\n"
        W_str += '\t' + '\t'.join(self.alphabet) + '\n'
        for a, Wi in zip(self.alphabet, self.W):
            W_str += '  ' + a + '\t' + '\t'.join([str(np.round(max(0, p), 4)) for p in Wi]) + '\n'

        if not multi_site:
            Q_str = "\nActual rates from j->i (Q_ij):\n"
            Q_str += '\t' + '\t'.join(self.alphabet) + '\n'
            for a, Qi in zip(self.alphabet, self.Q):
                Q_str += '  ' + a + '\t' + '\t'.join([str(np.round(max(0, p), 4)) for p in Qi]) + '\n'

        return eq_freq_str + W_str + Q_str

    @staticmethod
    def from_file(gtr_fname):
        try:
            with open(gtr_fname) as f:
                alphabet = []
                pi = []
                while True:
                    line = f.readline()
                    if not line:
                        break
                    if line.strip().startswith("Substitution rate (mu):"):
                        mu = float(line.split(":")[1].strip())
                    elif line.strip().startswith("Equilibrium frequencies (pi_i):"):
                        line = f.readline()
                        while line.strip() != "":
                            alphabet.append(line.split(":")[0].strip())
                            pi.append(float(line.split(":")[1].strip()))
                            line = f.readline()
                        if not np.any([len(alphabet) == len(a) and np.all(np.array(alphabet) == a) for a in alphabets.values()]):
                            raise ValueError("GTR: was unable to read custom GTR model in " + str(gtr_fname) + " - Alphabet not recognized")
                    elif line.strip().startswith("Symmetrized rates from j->i (W_ij):"):
                        line = f.readline()
                        line = f.readline()
                        n = len(pi)
                        W = np.ones((n, n))
                        j = 0
                        while line.strip() != "":
                            values = line.split()
                            for i in range(n):
                                W[j, i] = float(values[i + 1])
                            j += 1
                            line = f.readline()
                        if j != n:
                            raise ValueError("GTR: was unable to read custom GTR model in " + str(gtr_fname) + " - Number of lines in W matrix does not match alphabet length")
                gtr = GTR.custom(mu, pi, W, alphabet=alphabet)
                return gtr
        except:
            raise ('GTR: was unable to read custom GTR model in ' + str(gtr_fname))

    def assign_rates(self, mu=1.0, pi=None, W=None):
        n = len(self.alphabet)
        self._mu = mu
        self.is_site_specific = False

        if pi is not None and len(pi) == n:
            Pi = np.array(pi)
        else:
            if pi is not None and len(pi) != n:
                self.logger("length of equilibrium frequency vector does not match alphabet length", 4, warn=True)
                self.logger("Ignoring input equilibrium frequencies", 4, warn=True)
            Pi = np.ones(shape=(n,))

        self._Pi = Pi / np.sum(Pi)

        if W is None or W.shape != (n, n):
            if (W is not None) and W.shape != (n, n):
                self.logger("Substitution matrix size does not match alphabet size", 4, warn=True)
                self.logger("Ignoring input substitution matrix", 4, warn=True)
            W = np.ones((n, n))
            np.fill_diagonal(W, 0.0)
            np.fill_diagonal(W, -W.sum(axis=0))
        else:
            W = np.array(W)

        self._W = 0.5 * (W + W.T)
        np.fill_diagonal(W, 0)
        average_rate = avg_transition(W, self.Pi, gap_index=self.gap_index)
        self._W = W / average_rate
        self._mu *= average_rate

        self._eig()

    @classmethod
    def custom(cls, mu=1.0, pi=None, W=None, **kwargs):
        gtr = cls(**kwargs)
        gtr.assign_rates(mu=mu, pi=pi, W=W)
        return gtr

    @staticmethod
    def standard(model, **kwargs):
        from nuc_models import JC69, K80, F81, HKY85, T92, TN93
        from aa_models import JTT92

        if model.lower() in ['jc', 'jc69', 'jukes-cantor', 'jukes-cantor69', 'jukescantor', 'jukescantor69']:
            model = JC69(**kwargs)
        elif model.lower() in ['k80', 'kimura80', 'kimura1980']:
            model = K80(**kwargs)
        elif model.lower() in ['f81', 'felsenstein81', 'felsenstein1981']:
            model = F81(**kwargs)
        elif model.lower() in ['hky', 'hky85', 'hky1985']:
            model = HKY85(**kwargs)
        elif model.lower() in ['t92', 'tamura92', 'tamura1992']:
            model = T92(**kwargs)
        elif model.lower() in ['tn93', 'tamura_nei_93', 'tamuranei93']:
            model = TN93(**kwargs)
        elif model.lower() in ['jtt', 'jtt92']:
            model = JTT92(**kwargs)
        else:
            raise KeyError("The GTR model '{}' is not in the list of available models.".format(model))

        model.mu = kwargs['mu'] if 'mu' in kwargs else 1.0
        return model

    @classmethod
    def random(cls, mu=1.0, alphabet='nuc', rng=None):
        if rng is None:
            rng = np.random.default_rng()

        alphabet = alphabets[alphabet]
        gtr = cls(alphabet)
        n = gtr.alphabet.shape[0]
        pi = 1.0 * rng.randint(0, 100, size=(n))
        W = 1.0 * rng.randint(0, 100, size=(n, n))

        gtr.assign_rates(mu=mu, pi=pi, W=W)
        return gtr

    @classmethod
    def infer(cls, nij, Ti, root_state, fixed_pi=None, pc=1.0, gap_limit=0.01, **kwargs):
        from scipy import linalg as LA
        gtr = cls(**kwargs)
        gtr.logger("GTR: model inference ", 1)
        dp = 1e-5
        Nit = 40
        pc_mat = pc * np.ones_like(nij)
        np.fill_diagonal(pc_mat, 0.0)
        np.fill_diagonal(nij, 0.0)
        count = 0
        pi_old = np.zeros_like(Ti)
        
        if fixed_pi is None:
            pi = np.ones_like(Ti)
        else:
            pi = np.copy(fixed_pi)
        pi /= pi.sum()
        
        W_ij = np.ones_like(nij)
        mu = (nij.sum() + pc) / (Ti.sum() + pc)
        
        while LA.norm(pi_old - pi) > dp and count < Nit:
            gtr.logger(' '.join(map(str, ['GTR inference iteration', count, 'change:', LA.norm(pi_old - pi)])), 3)
            count += 1
            pi_old = np.copy(pi)
            W_ij = (nij + nij.T + 2 * pc_mat) / mu / (np.outer(pi, Ti) + np.outer(Ti, pi) + ttconf.TINY_NUMBER + 2 * pc_mat)

            np.fill_diagonal(W_ij, 0)
            scale_factor = avg_transition(W_ij, pi, gap_index=gtr.gap_index)

            W_ij = W_ij / scale_factor
            if fixed_pi is None:
                pi = (np.sum(nij + pc_mat, axis=1) + root_state) / (ttconf.TINY_NUMBER + mu * np.dot(W_ij, Ti) + root_state.sum() + np.sum(pc_mat, axis=1))
                pi /= pi.sum()
                mu = (nij.sum() + pc) / (np.sum(pi * (W_ij.dot(Ti))) + pc)
            else:
                mu = (nij.sum() + pc) / (np.sum(pi * (W_ij.dot(pi))) * Ti.sum() + pc)

        if count >= Nit:
            gtr.logger('WARNING: maximum number of iterations has been reached in GTR inference', 3, warn=True)
            if LA.norm(pi_old - pi) > dp:
                gtr.logger('the iterative scheme has not converged', 3, warn=True)
            elif np.abs(1 - np.max(pi.sum(axis=0))) > dp:
                gtr.logger('the iterative scheme has converged, but proper normalization was not reached', 3, warn=True)
        
        if gtr.gap_index is not None:
            if pi[gtr.gap_index] < gap_limit:
                gtr.logger('The model allows for gaps which are estimated to occur at a low fraction of %1.3e' % pi[gtr.gap_index] +
                          ' this can potentially result in artificats.' +
                          ' gap fraction will be set to %1.4f' % gap_limit, 2, warn=True)
            pi[gtr.gap_index] = gap_limit
            pi /= pi.sum()

        gtr.assign_rates(mu=mu, W=W_ij, pi=pi)
        return gtr

    def _eig(self):
        self.eigenvals, self.v, self.v_inv = self._eig_single_site(self.W, self.Pi)

    def _eig_single_site(self, W, p):
        assert np.abs(np.diag(W).sum()) < 1e-10

        tmpp = np.sqrt(p)
        symQ = W * np.outer(tmpp, tmpp)
        np.fill_diagonal(symQ, -np.sum(W * p, axis=1))

        eigvals, eigvecs = np.linalg.eigh(symQ)
        tmp_v = eigvecs.T * tmpp
        one_norm = np.sum(np.abs(tmp_v), axis=1)
        return eigvals, tmp_v.T / one_norm, (eigvecs * one_norm).T / tmpp

    def state_pair(self, seq_p, seq_ch, pattern_multiplicity=None, ignore_gaps=False):
        if pattern_multiplicity is None:
            pattern_multiplicity = np.ones_like(seq_p, dtype=float)

        from collections import Counter
        if seq_ch.shape != seq_p.shape:
            raise ValueError("GTR.state_pair: Sequence lengths do not match!")

        if len(self.alphabet) < 10:
            pair_count = []
            bool_seqs_p = []
            bool_seqs_ch = []
            for seq, bs in [(seq_p, bool_seqs_p), (seq_ch, bool_seqs_ch)]:
                for ni, nuc in enumerate(self.alphabet):
                    bs.append(seq == nuc)

            for n1, nuc1 in enumerate(self.alphabet):
                if (self.gap_index is None) or (not ignore_gaps) or (n1 != self.gap_index):
                    for n2, nuc2 in enumerate(self.alphabet):
                        if (self.gap_index is None) or (not ignore_gaps) or (n2 != self.gap_index):
                            count = ((bool_seqs_p[n1] & bool_seqs_ch[n2]) * pattern_multiplicity).sum()
                            if count:
                                pair_count.append(((n1, n2), count))
        else:
            num_seqs = []
            for seq in [seq_p, seq_ch]:
                tmp = np.ones_like(seq, dtype=int)
                for ni, nuc in enumerate(self.alphabet):
                    tmp[seq == nuc] = ni
                num_seqs.append(tmp)
            pair_count = defaultdict(int)
            if ignore_gaps:
                for i in range(len(seq_p)):
                    if self.gap_index != num_seqs[0][i] and self.gap_index != num_seqs[1][i]:
                        pair_count[(num_seqs[0][i], num_seqs[1][i])] += pattern_multiplicity[i]
            else:
                for i in range(len(seq_p)):
                    pair_count[(num_seqs[0][i], num_seqs[1][i])] += pattern_multiplicity[i]
            pair_count = pair_count.items()

        return (np.array([x[0] for x in pair_count], dtype=int),
                np.array([x[1] for x in pair_count], dtype=int))

    def prob_t_compressed(self, seq_pair, multiplicity, t, return_log=False):
        if t < 0:
            logP = -ttconf.BIG_NUMBER
        else:
            tmp_eQT = self.expQt(t)
            logQt = np.log(np.maximum(tmp_eQT, ttconf.SUPERTINY_NUMBER))
            logP = np.sum(logQt[seq_pair[:, 1], seq_pair[:, 0]] * multiplicity)

        return logP if return_log else np.exp(logP)

    def prob_t(self, seq_p, seq_ch, t, pattern_multiplicity=None, return_log=False, ignore_gaps=True):
        seq_pair, multiplicity = self.state_pair(seq_p, seq_ch,
                                                 pattern_multiplicity=pattern_multiplicity, ignore_gaps=ignore_gaps)
        return self.prob_t_compressed(seq_pair, multiplicity, t, return_log=return_log)

    def optimal_t(self, seq_p, seq_ch, pattern_multiplicity=None, ignore_gaps=False):
        seq_pair, multiplicity = self.state_pair(seq_p, seq_ch,
                                                 pattern_multiplicity=pattern_multiplicity,
                                                 ignore_gaps=ignore_gaps)
        return self.optimal_t_compressed(seq_pair, multiplicity)

    def optimal_t_compressed(self, seq_pair, multiplicity, profiles=False, tol=1e-10):
        def _neg_prob(t, seq_pair, multiplicity):
            if profiles:
                res = -1.0 * self.prob_t_profiles(seq_pair, multiplicity, t**2, return_log=True)
                return res + np.exp(t**4 / 10000)
            else:
                return -1.0 * self.prob_t_compressed(seq_pair, multiplicity, t**2, return_log=True)

        hamming_distance = np.sum(multiplicity[seq_pair[:, 1] != seq_pair[:, 0]]) / np.sum(multiplicity)
        try:
            from scipy.optimize import minimize_scalar
            opt = minimize_scalar(_neg_prob,
                                  bracket=[-np.sqrt(ttconf.MAX_BRANCH_LENGTH), np.sqrt(hamming_distance), np.sqrt(ttconf.MAX_BRANCH_LENGTH)],
                                  args=(seq_pair, multiplicity), tol=tol, method='brent')
            new_len = opt["x"]**2
            if 'success' not in opt:
                opt['success'] = True
                self.logger("WARNING: the optimization result does not contain a 'success' flag:" + str(opt), 4, warn=True)
        except ImportError:
            import scipy
            print('legacy scipy', scipy.__version__)
            from scipy.optimize import fminbound
            new_len = fminbound(_neg_prob,
                                -np.sqrt(ttconf.MAX_BRANCH_LENGTH), np.sqrt(ttconf.MAX_BRANCH_LENGTH),
                                args=(seq_pair, multiplicity))
            new_len = new_len**2
            opt = {'success': True}

        if new_len > .9 * ttconf.MAX_BRANCH_LENGTH:
            self.logger("WARNING: GTR.optimal_t_compressed -- The branch length seems to be very long!", 4, warn=True)

        if opt["success"] != True:
            new_len = hamming_distance

        return new_len

    def prob_t_profiles(self, profile_pair, multiplicity, t, return_log=False, ignore_gaps=True):
        if t < 0:
            logP = -ttconf.BIG_NUMBER
        else:
            Qt = self.expQt(t)
            if len(Qt.shape) == 3:
                res = np.einsum('ai,ija,aj->a', profile_pair[1], Qt, profile_pair[0])
            else:
                res = np.einsum('ai,ij,aj->a', profile_pair[1], Qt, profile_pair[0])
            if ignore_gaps and (self.gap_index is not None):
                non_gap_frac = (1 - profile_pair[0][:, self.gap_index]) * (1 - profile_pair[1][:, self.gap_index])
                logP = np.sum(multiplicity * np.log(res + ttconf.SUPERTINY_NUMBER) * non_gap_frac)
            else:
                logP = np.sum(multiplicity * np.log(res + ttconf.SUPERTINY_NUMBER))

        return logP if return_log else np.exp(logP)

    def propagate_profile(self, profile, t, return_log=False):
        Qt = self.expQt(t)
        res = profile.dot(Qt)
        return np.log(res) if return_log else res

    def evolve(self, profile, t, return_log=False):
        Qt = self.expQt(t).T
        res = profile.dot(Qt)
        return np.log(res) if return_log else res

    def _exp_lt(self, t):
        log_val = self.mu * t * self.eigenvals
        if any(i > 10 for i in log_val):
            raise ValueError("Error in computing exp(Q * t): Q has positive eigenvalues or the branch length t is too large. "
                             "This is most likely caused by incorrect input data.")
        return np.exp(log_val)

    def expQt(self, t):
        eLambdaT = np.diag(self._exp_lt(t))
        Qs = self.v.dot(eLambdaT.dot(self.v_inv))
        return np.maximum(0, Qs)

    def expQs(self, s):
        return self.expQt(s**2)

    def expQsds(self, s):
        lambda_eLambdaT = np.diag(2.0 * self._exp_lt(s**2) * self.eigenvals * s)
        return self.v.dot(lambda_eLambdaT.dot(self.v_inv))

    def sequence_logLH(self, seq, pattern_multiplicity=None):
        if pattern_multiplicity is None:
            pattern_multiplicity = np.ones_like(seq, dtype=float)
        return np.sum([np.sum((seq == state) * pattern_multiplicity * np.log(self.Pi[si]))
                       for si, state in enumerate(self.alphabet)])

    def average_rate(self):
        return self.mu * avg_transition(self.W, self.Pi, gap_index=self.gap_index)

    def save_to_npz(self, outfile):
        full_gtr = self.mu * np.dot(self.Pi, self.W)
        desc = np.array(["GTR matrix description\n", "Substitution rate: " + str(self.mu)])
        np.savez(outfile, description=desc, full_gtr=full_gtr, char_dist=self.Pi, flow_matrix=self.W)

if __name__ == "__main__":
    pass