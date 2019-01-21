import numpy as np
from collections import Counter
import Pbas.Reader as rd
from scipy.misc import logsumexp

# -------------------- Constants --------------------
MAX_EM_ITERS = 100
MOTIF_STATE = 'M'
BEGIN_STATE = 'B'
END_STATE = 'E'
DELETION_STATE = 'D'
INSERTION_STATE = 'I'
GAP = '.'
LETTERS = ['A', 'T', 'C', 'G']
BACKGROUND_SYMBOL = '0'
MOTIF_SYMBOL = '1'
M = 1
I = 0

class HMMState:
    def __init__(self, emission_probabilities, type=None):
        """
        :param emission_probabilities: A dictionary of the form
        {l: e(l)} for each l=possible emitted letter and e(l) the
        emission probability in log
        """
        self.log_emission_probabilities = emission_probabilities  # dict of AA probability
        self.log_transition_probabilities = [-np.inf, -np.inf, -np.inf, -np.inf]  # array size 4 - [M, I, D, B]
        self.type = type  # B/M/I/D

    def update_log_emissions(self, letter, new_log_prob):
        self.log_emission_probabilities[letter] = new_log_prob

    def update_log_transition(self, other, new_log_prob):
        self.log_transition_probabilities[other] = new_log_prob

    def check_prob(self):
        print(np.exp(logsumexp(self.log_transition_probabilities)))

class HMM_Model:
    def __init__(self, states):
        self.states = states
        self.start_state = self.states[0][0]

# -------------------- Helper methods --------------------


def build_states(motif_matrix, motif_s, motif_i, log_bg_e):
    """
    Initialize the state and calculate their emissions
    :param motif_matrix: matrix of all the motif
    :param motif_s: vector of T or F. T in i if the position i is in the motif, F otherwise
    :return: array of HMMState object represent states of the model
    """
    deletion_emition = {GAP: 0}
    states = [[HMMState(log_bg_e, BEGIN_STATE)]]
    for i in range(len(motif_i)-1):
        emi = Counter(np.ravel(motif_matrix[:, motif_i[i]:motif_i[i+1]]))
        if GAP in emi.keys():
            emi.pop(GAP)
        ls = np.log(sum(emi.values()))
        emi = {x: emi[x]-ls for x in emi.keys()}
        states.append([HMMState(emi, MOTIF_STATE), HMMState(emi, INSERTION_STATE), HMMState(deletion_emition, DELETION_STATE)])
    emi = Counter(np.ravel(motif_matrix[:, motif_i[-1]::]))
    if GAP in emi.keys():
        emi.pop(GAP)
    ls = np.log(sum(emi.values()))
    emi = {x: emi[x]-ls for x in emi.keys()}
    states.append([HMMState(emi, MOTIF_STATE), HMMState(emi, INSERTION_STATE), HMMState(deletion_emition, DELETION_STATE)])
    return states


def build_transition(motif_matrix, motif_s, motif_i, log_bg_t, states):
    """
    Update the transition of each state
    :param motif_matrix: matrix of all the motif
    :param motif_s: vector of T or F. T in i if the position i is in the motif, F otherwise
    :param motif_i: vector represent the indexes of the motif
    :param log_bg_t: the log transition to stay on background state
    :param states: array of HMMState objects
    """
    states[0][0].update_log_transition(3, log_bg_t)
    states[0][0].update_log_transition(0,  np.log(1 - np.exp(log_bg_t)) + np.log(sum(motif_matrix[:, motif_i[0]] != GAP)))
    states[0][0].update_log_transition(2, np.log(1 - np.exp(log_bg_t)) + np.log(sum(motif_matrix[:, motif_i[0]] == GAP)))
    for i in range(len(motif_i)):
        init_pro = np.log(1/(motif_matrix.shape[0]+3))
        in_m, in_i, in_d = 3, 3, 3
        m_to_m, m_to_i, m_to_d = init_pro, init_pro, init_pro
        i_to_i, i_to_d, i_to_m = init_pro, init_pro, init_pro
        d_to_i, d_to_d, d_to_m = init_pro, init_pro, init_pro
        next_i = -1 if i == len(motif_i) - 1 else motif_i[i+1]
        con1 = motif_matrix[:, motif_i[i]] != GAP
        con3 = motif_matrix[:, next_i] != GAP
        con2 = np.all(motif_matrix[:, motif_i[i]+1:next_i] == GAP, axis=1) if motif_i[i]+1 != next_i else np.array([True]*len(con1))
        ncon2 = np.any(motif_matrix[:, motif_i[i]+1:next_i] != GAP, axis=1) if motif_i[i]+1 != next_i else np.array([False]*len(con1))
        if motif_i[i]+1 != next_i:
            temp_mat = motif_matrix[:, motif_i[i]+1:next_i]
            temp_mat = temp_mat[np.sum(temp_mat != GAP, axis=1) > 1, :]
            m_to_i = np.log((np.sum(con1&ncon2)+1))
            i_to_i = np.log((np.sum(temp_mat != GAP) - len(temp_mat)+1)/
                            (motif_matrix.shape[0]*(next_i - motif_i[i] - 1)+3))
            d_to_i = np.log((np.sum(~con1&con2)+1))
            in_m += np.sum(con1&ncon2)
            in_i += np.sum(temp_mat != GAP) - len(temp_mat)
            in_d += np.sum(~con1&con2)
            if i != len(motif_i) - 1:
                i_to_d = np.log((np.sum(ncon2&~con3)+1))
                i_to_m = np.log((np.sum(ncon2&con3)+1))
                in_i += np.sum(ncon2&~con3)
                in_i += np.sum(ncon2&con3)
        if i != len(motif_i) - 1:
            m_to_m = np.log((np.sum(con1&con2&con3)+1))
            m_to_d = np.log((np.sum(con1&con2&~con3)+1))
            in_m += np.sum(con1&con2&con3) + np.sum(con1&con2&~con3)
            d_to_d = np.log((np.sum(~con1&con2&~con3)+1))
            d_to_m = np.log((np.sum(~con1&con2&con3)+1))
            in_d += np.sum(~con1&con2&~con3) + np.sum(~con1&con2&con3)
        else:  # we got to the last position in the motif
            states[i+1][0].update_log_transition(3, np.log(1 - np.exp(m_to_i)))
            states[i+1][1].update_log_transition(3, np.log(1 - np.exp(i_to_i)))
            states[i+1][2].update_log_transition(3, np.log(1 - np.exp(d_to_i)))
        states[i+1][0].update_log_transition(0, m_to_m - np.log(in_m))
        states[i+1][0].update_log_transition(1, m_to_i - np.log(in_m))
        states[i+1][0].update_log_transition(2, m_to_d - np.log(in_m))
        states[i+1][1].update_log_transition(0, i_to_m - np.log(in_m))
        states[i+1][1].update_log_transition(1, i_to_i - np.log(in_i))
        states[i+1][1].update_log_transition(2, i_to_d - np.log(in_i))
        states[i+1][2].update_log_transition(0, d_to_m - np.log(in_d))
        states[i+1][2].update_log_transition(1, d_to_i - np.log(in_d))
        states[i+1][2].update_log_transition(2, d_to_d - np.log(in_d))
        # test
        for j in range(3):
            if j == 0:
                print("i:", i, "motiv_i:", motif_i[i])
                if i != len(motif_i) - 1:
                    print("i+1:", i, "motiv_i+1:", motif_i[i+1])
            print("state: ", j)
            states[i+1][j].check_prob()


def init_hmm_model(motif_matrix, log_bg_e, log_bg_t):
    """
    Initialize a new HMMModel object in accordance with the
    requirements of 76558 ex2
    :param p: Initial transfer probability into a motif state
    :param initial_emissions_tsv_path:
    :return: The HMMModel
    """
    motif_s = motif_matrix != GAP
    motif_s = np.sum(motif_s, axis=0)
    motif_s = motif_s > len(motif_matrix) / 2
    motif_i = np.where(motif_s == 1)[0]
    states = build_states(motif_matrix, motif_s, motif_i, log_bg_e)
    build_transition(motif_matrix, motif_s, motif_i,log_bg_t, states)
    return HMM_Model(states)

# -------------------- Main --------------------


def main():
    """
    Run HMM
    """
    reader = rd.Reader("PF00096", "full", offline=True, motif=True, save_file=True)
    motif_matrix = reader.get_fasta()
    log_bg_e = reader.get_bg_e()
    log_bg_t = reader.get_bg_t()
    init_hmm_model(motif_matrix, log_bg_e, log_bg_t)


if __name__ == '__main__':
    main()
