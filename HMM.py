import numpy as np
from collections import Counter
import Pbas.Reader as rd

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
    def __init__(self, emission_probabilities, type=None, transition_probabilities=[-np.inf, -np.inf, -np.inf, -np.inf]):
        """
        :param emission_probabilities: A dictionary of the form
        {l: e(l)} for each l=possible emitted letter and e(l) the
        emission probability in log
        """
        self.log_emission_probabilities = emission_probabilities  # dict of AA probability
        self.log_transition_probabilities = transition_probabilities  # array size 4 - [M, I, D, B]
        self.type = type  # B/M/I/D

    def update_log_emissions(self, letter, new_log_prob):
        self.log_emission_probabilities[letter] = new_log_prob

    def update_log_transition(self, other, new_log_prob):
        self.log_transition_probabilities[other] = new_log_prob

class HMM_Model:
    def __init__(self, states):
        self.states = states
        self.start_state = self.states[0, 0]

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
    num_states = (len(states) - 2) * 3 + 2
    motif_gap_or_nuc = motif_matrix != GAP  # matrix with size like motif_mat represent if there is nuc or gap
    trans_mat = np.zeros((num_states, num_states), dtype=np.float64)
    states[0, 0].update_log_transition(3, log_bg_t)
    states[0, 0].update_log_transition(0, log_bg_t * sum(motif_matrix[:, 0] != GAP))
    states[0, 0].update_log_transition(2, 1 - log_bg_t * sum(motif_matrix[:, 0] == GAP))
    s = motif_s.find(motif_s == 1)
    trans_mat[0, 1] = np.sum(motif_s[:, s], axis=1)
    for i in range(len(motif_i)):
        m_to_m, m_to_i, m_to_d = -np.inf, -np.inf, -np.inf
        i_to_i, i_to_d, i_to_m = -np.inf, -np.inf, -np.inf
        d_to_i, d_to_d, d_to_m = -np.inf, -np.inf, -np.inf
        if i != len(motif_i) - 1:
            con1 = motif_matrix[:, motif_i[i]] != GAP
            con2 = motif_matrix[:, motif_i[i]+1:motif_matrix[i+1]] == GAP
            con3 = motif_matrix[:, motif_i[i+1]] != GAP
            m_to_m = np.log((np.sum(con1&con2&con3)+1)/(motif_matrix.shape[0]+1))
            m_to_i = np.log((np.sum(con1&~con2)+1)/(motif_matrix.shape[0]+1))
            m_to_d = np.log((np.sum(con1&con2&~con3)+1)/(motif_matrix.shape[0]+1))
        if motif_i[i]+1 != motif_i[i+1]:
            temp_mat = motif_matrix[motif_i[i]+1:motif_i[i+1]]
            temp_mat = temp_mat[np.sum(temp_mat != GAP, 0) > 1, :]
            i_to_i = np.sum(temp_mat != GAP) - len(temp_mat)
            i_to_d = np.log((np.sum(~con2&~con3)+1)/(motif_matrix.shape[0]+1))
            i_to_m = np.log((np.sum(~con2&con3)+1)/(motif_matrix.shape[0]+1))
        if i != len(motif_i) - 1:
            d_to_d = np.log((np.sum(~con1&con2&~con3)+1)/(motif_matrix.shape[0]+1))
            d_to_m = np.log((np.sum(~con1&con2&con3)+1)/(motif_matrix.shape[0]+1))
            d_to_i = np.log((np.sum(~con1&con2)+1)/(motif_matrix.shape[0]+1))
        else:  # we got to the last position in the motif
            states[i, 0].update_log_transition(3, np.log(1 - np.exp(m_to_i)))
            states[i, 1].update_log_transition(3, np.log(1 - np.exp(i_to_i)))
            states[i, 2].update_log_transition(3, np.log(1 - np.exp(d_to_i)))
        states[i, 0].update_log_transition(0, m_to_m)
        states[i, 0].update_log_transition(1, m_to_i)
        states[i, 0].update_log_transition(2, m_to_d)
        states[i, 1].update_log_transition(0, i_to_m)
        states[i, 1].update_log_transition(1, i_to_i)
        states[i, 1].update_log_transition(2, i_to_d)
        states[i, 2].update_log_transition(0, d_to_m)
        states[i, 2].update_log_transition(1, d_to_i)
        states[i, 2].update_log_transition(2, d_to_d)


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
    model = init_hmm_model(motif_matrix, log_bg_e, log_bg_t)


if __name__ == '__main__':
    main()