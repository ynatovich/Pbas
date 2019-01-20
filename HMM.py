import numpy as np
from collections import Counter
import Reader

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


class HMMState:
    def __init__(self, emission_probabilities, type=None):
        """
        :param emission_probabilities: A dictionary of the form
        {l: e(l)} for each l=possible emitted letter and e(l) the
        emission probability in log
        """
        self.log_emission_probabilities = emission_probabilities
        self.type = type

    def update_log_emissions(self, letter, new_log_prob):
        self.log_emission_probabilities[letter] = new_log_prob


class HMM_Model:
    def __init__(self, states):
        self.states = states
        self.start_state = self.states[0]
        self.transition_table = None
        self.log_transition_table = None

    def set_transition_table(self, table):
        self.log_transition_table = np.log(table)

    def update_log_transition_table(self, i, j, value):
        self.log_transition_table[i, j] = value

# -------------------- Helper methods --------------------


def build_states(motif_matrix, motif_s, log_back_emi):
    """

    :param motif_matrix:
    :param motif_s:
    :return:
    """
    deletion_emition = {GAP: 0}
    motif_i = np.where(motif_s == 1)
    states = [[HMMState(log_back_emi, BEGIN_STATE)]]
    for i in range(len(motif_i)):
        emi = Counter(motif_matrix[:, motif_i[i]:motif_i[i+1]])
        emi.pop(GAP)
        s = sum(emi.values())
        emi = {x: np.log(emi[x]/s) for x in emi.keys()}
        states.append([HMMState(emi, MOTIF_STATE), HMMState(emi, INSERTION_STATE), HMMState(deletion_emition, DELETION_STATE)])
    states.append([HMMState(log_back_emi, END_STATE)])
    return states


def init_hmm_model(motif_matrix, log_back_emi):
    """
    Initialize a new HMMModel object in accordance with the
    requirements of 76558 ex2
    :param p: Initial transfer probability into a motif state
    :param initial_emissions_tsv_path:
    :return: The HMMModel
    """
    motif_s = motif_matrix != GAP
    motif_s = np.sum(motif_s, axis=1)
    motif_s = motif_s > len(motif_matrix) / 2
    states = build_states(motif_matrix, motif_s, log_back_emi)
    



    # emission_table
    # motif_states = [HMMState(emission_probabilities=emissions) for emissions in emission_table]
    # model = HMM_Model(states=[B] + motif_states)
    # transition_table = np.zeros(shape=(len(motif_states) + 1, len(motif_states) + 1))
    # transition_table[0, 0] = 1 - p
    # transition_table[0, 1] = p
    # for i in range(1, len(model.states) - 1):
    #     transition_table[i, i + 1] = 1
    # transition_table[-1, 0] = 1
    # model.set_transition_table(transition_table)
    # return model


# -------------------- Main --------------------

def main():
    """
    Run 76558 ex2 part I
    """
    reader = Reader("PF00096", "full")
    motif_matrix = reader.get_fasta()
    # log_back_emi =
    model = init_hmm_model(motif_matrix, log_back_emi)

if __name__ == '__main__':
    main()
