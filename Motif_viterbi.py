import argparse
import numpy as np
import pandas as pd
import HMM
from HMM import HMM_Model, HMMState
import TCM_EM

# -------------------- Constants --------------------
LETTERS_PER_LINE = 50
DEBUG = False
# -------------------- Helper methods --------------------
def init_hmm_model(p, initial_emissions_tsv_path):
    """
    Initialize a new HMMModel object in accordance with the
    requirements of 76558 ex2
    :param p: Initial transfer probability into a motif state
    :param initial_emissions_tsv_path:
    :return: The HMMModel
    """
    B = HMMState(emission_probabilities=dict(A=0.25, T=0.25, C=0.25, G=0.25))
    emission_table = pd.read_csv(initial_emissions_tsv_path, sep="\t").to_dict(
        'records')
    motif_states = [HMMState(emission_probabilities=emissions) for emissions in
                    emission_table]
    model = HMM_Model(states=[B] + motif_states)
    transition_table = np.zeros(
        shape=(len(motif_states) + 1, len(motif_states) + 1))
    transition_table[0, 0] = 1 - p
    transition_table[0, 1] = p
    for i in range(1, len(model.states) - 1):
        transition_table[i, i + 1] = 1
    transition_table[-1, 0] = 1
    model.set_transition_table(transition_table)
    return model


def print_paths(seq, *paths):
    """
    Given a source sequence and 1 or more paths, print them according
     to the requirements of 76558 ex2
    """
    layout_paths = [None] * len(paths)
    for i, path in enumerate(paths):
        path = [HMM.BACKGROUND_SYMBOL if l == 0 else HMM.MOTIF_SYMBOL for l in path]
        layout_paths[i] = [path[i:i + LETTERS_PER_LINE] for i in
                           range(0, len(path), LETTERS_PER_LINE)]
    s = [seq[i:i + LETTERS_PER_LINE] for i in
         range(0, len(seq), LETTERS_PER_LINE)]
    for i in range(len(s)):
        print("".join(s[i]))
        for path in layout_paths:
            print("".join(path[i]))
        print()


# -------------------- Main --------------------

def main():
    """
    Run 76558 ex2 part I
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('seq',
        help='Sequence to analyze')
    parser.add_argument('initial_emissions',
        help='.tsv file with the initial emission probabilities')
    parser.add_argument('p',
        help='Transition probability from background state to motif state')
    command_args = parser.parse_args()
    model = init_hmm_model(float(command_args.p),
                           command_args.initial_emissions)

    viterbi_path, viterbi_log_likelihood = model.viterbi(command_args.seq)
    map_path, map_log_likelihood = model.MAP_path(command_args.seq)

    print("viterbi_log_likelihood={}\nMAP_log_likelihood={}"
          .format(viterbi_log_likelihood, map_log_likelihood))
    print_paths(command_args.seq, viterbi_path, map_path)

    # -------------------- Test --------------------
    if DEBUG:
        import random
        motif = "TAAGCTT"
        S = ['A', 'C', 'G', 'T']
        seq = ''.join([''.join(random.choices(S, k=43)) + motif for i in range(5)])
        p = 0.005
        model = init_hmm_model(p, "initial_emission.tsv")
        viterbi_path, viterbi_log_likelihood = model.viterbi(seq)
        map_path, map_log_likelihood = model.MAP_path(seq)
        print("viterbi_log_likelihood={}\nMAP_log_likelihood={}"
              .format(viterbi_log_likelihood, map_log_likelihood))
        print_paths(seq, viterbi_path, map_path)


if __name__ == '__main__':
    main()
