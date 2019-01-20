import argparse
import numpy as np
import HMM
from HMM import HMM_Model, HMMState, MOTIF_STATE, BACKGROUND_STATE
from itertools import groupby
from collections import Counter

PLOT_MODE = True

def fastaread(fasta_name):
    """ Read fasta files """
    f = open(fasta_name)
    faiter = (x[1] for x in groupby(f, lambda line: line.startswith(">")))
    for header in faiter:
        header = next(header)[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq


def most_occurring_words(fasta_name, word_len, n_words):
    """
    :param fasta_name:
    :param word_len: word length
    :param n_words: number of words
    :return: The n most occurring k-long words
    """
    cntr = Counter()
    n_possible = 0
    for _, seq in fastaread(fasta_name):
        cntr.update([seq[i:i + word_len] for i in range(len(seq) - word_len + 1)])
        n_possible += len(seq) - word_len + 1
    return cntr.most_common(n_words), n_possible




def init_hmm_model(motif, p, alpha):
    """
    Initialize a new HMMModel object in accordance with the
    requirements of 76558 ex2 part II
    :param motif:
    :param p: Initial transfer probability into a motif state
    :param alpha:
    :return: The created HMM model
    """
    B = HMMState(emission_probabilities=dict(A=0.25, T=0.25, C=0.25, G=0.25),
                 type=BACKGROUND_STATE)
    motif_states = [None]*len(motif)
    motif_state_emissions = dict(A=0.25, T=0.25, C=0.25, G=0.25)
    for i, l in enumerate(motif):
        motif_state_emissions = dict.fromkeys(motif_state_emissions, alpha)
        motif_state_emissions[l] = 1 - 3*alpha
        motif_states[i] = HMMState(emission_probabilities=motif_state_emissions,
                                   type=MOTIF_STATE)

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


def output_motif_txt(model, fasta_seqs, index):
    """
    Output a .txt file named motif<index>.txt which describes the HMM model after running TCM_EM
     for the given motif (Full format specifications  in 76558 ex2 part II)
    :param model:
    :param fasta_seqs: An array of sequence tuples to evaluate the model on, each consists of
    the sequence name followed by the sequence itself
    :param index:
    """
    txt = open("motif" + str(index) + ".txt", "w")
    for letter in HMM.LETTERS:
        for state in model.states:
            if state.type != HMM.MOTIF_STATE:
                continue
            txt.write("%.2f  " % state.emission_probabilities[letter])
        txt.write('\n')
    for seq in fasta_seqs:
        path = model.viterbi(seq[1])
        starting_locations = [i for (i, x) in enumerate(path[0]) if x == 1]
        starting_locations_string = "0" if len(starting_locations) == 0 else str(starting_locations)[1:-1]
        txt.write(seq[0] + " "*(10 - len(seq[0])))
        txt.write(starting_locations_string)
        txt.write('\n')
    txt.close()



def output_histories(ll_histories, motifs):
    """
    Output a .tab file with the log_likelihood as a function of number of iterations for each of
    the given motifs (Full format specifications in 76558 ex2 part II)
    :param ll_histories:
    :param motifs:
    :return:
    """
    txt = open("History.tab", "w")
    txt.write(" ")
    for motif in motifs:
        txt.write(motif[0] + " "*(9 - len(motif[0])))
    txt.write("\n")
    max_len = max([len(h[1]) for h in ll_histories])
    for i in range(max_len):
        for hist in ll_histories:
            if i < len(hist[1]):
                txt.write(("%.2f  " % hist[1][i]))
            else:
                txt.write("         ")
        txt.write("\n")
    txt.close()

def main():
    """
    Run 76558 ex2 part II
    """
    # --------------- Command-line Argument handling ---------------
    parser = argparse.ArgumentParser()
    parser.add_argument('trainSeq',
                        help='File name of the list of sequences in which to '
                             'search for the motifs')
    parser.add_argument('k',
                        help='The length of the motifs')
    parser.add_argument('convergenceThr',
                        help='The stopping condition for the EM is a threshold'
                             ' on the improvement of the log likelihood.')
    parser.add_argument('L1',
                        help='Number of seeds used for initialization.')
    parser.add_argument('alpha',
                        help='Softening parameter used to create the initial '
                             'profile')
    command_args = parser.parse_args()


    # --------------- Get the sequences and motifs ---------------
    motifs, n_possible = most_occurring_words(fasta_name=command_args.trainSeq,
                                  word_len=int(command_args.k),
                                  n_words=int(command_args.L1))
    fasta_seqs = list(fastaread(command_args.trainSeq))
    n = len(fasta_seqs)
    seq_names = [None] * n
    seqs = [None] * n
    for i in range(n):
        seq_names[i] = fasta_seqs[i][0]
        seqs[i] = fasta_seqs[i][1]
    ll_histories = [None] * int(command_args.L1)

    # --------------- For each motif, setup a HMM model and perform TCM-EM---------------
    for i, motif_t in enumerate(motifs):
        motif, n_occurrences = motif_t
        p = n_occurrences / n_possible
        model = init_hmm_model(motif=motif, p=p, alpha=float(command_args.alpha))
        log_likelihoods = model.TCM_EM(seqs, float(command_args.convergenceThr))
        ll_histories[i] = (motif, log_likelihoods)
        output_motif_txt(model, fasta_seqs, i+1)
    output_histories(ll_histories, motifs)

    # --------------- Plot the required graphs ---------------
    if PLOT_MODE:
        import matplotlib.pyplot as plt
        ll = ll_histories[0][1]
        iterations = np.arange(1, len(ll) + 1, 1)
        plt.title("EM log likelihood as a function of #iterations")
        plt.plot(iterations, ll)
        plt.show()



if __name__ == '__main__':
    main()