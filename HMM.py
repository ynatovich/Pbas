import numpy as np

# -------------------- Constants --------------------
MAX_EM_ITERS = 100
MOTIF_STATE = 'M'
BACKGROUND_STATE = 'B'
LETTERS = ['A', 'T', 'C', 'G']
BACKGROUND_SYMBOL = '0'
MOTIF_SYMBOL = '1'

class HMMState:
    def __init__(self, emission_probabilities, type=None):
        """
        :param emission_probabilities: A dictionary of the form
        {l: e(l)} for each l=possible emitted letter and e(l) the
        emission probability
        """
        self.emission_probabilities = emission_probabilities
        self.log_emission_probabilities = {key: np.log(p) for key, p in
                                           emission_probabilities.items()}
        self.type = type

    def update_log_emissions(self, letter, new_log_prob):
        self.log_emission_probabilities[letter] = new_log_prob
        self.emission_probabilities[letter] = np.exp(new_log_prob)

class HMM_Model:
    def __init__(self, states):
        self.states = states
        self.start_state = self.states[0]
        self.transition_table = None
        self.log_transition_table = None

    def set_transition_table(self, table):
        self.transition_table = table
        self.log_transition_table = np.log(table)

    def update_log_transition_table(self, i, j, value):
        self.transition_table[i, j] = np.exp(value)
        self.log_transition_table[i, j] = value

    def forward(self, seq):
        """
        Implement the forward algorithm for the given HMM model and sequence
        :return: The forward table (dim = |model states| X |sequence length| )
         and the computed log likelihood
        """
        # -------------------- Initialization --------------------
        L = len(seq)
        n_states = len(self.states)
        f = np.zeros(shape=(len(self.states), L), dtype=np.float64)
        f[0, 0] = 1
        # -------------------- Recursion --------------------
        for i in range(1, L):
            for l, state in enumerate(self.states):
                e_l_xi = state.log_emission_probabilities[seq[i]]
                f[l, i] = e_l_xi + np.logaddexp.reduce(
                    ([f[k, i - 1] + self.log_transition_table[k, l]
                      for k in range(n_states)]))
        # -------------------- Termination --------------------
        p = f[0, -1]
        # p = np.logaddexp.reduce([f[k, L - 1] + self.log_transition_table[k, 0]
        #                          for k in range(n_states)])
        return f, p

    def backwards(self, seq):
        """
        Implement the backwards algorithm for the given HMM model and sequence
        :return: The backwards table (dim = |model states| X |sequence length|)
         and the computed log likelihood
        """
        # -------------------- Initialization --------------------
        L = len(seq)
        n_states = len(self.states)
        b = np.zeros(shape=(len(self.states), L), dtype=np.float64)
        for k, state in enumerate(self.states):
            b[k, L - 1] = self.transition_table[k, 0]
        # -------------------- Recursion --------------------
        for i in reversed(range(L - 1)):
            for k, state in enumerate(self.states):
                arr = [0] * n_states
                for l, state in enumerate(self.states):
                    e_l_xi = state.log_emission_probabilities[seq[i + 1]]
                    arr[l] = e_l_xi + self.log_transition_table[k, l] + b[l, i + 1]
                b[k, i] = np.logaddexp.reduce(arr)
        # -------------------- Termination --------------------
        arr = [0] * n_states
        for l, state in enumerate(self.states):
            arr[l] = \
                self.log_transition_table[0, l] + \
                state.log_emission_probabilities[seq[0]] + b[l, 0]
        p = np.logaddexp.reduce(arr)
        return b, p

    def MAP_path(self, seq):
        """
        Maximum a-posteriori path. Computed for each observation by taking
        pi[i] = argmax_k p(pi[i]=k | seq)
        where
        p(pi[i]=k | seq) = (f[k, i] * b[k, i]) / p(seq)
        f,b are the forward/backward matrices and p(seq) is obtained during the
        forward run
        :return: The computed MAP-path and the forward log likelihood
        """
        f, pf = self.forward(seq)
        b, pb = self.backwards(seq)
        L = len(seq)
        posteriors = np.zeros(shape=(L, len(self.states)), dtype=np.float64)
        for i in range(1, len(seq)):
            for k, state in enumerate(self.states):
                posteriors[i, k] = (f[k, i] * b[k, i]) / pf
        path = [0] * L
        for i in range(L):
            path[i] = np.argmax([posteriors[i, k] for k in range(len(self.states))])
        return path, pf

    def viterbi(self, seq):
        """
        Implement the Viterbi algorithm for the given HMM model and sequence
        :return: The computed path and log likelihood
        """
        # -------------------- Initialization --------------------
        L = len(seq)
        n_states = len(self.states)
        v = np.zeros(shape=(len(self.states), L), dtype=np.float64)
        ptr = np.zeros(shape=(len(self.states), L), dtype=int)
        v[0, 0] = 1
        max_p_k = arg_max_p_k = 0
        # -------------------- Recursion --------------------
        for i in range(1, len(seq)):
            for l, state in enumerate(self.states):
                e_l_xi = state.log_emission_probabilities[seq[i]]
                arr = ([v[k, i - 1] + self.log_transition_table[k, l]
                        for k in range(n_states)])
                max_p_k = np.max(arr)
                arg_max_p_k = np.argmax(arr)
                v[l, i] = e_l_xi + max_p_k
                ptr[l, i] = arg_max_p_k
        # -------------------- Traceback --------------------
        path = self.traceback(ptr, arg_max_p_k)
        return path, max_p_k

    @staticmethod
    def traceback(ptr_matrix, last_state):
        """
        Given a pointer matrix and the last state :return: the path
        """
        n_states, L = ptr_matrix.shape
        path = [0] * L
        path[-1] = last_state
        ptr = ptr_matrix[last_state, -1]
        for i in reversed(range(L - 1)):
            path[i] = ptr
            ptr = ptr_matrix[ptr, i]
        return path



    def TCM_EM(self, seqs, convergenceThr):
        """
        Perform the EM algorithm to find plausible parameters for this model given the TCM
        architecture described in 76558 ex2
        :param seqs: An array of sequences on which to perform the algorithm
        :return:
        """
        iters = 1
        log_likelihoods = [-np.inf]  # The log likelihoods after each iteration
        forward_matrices = [np.zeros(shape=(len(self.states), len(seq)), dtype=np.float64) for seq in seqs]
        backward_matrices = [np.zeros(shape=(len(self.states), len(seq)), dtype=np.float64) for seq in seqs]
        seq_log_likelihoods = [0] * len(seqs)
        # -------------------- EM-loop --------------------
        while True:
            print(iters)
        # -------------------- E-step --------------------
            # -------------------- f/b matrices--------------------
            for j, seq_j in enumerate(seqs):
                f, pf = self.forward(seq_j)
                b, pb = self.backwards(seq_j)
                forward_matrices[j] = f
                backward_matrices[j] = b
                seq_log_likelihoods[j] = pf
            iteration_log_likelihood = np.logaddexp.reduce(seq_log_likelihoods)
            # -------------------- transition --------------------
            E_transitions = []
            for transfer in [[0, 0], [0, 1]]:
                arr = [np.zeros(len(seq_j) - 1, dtype=np.float64) for seq_j in seqs]
                k, l = transfer
                for j, seq_j in enumerate(seqs):
                    for i in range(len(seq_j) - 1):
                        arr[j][i] = \
                            (forward_matrices[j][k, i] +
                             self.log_transition_table[k, l] +
                             self.states[l].log_emission_probabilities[seq_j[i + 1]] +
                             backward_matrices[j][l, i + 1] -
                             seq_log_likelihoods[j])
                arr = [x for s in arr for x in s]
                E_transitions.append(np.logaddexp.reduce(arr))
            # -------------------- emissions --------------------
            E_emissions = [dict() for i in range(len(self.states))]
            for k, state in enumerate(self.states):
                if state.type == BACKGROUND_STATE:
                    continue
                for x in LETTERS:
                    arr = [[-np.inf]*(len(seq_j)) for seq_j in seqs]
                    for j, seq_j in enumerate(seqs):
                        for i in range(len(seq_j)):
                            if seq_j[i] == x:
                                arr[j][i] = \
                                    (forward_matrices[j][k, i] +
                                     backward_matrices[j][k, i] -
                                     seq_log_likelihoods[j])
                    arr = [x for s in arr for x in s]
                    E_emissions[k][x] = np.logaddexp.reduce(arr)
        # -------------------- M-step --------------------
            # -------------------- transition --------------------
            transition_log_sum = np.logaddexp.reduce(E_transitions)
            self.update_log_transition_table(0, 0, E_transitions[0] - transition_log_sum)
            self.update_log_transition_table(0, 1, E_transitions[1] - transition_log_sum)
            # -------------------- emissions --------------------
            for k, state in enumerate(self.states):
                if state.type == BACKGROUND_STATE:
                    continue
                em_log_sum = np.logaddexp.reduce([E_emissions[k][y] for y in LETTERS])
                for x in LETTERS:
                    state.update_log_emissions(letter=x, new_log_prob=E_emissions[k][x] - em_log_sum)

        # --------------- Convergence test ---------------
            if iters > MAX_EM_ITERS or abs(iteration_log_likelihood - log_likelihoods[-1]) < convergenceThr:
                break
            # if iteration_log_likelihood > log_likelihoods[-1]:
            log_likelihoods.append(iteration_log_likelihood)
            iters += 1
        return log_likelihoods[1:]




    def __str__(self):
        s = "----------------------------------\n"
        s += "HMM_Model\n"
        s += "Transition matrix:\n"
        s += str(self.transition_table)
        s += "\nStates:\n"
        for i, state in enumerate(self.states):
            s +="state_" + str(i) + "emissions: " + str(state.emission_probabilities) + "\n"
        s += "----------------------------------\n"
        return s










