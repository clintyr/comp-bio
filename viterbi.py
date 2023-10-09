# VITERBI Algorithm for HMM
# + refers to High-GC, and - refers to Low-GC in this example
# When indexing 0 is + and 1 is -.

from math import log
state_idx = { '+' : 0, '-' : 1 }

# initial distribution over states
init_dist = [0.5,0.5]

# transition probabilities
tr = [
    #  to+   to-
    [ 0.99, 0.01 ], # from+
    [ 0.01, 0.99 ]  # from-
]

# emission probabilities
em = [
    #    A     G     C     T
    [ 0.20, 0.30, 0.30, 0.20], # +
    [ 0.30, 0.20, 0.20, 0.30]  # -
]

def viterbi(X):
    """Returns the Viterbi path for the emission sequence X.
    X should be a list of integers, 0=A, 1=G, 2=C, 3=T.
    The returned Y is a list of integers, 0=High-GC, 1=Low-GC.
    """
    N = len(tr)
    L = len(X)
    assert len(em) == N
    V = [[0]*N for _ in range(L)]
    TB = [[0]*N for _ in range(L)]
    for i in range(L):
        Vprev = []
        if i == 0:
            Vprev = [log(pk0) for pk0 in init_dist]
        else:
            Vprev = V[i-1]
        for k in range(N):
            max_tr_prob = Vprev[0] + log(tr[0][k])
            last_state = 0
            for prev_state in range(N)[1:]:
                tr_prob = Vprev[prev_state] + log(tr[prev_state][k])
                if tr_prob > max_tr_prob:
                    max_tr_prob = tr_prob
                last_state = prev_state
                max_prob = max_tr_prob + log(em[k][X[i]])
            V[i][k] = max_prob
            TB[i][k] = last_state

    # perform traceback and return the predicted hidden state sequence
    Y = [-1 for i in range(L)]
    _, yL = max([(V[L - 1][k], k) for k in range(N)])
    Y[L - 1] = yL
    for i in range(L - 2, -1, -1):
        Y[i] = TB[i + 1][Y[i + 1]]
    return Y
