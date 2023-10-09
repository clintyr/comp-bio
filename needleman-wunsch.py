base_idx = { 'A' : 0, 'G' : 1, 'C' : 2, 'T' : 3 } # link sequence to similarity matrix
S = [
	# A  G   C   T
	[3, -1, -2, -2], # A
	[-1, 3, -2, -2], # G
	[-2, -2, 3, -1], # C
	[-2, -2, -1, 3]  # T
	] # similarity matrix for scoring alignment
gap_pen = 4 # gap penalty
PTR_NONE, PTR_GAP1, PTR_GAP2, PTR_BASE = 0, 1, 2, 3 # useful for traceback

def seqalign_DP(seq1, seq2, sim_matrix, gap_pen):
	"""returns the score of the optimal Needleman-Wunsch alignment for seq1 and seq2
	with some alignment similarity matrix and gap_pen > 0; 
	also returns score matrix and DP table for traceback 
	"""
	F = [[0 for j in range(len(seq2) + 1)] for i in range(len(seq1) + 1)]
	TB = [[PTR_NONE for j in range(len(seq2) + 1)] for i in range(len(seq1) + 1)]

	# initialize dynamic programming table for Needleman-Wunsch alignment
	for i in range(1, len(seq1)+1):
		F[i][0], TB[i][0] = -i*gap_pen, PTR_GAP2 # For gaps in the seq2 (columns)
	
	for j in range(1, len(seq2)+1):
		F[0][j], TB[0][j] = -j*gap_pen, PTR_GAP1 # indicates a gap in seq1 (rows)
    
    # populate the table while moving through sequences
	for i in range(1, len(seq1) + 1):
		for j in range(1, len(seq2) + 1):
			res = [float('-inf'), F[i][j-1] - gap_pen, F[i-1][j] - gap_pen, F[i-1][j-1]
			       + sim_matrix[base_idx[seq1[i-1]]][base_idx[seq2[j-1]]]]
			opt_score = max(res)
			opt_index = res.index(opt_score)
			F[i][j], TB[i][j] = opt_score, opt_index

	return F[len(seq1)][len(seq2)], F, TB


def traceback(seq1, seq2, TB):
	"""traces back the optimal Needleman-Wunsch alignment for seq1, seq2 from DP matrix"""
	s1, s2 = "", ""
	i, j = len(seq1), len(seq2)
	while TB[i][j] != PTR_NONE:
		if TB[i][j] == PTR_BASE: # base to base alignment
			s1 = seq1[i-1] + s1; s2 = seq2[j-1] + s2
			i -= 1; j-= 1
			
		elif TB[i][j] == PTR_GAP1: # gap in seq1
			s1 = '-' + s1 ; s2 = seq2[j-1] + s2
			j -= 1
			
		elif TB[i][j] == PTR_GAP2: # gap in seq2
			s1 = seq1[i-1] + s1; s2 = '-' + s2
			i -= 1
			
		else: assert False

	return s1, s2

def getAlign(seq1,seq2):
	score, F, TB = seqalign_DP(seq1, seq2, S, gap_pen)
	s1, s2 = traceback(seq1, seq2, TB)
	return score, F, TB, s1, s2


seq1, seq2 = "CTAAGTACT", "CATTA"
score, F, TB, s1, s2 = getAlign(seq1, seq2)

# print("\n".join([str(score), s1, s2]))
# for line in F:
# 	print(line)