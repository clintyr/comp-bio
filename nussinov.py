				 #  A   U   C   G
scoring_matrix = [[ 0, -1,  0,  0], # A
				  [-1,  0,  0, -1], # U
				  [ 0,  0,  0, -1], # C
				  [ 0, -1, -1,  0]] # G

# the above matrix is scoring A-U, G-U, and C-G pairs as -1 and all other pairs as 0.

b_map = {'A': 0, 'U': 1, 'C': 2, 'G': 3} # base map

def nussinov(seq):
	"""returns the score of nussinov algorithm given an RNA sequence"""
    
	seq = [b for b in seq]
	N = len(seq)
	Nussi = [[0 for j in range(N)] for i in range(N+1)]
	for i in range(N-1, 0, -1):
		for j in range(N-1, 0, -1):
			k = j - (N-i)
			if k < 0 or j < 0:
				continue
			best = k+1
			# print(best)
			res1 = Nussi[best][j]
			res2 = Nussi[best][best-1] + Nussi[best+1][j] + scoring_matrix[b_map[seq[best-1]]][b_map[seq[best]]]
			for l in range(k+1, j+1):
				res = Nussi[k+1][l-1] + Nussi[l+1][j] + scoring_matrix[b_map[seq[k]]][b_map[seq[l]]]
				res2 = min(res, res2)
                
			Nussi[k][j] = min(res1, res2)
   
	score = Nussi[0][N-1]

	return score

seq = 'GCACGACGAA'
print(nussinov(seq))