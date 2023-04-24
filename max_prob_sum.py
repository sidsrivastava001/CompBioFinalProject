"""Enhanced Nussinov algorithm for RNA secondary structure prediction,
    as described in the paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6540740/ 

Args:
    probL[i]: The probability of the ith base of the RNA sequence being "("
    probR[i]: The probability of the ith base of the RNA sequence being ")" 
    probP[i]: The probability of the ith base of the RNA sequence being "."
    seq: the RNA sequence for which the secondary structure is to be predicted

    Note that the probL, probR, probP list values come form the results of the 
    CNN model, which tells us the likelihood of each base being a "(" or ")" or "."
    That is, the model will tell us the probability of each RNA base pairing up with
    another base, or not pairing up with another base. Using this information, we need to run a version of the Nussinov algorithm in order to find a valid dot bracket 
representation that maximizes the the cumulative probabilities of each base pair.

Returns:
    The maximum probability sum of the given sequence. 
"""


# Helper function to compute the score of a base pair (Ri, Rj)
def computePairScores(Rj, Ri, i, j, probL, probR, probP):
  # set of possible rna base pairs
  possibleBasePairs = set(['AU', 'UA', 'GC', 'CG', 'GU', 'UG'])

  if Ri + Rj in possibleBasePairs:
    return probL[i] + probR[j]
  else:
    return probP[i] + probP[j]


def MaximumProbabilitySum(probL, probR, probP, seq):
  n = len(seq)
  N = [[0 for i in range(n)] for j in range(n)]
  for i in range(n):
    for j in range(n):
      if i == j:
        N[i][j] = probP[i]
      elif j < i:
        N[i][j] = 0
      else:
        N[i][j] = max(
          N[i + 1][j] + probP[i],
          N[i][j - 1] + probP[j], 
          N[i + 1][j - 1] + computePairScores(seq[i], seq[j], i, j, probL, probR, probP),
          max([N[i][k] + N[k + 1][j] for k in range(i, j)]))
  return N

# Backtracking function to find the optimal dot bracket representation of the RNA sequence
# from the DP table
def BackTrack(i, j, fold, N, seq):
  if i >= j:
    return
  elif N[i][j] == N[i + 1][j] + probP[i]:
    BackTrack(i + 1, j, fold, N, seq)
  elif N[i][j] == N[i][j - 1] + probP[j]:
    BackTrack(i, j - 1, fold, N, seq)
  elif N[i][j] == (N[i + 1][j - 1] + computePairScores(seq[i], seq[j], i, j, probL, probR, probP)):
    fold[i] = '('
    fold[j] = ')'
    BackTrack(i + 1, j - 1, fold, N, seq)
  else:
    for k in range(i + 1, j - 1):
      if N[i][j] == N[i][k] + N[k + 1][j]:
        BackTrack(i, k, fold, N, seq)
        BackTrack(k + 1, j, fold, N, seq)
        break

  return fold

# Main function to run the Enhanced Nussinov algorithm
def EnhancedNussinov(probL, probR, probP, seq):
  n = len(seq)
  N = MaximumProbabilitySum(probL, probR, probP, seq)
  fold = ['.' for i in range(n)]
  fold = BackTrack(0, n - 1, fold, N, seq)
  return ''.join(fold)


    
