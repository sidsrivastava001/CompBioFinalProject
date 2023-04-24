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
# def computePairScores(Rj, Ri, i, j, probL, probR, probP):
#   # set of possible rna base pairs
#   possibleBasePairs = set(['AU', 'UA', 'GC', 'CG', 'GU', 'UG'])

#   if Ri + Rj in possibleBasePairs:
#     return probL[i] + probR[j]
#   else:
#     return probP[i] + probP[j]


# def MaximumProbabilitySum(probL, probR, probP, seq):
#   n = len(seq)
#   N = [[0 for i in range(n)] for j in range(n)]
#   for i in range(n):
#     for j in range(n):
#       if i == j:
#         N[i][j] = probP[i]
#       elif j < i:
#         N[i][j] = 0
#       else:
#         N[i][j] = max(
#           N[i + 1][j] + probP[i],
#           N[i][j - 1] + probP[j], 
#           N[i + 1][j - 1] + computePairScores(seq[i], seq[j], i, j, probL, probR, probP),
#           max([N[i][k] + N[k + 1][j] for k in range(i, j)]))
#   return N

# # Backtracking function to find the optimal dot bracket representation of the RNA sequence
# # from the DP table
# def BackTrack(i, j, fold, N, seq):
#   if j <= i:
#     return
#   elif N[i][j] == N[i + 1][j] + probP[i]:
#     BackTrack(i + 1, j, fold, N, seq)
#   elif N[i][j] == N[i][j - 1] + probP[j]:
#     BackTrack(i, j - 1, fold, N, seq)
#   elif N[i][j] == N[i + 1][j - 1] + computePairScores(seq[i], seq[j], i, j, probL, probR, probP):
#     fold[i] = '('
#     fold[j] = ')'
#     BackTrack(i + 1, j - 1, fold, N, seq)
#   else:
#     for k in range(i, j):
#       if N[i][j] == N[i][k] + N[k + 1][j]:
#         BackTrack(i, k, fold, N, seq)
#         BackTrack(k + 1, j, fold, N, seq)
#         break
  

# # Main function to run the Enhanced Nussinov algorithm
# def EnhancedNussinov(probL, probR, probP, seq):
#   n = len(seq)
#   N = MaximumProbabilitySum(probL, probR, probP, seq)
#   fold = ['.' for i in range(n)]
#   BackTrack(0, n-1, fold, N, seq)
#   return fold

import numpy as np
import sys



# Checks if pair is complementary and assigns the appropriate score as defined in the paper
def isPairScore(seq, i, j, probL, probR, probP):
  possibleBasePairs = set(['AU', 'UA', 'GC', 'CG', 'GU', 'UG']) # includes wobble pairs
  Ri = seq[i]
  Rj = seq[j]

  if Ri + Rj in possibleBasePairs:
    return probL[i] + probR[j]
  else:
    return probP[i] + probP[j]


# Initializes an n x n matrix
def initMatrix(Text : str) -> list[list[int]]:
  n = len(Text)

  M = np.empty([n, n])
  M[:] = np.NAN

  M[range(n), range(n)] = 0
  M[range(1, n), range(n - 1)] = 0

  return M

# Forward direction for Updated Nussinov algorithm
def Forward(probL, probR, probP, M, Text, min_loop_length = 0):
  n = len(Text)

  for k in range(1, n):
    for i in range(n - k):
      j = i + k

      if (j - i >= min_loop_length):
        down = M[i + 1][j]
        left = M[i][j - 1]
        diag = M[i + 1][j - 1] + isPairScore(Text, i, j, probL, probR, probP)
        other = max([M[i][t] + M[t + 1][j] for t in range(i, j)])

        M[i][j] = max(down, left, diag, other)
      else:
        M[i][j] = 0

  return M

# Back direction of Updated Nussinov algorithm, backtrack
def Back(M, Text, pairs, i, length, probL, probR, probP):
  j = length

  if i < j:
    if M[i][j] == M[i + 1][j]:
      Back(M, Text, pairs, i + 1, j, probL, probR, probP)
    elif M[i][j] == M[i][j - 1]:
      Back(M, Text, pairs, i, j - 1, probL, probR, probP)
    elif M[i][j] == M[i + 1][j - 1] + isPairScore(Text, i, j, probL, probR, probP):
      pairs.append((i, j))
      Back(M, Text, pairs, i + 1, j - 1, probL, probR, probP)
    else:
      for k in range(i + 1, j - 1):
        if M[i][j] == M[i, k] + M[k + 1][j]:
          Back(M, Text, pairs, i, k, probL, probR, probP)
          Back(M, Text, pairs, k + 1, j, probL, probR, probP)
          break


# Convert to dotBracket
def dotBracket(Text : str, pairs : list[tuple[int, int]]) -> str:
  dot = ["." for i in range(len(Text))]

  for s in pairs:
    dot[max(s)] = ")"
    dot[min(s)] = "("

  return "".join(dot)

# Run Updated Nussinov algorithm
def UpdatedNussinov(Text, probL, probR, probP):
  n = len(Text)

  initial_matrix = initMatrix(Text)
  forward_matrix = Forward(probL, probR, probP, initial_matrix, Text)

  pairs = []
  Back(forward_matrix, Text, pairs, 0, n - 1, probL, probR, probP)

  return dotBracket(Text, pairs)

# probL = [0.5, 0.5, 0.5, 0.5, 0.5]
# probR = [0.5, 0.5, 0.5, 0.5, 0.5]
# probP = [0.5, 0.5, 0.5, 0.5, 0.5]
# seq = "GGGGG"
# print(UpdatedNussinov(seq, probL, probR, probP))
