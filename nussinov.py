import numpy as np
import sys

# Constants
min_loop_length = 0
pairs = {"A": "U", "U": "A", "G": "C", "C": "G"}

# Checks if pair is complementary
def isPair(pair : str) -> bool:
  return pairs[pair[0]] == pair[1]

# Initializes an n x n matrix
def initMatrix(Text : str) -> list[list[int]]:
  n = len(Text)

  M = np.empty([n, n])
  M[:] = np.NAN

  M[range(n), range(n)] = 0
  M[range(1, n), range(n - 1)] = 0

  return M

# Forward direction for Nussinov algorithm
def Forward(M : list[list[int]], Text : str) -> list[list[int]]:
  n = len(Text)

  for k in range(1, n):
    for i in range(n - k):
      j = i + k

      if (j - i >= min_loop_length):
        down = M[i + 1][j]
        left = M[i][j - 1]
        diag = M[i + 1][j - 1] + isPair((Text[i], Text[j]))
        other = max([M[i][t] + M[t + 1][j] for t in range(i, j)])

        M[i][j] = max(down, left, diag, other)
      else:
        M[i][j] = 0

  return M

# Ming direction of Nussinov algorithm
def Back(M : list[list[int]], Text : str, pairs : list[tuple[int, int]], i : int, length : int) -> list[tuple[int, int]]: 
  j = length

  if i < j:
    if M[i][j] == M[i + 1][j]:
      Back(M, Text, pairs, i + 1, j)
    elif M[i][j] == M[i][j - 1]:
      Back(M, Text, pairs, i, j - 1)
    elif M[i][j] == M[i + 1][j - 1] + isPair((Text[i], Text[j])):
      pairs.append((i, j))
      Back(M, Text, pairs, i + 1, j - 1)
    else:
      for k in range(i + 1, j - 1):
        if M[i][j] == M[i, k] + M[k + 1][j]:
          Back(M, Text, pairs, i, k)
          Back(M, Text, pairs, k + 1, j)
          break


# Convert to dotBracket
def dotBracket(Text : str, pairs : list[tuple[int, int]]) -> str:
  dot = ["." for i in range(len(Text))]

  for s in pairs:
    dot[max(s)] = ")"
    dot[min(s)] = "("

  return "".join(dot)

# Run Nussinov algorithm
def Nussinov(Text : str) -> str:
  n = len(Text)

  initial_matrix = initMatrix(Text)
  forward_matrix = Forward(initial_matrix, Text)

  pairs = []
  Back(forward_matrix, Text, pairs, 0, n - 1)

  return dotBracket(Text, pairs)

if __name__ == "__main__":
  open("output.txt", "w").write(Nussinov(sys.argv[1]))

