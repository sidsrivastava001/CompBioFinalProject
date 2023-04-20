import numpy as np
import pandas as pd
import os
import math

pairs = {'G': 'C', 'C': 'G', 'A': 'U', 'U': 'A'}
wobble = {'G': 'U', 'U': 'G'}

def initialize_matrix(seq, dist, n):
  M = np.empty((n, n))
  M[:] = np.NaN

  # for the condition that j + 5 > i -> inf
  for x in range(M.shape[0]):
      for y in range(M.shape[1]):
        if y + dist > x:
          M[x, y] = np.inf
  return M

# h(x, y) = 2(i - j + 5)
def calc_hairpin_energy(x, y):
  return (x - y + 5) * 2

# stem energy: initialized to -4 (pairs), 0 (wobble), 4 (otherwise)
def calc_stem_energy(seq, x, y):
  if (pairs[seq[x]] == seq[y]):
    return -4
  elif ((seq[x] == "G" or seq[x] == "U") and wobble[seq[x]] == seq[y]):
    return 0
  else:
    return 4

# calculate energy for structures
# hairpins: V(x, y) = s(x, y) + h(x - 1, y + 1)
# match: V(x, y) = s(x, y) + W(x - 1, y + 1)
def calc_V(x, y, V, W, seq):
  hairpin_val = calc_hairpin_energy(x, y) + calc_stem_energy(seq, x, y)
  match_val = calc_stem_energy(seq, x, y) + W[x - 1, y + 1]
  if (hairpin_val < match_val):
    V[x, y] = hairpin_val
    return ('H', V)
  else:
    V[x, y] = match_val
    return ('M', V)

# calculate energy of optimal structure
def calc_W(x, y, V, W, path, seq):
  structure, V = calc_V(x, y, V, W, seq)

  k_dict = {}
  for k in range(y + 2, x):
    if k not in k_dict:
      k_dict[k] = W[x, k] + W[k - 1, y]
  min_k = min(k_dict.keys(), key=(lambda k: k_dict[k]))

  # update matrix W with the minimum choice
  W[x, y] = min(W[x - 1, y], W[x, y + 1], V[x, y], k_dict[min_k])
  if (W[x, y] == W[x - 1, y]):
    path[x, y] = 'L'
  if (W[x, y] == W[x, y + 1]):
    path[x, y] = 'D'
  if (W[x, y] == V[x, y]):
    path[x, y] = structure
  if (W[x, y] == k_dict[min_k]):
    path[x, y] = str(min_k)

  return (V, W, path)

def Zuker(V, W, n, dist, path, seq):
  final = W[n - 1, 0]

  while math.isnan(final):
    for x in range(n):
      for y in range(n):
        if y + dist <= x:
          left = W[x - 1, y]
          down = W[x, y + 1]
          if math.isnan(down) or math.isnan(left):
            continue
          V, W, path = calc_W(x, y, V, W, path, seq)
      final = W[n - 1, 0]

  return (V, W, path)

def Backtrack(path, backtrack, n):
  match = []

  curr_x = 0
  curr_y = n - 1
  curr_cell = path[curr_x, curr_y]
  path = path.transpose()
  backtrack[curr_x, curr_y] = "#"

  while curr_cell != 'H':
    if curr_cell == 'M':
      match.append((curr_x, curr_y))
      curr_x += 1
      curr_y -= 1
    elif curr_cell == "L":
      curr_y -= 1
    elif curr_cell == "D":
      curr_x += 1
    curr_cell = path[curr_x, curr_y]
    backtrack[curr_x, curr_y] = "#"

  hairpin = (curr_x, curr_y)
  return (match, hairpin, path, backtrack)

def dot_bracket(path, backtrack, n):
    matches, hairpin, path, backtrack = Backtrack(path, backtrack, n)
  
    dot_bracket = ["." for _ in range(n)]
    for match in matches:
      first_paren = match[0]
      second_paren = match[1]
      dot_bracket[first_paren] = '('
      dot_bracket[second_paren] = ')'
    dot_bracket = "".join(dot_bracket)

    return dot_bracket

seq = "AAUACUCCGUUGCAGCAU"
distance_constant = 5
n = len(seq)
W = initialize_matrix(seq, distance_constant, n)
V = initialize_matrix(seq, distance_constant, n)
path = np.zeros((n, n), 'U1')
backtrack = np.zeros((n, n), 'U1')

V, W, path = Zuker(V, W, n, distance_constant, path, seq)
db = dot_bracket(path, backtrack, n)
print(seq)
print(db)
