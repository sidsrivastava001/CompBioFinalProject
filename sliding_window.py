import numpy as np

sliding_param = 19
n = 20

def create_rand_matrix(n):
    return [[(i + j) for i in range(n)] for j in range(n)]
    
M = create_rand_matrix(n)
numpy_M = np.array(M, dtype=int)
zeros = np.zeros((sliding_param, n), dtype=int)
numpy_M = np.append(numpy_M, zeros, axis=0)
sliding_mats = [numpy_M[i: i + sliding_param, :].tolist() for i in range(n)]
