"""
Python script to test the accuracy of nussinov algorithm
against the RNA sequences downloaded from bprna database 
"""

import os
import sys
import numpy as np
from nussinov import Nussinov


def read_dbn(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        seq = lines[3].strip()
        dbn = lines[4].strip()
    return seq, dbn


def compare_dbn(dbn1, dbn2):
    if len(dbn1) != len(dbn2):
        return -1 

    n = len(dbn1)
    correct = 0

    for i in range(n):
        if dbn1[i] == dbn2[i]:
            correct += 1

    return correct


def run_tests():
    test_dir = "/Users/krishsuraparaju/Downloads/dbnFiles" # Change this to the directory where the dbn files are stored

    files = sorted(os.listdir(test_dir))

    num_passing = 0
    num_tests_ran = 0
    threshold_score = 0 

    for filename in files:
        if num_tests_ran == 1000:
            break

        if filename.endswith('.dbn'):
            seq, dbn = read_dbn(os.path.join(test_dir, filename))
            n = len(seq)
            
            # Skip over any sequence that has 'N' (which indicates sequencing wasn't accurate)
            valid = True if (set(seq) == {'A', 'C', 'G', 'U'}) else False 

            # Only consider sequences with length <= 120 
            if(n > 120 or not valid): 
                continue

            num_tests_ran+=1

            pred_dbn = Nussinov(seq)

            threshold_score = 0.6*n # At least 60% of the pairs must be correct (arb value)

            correct = compare_dbn(dbn, pred_dbn)
            seq_score = (correct/n)

            print(f'{filename}: {seq_score}')
            print(f'Threshold score: {threshold_score}')

            if seq_score >= threshold_score:
                num_passing += 1
    
    with(open('test_nussinov_res.txt')) as f:
        f.write(f'{num_passing} out of {num_tests_ran} correct\n')
        f.write(f'Accuracy: {num_passing/num_tests_ran}\n')


if __name__ == '__main__':
    run_tests()
            