Scientific Problem: Inferring secondary structures of unknown RNA molecules with an ML and dynamic programming based prediction algorithm.

RNA secondary structure prediction is a vital stepping stone towards tertiary prediction, providing researchers with an understanding of protein binding mechanisms and relationships with other proteins. Through identification of hydrogen bond patterns, we discover important structures that govern tertiary binding, like pseudoknots, hairpins, and internal/external loops. Finding a more efficient mechanism for calculating these secondary structures can lead to exciting discoveries in the RNA field.

Approach and Feasibility:

We construct an algorithm with an RNA string as input and secondary structure as output:
We first convert the string to an RNA matrix representation based on a hydrogen bond recurrence relation that will generate a neighbor matrix between all the nucleotides.
We’ll generate several normalized matrices using a sliding window.
These values will be passed into a CNN that will output a dot bracket structure that shows which nucleotides base pair with each other.
We run a maximum probability sum model that uses dynamic programming to turn this matrix into a secondary structure (outlined in the paper).

We’ll compare this approach with the Nussinov algorithm using the data referenced below. This is feasible because it’s already been implemented and the CNN structure is easily implementable using existing libraries.

Data and External Resources: 

For this project, we will be referencing the paper “A New Method of RNA Secondary Structure Prediction Based on Convolutional Neural Network and Dynamic Programming” by Zhang et. al. to create our convolutional neural network model architecture as well as the dynamic programming preprocessing of the RNA datasets. 

We will also use Nussinov’s algorithm for a baseline comparison of structures. We will use the bpRNA-1m dataset that contains known RNA dot-bracket representations to train our model. Finally, we will run our algorithms on a new set of roughly 200 RNA sequences identified in the paper “Promiscuous splicing-derived hairpins are dominant substrates of tailing-mediated defense of miRNA biogenesis in mammals” that have not been experimentally verified. Finally, we will be using the Keras ML framework to construct our model.
