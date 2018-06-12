#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 12:35:54 2018

@author: caitlin
"""

## Practice at creating a split decomposition from a matrix

import numpy as np
np.set_printoptions(threshold=np.inf)  # Print whole matrix
import itertools

# Enter distance matrix and taxa labels for practicing
d = np.matrix([[0,0,0,0,0,0,0],[4.0,0,0,0,0,0,0],[5.0,1.0,0,0,0,0,0],[7.0,3.0,2.0,0,0,0,0],
               [13.0,9.0,8.0,6.0,0,0,0],[8.0,12.0,13.0,11.0,5.0,0,0],[6.0,10.0,11.0,13.0,7.0,2.0,0]])
taxa = ["A","B","C","D","E","F","G"]

# Function - get all options for n choose 2 from a given set
def choose2(collection):
    #print("in choose2")
    combinations = list(itertools.combinations_with_replacement(collection,2))
    return(combinations)

# Function - calculate a pairwise distance
def pd(a,b,names,matrix):
    #print("in pd")
    a_ind = names.index(a)
    b_ind = names.index(b)
    a_pd = matrix[a_ind,b_ind] # Check to see whether the pairwise distance is stored under
    b_pd = matrix[b_ind,a_ind] # d[A,B] or d[B,A]
    dist = max(a_pd,b_pd) # If both entries are 0, will have a pairwise distance of 0
    return(dist)

# Function - calculate the inequality for a single quartet
def inequality(i,j,k,l,names,matrix):
    #print("in inequality")
    # Calculate d_ij + d_kl < max{d_ik + d_jl, d_il + d_jk}
    d1 = pd(i,j,names,matrix) + pd(k,l,names,matrix)
    d2 = pd(i,k,names,matrix) + pd(j,l,names,matrix)
    d3 = pd(i,l,names,matrix) + pd(j,k,names,matrix)
    #print(d1,d2,d3)
    if d1 < max(d2,d3):
        #result = True
        result = 0
        #print("min")
    elif (d1 > max(d2,d3)):
        #result = False
        result = 1
        #print("max")
    elif (d1 >= max(d2,d3)):
        #print("equal")
        # result = False
        result = 1
    return(result)
    
# Function - for one subset, identify whether it's a d-split or not
def issplit(partition,names,matrix):
    #print("in issplit")
    # Create a list to store 0 if the inequality holds and 1 if it doesn't
    its = [] # store one number for each iteration
    # Find all combinations for the first subset
    combs1 = choose2(partition[0])
    # Find all combinations for the second subset
    combs2 = choose2(partition[1])
    # Iterate through first subset combinations
    for x in combs1:
        #print(x)
        # Iterate through second subset combinations
        for y in combs2:
            #print(y)
            # Calculate inequality for those 4 taxa
            i = x[0]
            j = x[1]
            k = y[0]
            l = y[1]
            calculate_inequality = inequality(i,j,k,l,names,matrix)
            its.append(calculate_inequality)
    summed = sum(its)
    # If inequality holds, sum will be 0
    if summed == 0:
        split_result = True
    # If doesn't hold, sum larger than 0
    if summed > 0:
        split_result = False
    return(split_result)

# Function - Separate into partitions of two disjoint subsets
def partition2(collection):
    #print("in partition")
    parts = []
    # Need to go through the number of elements in each subset
    # n choose 1
    # n choose 2
    # ....
    # up to n choose (n/2) - for 6 need n choose 3, for 7 need n choose 3
    n = len(collection)
    k = range(1,n/2+1) # go from 1 to n/2
    t = set(collection) # se
    for i in k:
        sub = list(itertools.combinations(collection,i))
        for x in sub:
            x = set(x)
            y = (t-x) # Find disjoint subset 
            temp = [list(x),list(y)]
            parts.append(temp)
    return(parts)

def isolation_index(partition,names,matrix):
    #print("in isolation_index")
    # a = 1/2 min(max{d_ij+d_kl,d_ik+d_jl,d_il+d_jk}-d_ij-d_kl)
    mins = []
    # Find all combinations for the first subset
    combs1 = choose2(partition[0])
    # Find all combinations for the second subset
    combs2 = choose2(partition[1])
    # Try all quartets
    for x in combs1:
        for y in combs2:
            i = x[0]
            j = x[1]
            k = y[0]
            l = y[1]
            # Calculate all distances
            d1 = pd(i,j,names,matrix) + pd(k,l,names,matrix) # d_ij + d_kl
            d2 = pd(i,k,names,matrix) + pd(j,l,names,matrix) # d_ik + d_jl
            d3 = pd(i,l,names,matrix) + pd(j,k,names,matrix) # d_il + d_jk
            temp = max(d1,d2,d3) - d1
            mins.append(temp)
    minimum = min(mins)
    alpha = 0.5 * minimum
    return(alpha)

# Create a matrix based on the partition
def make_splitmatrix(partition,names,alpha):
    #print("in make_splitmatrix")
    # For a tree AB-CD, the distance AB and CD should 0, and AC=BC=AD=BC=alpha (alpha is the isolation index)
    # The matrix would then be [[0 0 a a], [ 0,0,a,a], [a,a,0,0], [a,a,0,0]]
    # Create a matrix of the same number of rows and columns as the original matrix
    matrix = make_matrix(names)
    pd = choose2(names)
    ss1 = partition[0]
    ss2 = partition[1]
    # For each possible combination of 2 taxa:
    for i in pd:
    # Check if they're in the same subset
        if (i[0] in ss1 and i[1] in ss1) or (i[0] in ss2 and i[1] in ss2):
            # If they are, leave that pairwise distance as a 0
            continue
        if (i[0] in ss1 and i[1] in ss2) or (i[0] in ss2 and i[1] in ss1):
            # If they're not, change that distance to alpha (the isolation index) 
            matrix[names.index(i[1]),names.index(i[0])] = alpha
    return(matrix)
                   
def make_matrix(names):
    #print("in make_matrix")
    row = [0] * len(names)
    rows = [row] * len(names)
    matrix = np.matrix(rows)
    return(matrix)
    

# Calculate all partitions
sets = partition2(taxa)
summed_matrix = make_matrix(taxa)
for part in sets:
    print(part)
    split = issplit(part,taxa,d)
    if split == True:
        # calculate isolation index
        ii = isolation_index(part,taxa,d)
        # create matrix 
        matrix = make_splitmatrix(part,taxa,ii)
        # print found a split, partition, matrix
        print("Found a split:", part)
        print(matrix)
        summed_matrix = summed_matrix + matrix
    if split == False:
        continue
print("final distance matrix:")
print(summed_matrix)
        


