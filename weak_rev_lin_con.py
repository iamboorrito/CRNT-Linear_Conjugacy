'''
Created on Mar 13, 2018

Based on the paper 'A Linear Programming Approach to Weak Reversibility 
and Linear Conjugacy of Chemical Reaction Networks' by Matthew D. Johnston, 
David Siegel and G ́abor Szederkenyi (2011) arXiv:1107.1659

@author: Evan Burton
'''

from crnpy.crn import from_react_strings
from pulp import LpProblem, LpVariable, lpSum, value
from itertools import product
import numpy as np
from pulp.solvers import CPLEX
from pulp.constants import LpMaximize

def weak_rev_lin_con(crn, eps, ubound):
    
    m = crn.n_complexes
    n = crn.n_species
    
    # Y is n by m
    Y = np.array(crn.complex_matrix).astype(np.float64)
    
    # Ak is m by m
    Ak = np.array(crn.kinetic_matrix).astype(np.float64)
    
    # M is n by m
    M = Y.dot(Ak)

    print("Input CRN defined by\nY =\n", Y)
    print("\nM =\n", M)    
    print("\nAk =\n", Ak)
    print("\nComputing linearly conjugate network with min deficiency ... START\n")

    # Ranges for iteration
    col_range = range(1, m+1)
    row_range = range(1, n+1)

    # Set up problem model object
    prob = LpProblem("Weakly Reversible Linearly Conjugate Network", LpMaximize)
    
    # Get a list of all off-diagonal entries to use later to greatly
    # simplify loops
    off_diag = [(i, j) for i,j in product(col_range, repeat=2) if i != j]
    
    # Decision variables for matrix A, only need off-diagonal since A
    # has zero-sum columns
    A = LpVariable.dicts("A", [(i, j) for i in col_range for j in col_range])
    Ah = LpVariable.dicts("Ah", [(i, j) for i in col_range for j in col_range])
     
    # Decision variables for the diagonal of T
    T = LpVariable.dicts("T", row_range, eps, ubound)
    
    # Binary variables for counting partitions used and assigning complexes to linkage classes
    delta = LpVariable.dicts("delta", off_diag, 0, 1, "Integer")
    
    # Objective
    prob += -lpSum(delta[i, j] for (i,j) in off_diag)
    
    # Y*A = T*M
    for i in row_range:
        for j in col_range:
            prob += lpSum( Y[i-1, k-1]*A[k, j] for k in col_range ) == M[i-1, j-1]*T[i]
            
    # A and Ah have zero-sum columns
    for j in col_range:
        prob += lpSum( A[i,j] for i in col_range ) == 0    
        prob += lpSum( Ah[i,j] for i in col_range ) == 0
        prob += lpSum( Ah[j,i] for i in col_range ) == 0
    
    # Off-diagonal entries are nonnegative and are switch on/off by delta[i,j]
    for (i,j) in off_diag:
        # A constraints
        prob += A[i, j] >= 0
        prob += A[i, j] - eps*delta[i, j] >= 0
        prob += A[i, j] - ubound*delta[i, j] <= 0
        
        # Ah constraints
        prob += Ah[i, j] >= 0
        prob += Ah[i, j] - eps*delta[i, j] >= 0
        prob += Ah[i, j] - ubound*delta[i, j] <= 0
        
    # Diagonal entries of A, Ah are non-positive
    for j in col_range:
        prob += A[j, j] <= 0
        prob += Ah[j, j] <= 0
        
    status = prob.solve(solver=CPLEX()) 
    #print(prob)
    
    # Problem successfully solved, report results
    if status == 1:
        
        # Get solutions to problem
        Tsol = np.zeros((n, n))
        Asol = np.zeros((m,m))
        
        for i in col_range:
            for j in col_range:
                Asol[i-1, j-1] = value(A[i, j])
        
        for i in row_range:
            Tsol[i-1, i-1] = value(T[i])
        
        print("\nA =\n", Asol)
        print("\nT =\n", Tsol)

    else:
        print("No solution found")
        
crn = from_react_strings(["X1 + 2 X2 ->(1.5) X1", 
                          "2 X1 + X2 ->(1) 3 X2",
                          "X1 + 3 X2 ->(1) X1 + X2",
                          "X1 + X2 ->(1) 3 X1 + X2"])

weak_rev_lin_con(crn, 0.6666, 20)