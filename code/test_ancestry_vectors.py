"""Test ancestry vectors module with simple 3-pop admixture case."""

import sys
sys.path.insert(0, "/home/yanlin/1KGothers/code")

import networkx as nx
from hapgraph.topology.greedy_search import AdmixGraph
from hapgraph.topology.ancestry_vectors import (
    compute_ancestry_vectors,
    expected_f2_matrix,
    expected_f3_array,
    expected_f4_array,
)

# Simple 3-population tree: root -> (A, B, C)
# Admixture: C gets 0.3 from A
tree = nx.DiGraph()
tree.add_edge("root", "A", weight=0.1)
tree.add_edge("root", "B", weight=0.1)
tree.add_edge("root", "C", weight=0.05)
tree.add_edge("root", "anc", weight=0.0)  # internal node

graph = AdmixGraph(
    tree=tree,
    admixture_edges=[("A", "C", 0.3)],
    populations=["A", "B", "C"],
)

Q, pops, branches = compute_ancestry_vectors(graph)
branch_lengths = [tree[u][v]["weight"] for u, v in branches]

print("Populations:", pops)
print("Branches:", branches)
print("Branch lengths:", branch_lengths)
print("\nAncestry matrix Q (rows=populations, cols=branches):")
print(Q)
print("\nRow sums (should be 1):", Q.sum(axis=1))

f2_mat = expected_f2_matrix(Q, branch_lengths)
print("\nF2 matrix:")
print(f2_mat)
print("\nF2(A,B) =", f2_mat[0, 1])
print("F2(A,C) =", f2_mat[0, 2])
print("F2(B,C) =", f2_mat[1, 2])
print("Expected: F2(A,B) > F2(A,C) since C is admixed from A ->", f2_mat[0, 1] > f2_mat[0, 2])

# F3 tests: F3(C; A, B) - C is admixed from A, so F3 should be negative
triples = [(2, 0, 1)]  # F3(C; A, B)
f3_arr = expected_f3_array(Q, branch_lengths, triples)
print("\nF3(C; A, B) =", f3_arr[0], "(negative indicates admixture)")

# F4 test
quartets = [(0, 1, 2, 2)]  # F4(A, B; C, C) = 0 (same pop)
quartets2 = [(0, 2, 1, 2)]  # F4(A, C; B, C)
f4_arr = expected_f4_array(Q, branch_lengths, quartets2)
print("F4(A, C; B, C) =", f4_arr[0])
