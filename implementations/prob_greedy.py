import numpy as np

adj_matrix = [	[0, 1, 1, 0, 1, 1, 0, 0, 0, 0],
				[1, 0, 1, 0, 0, 0, 1, 0, 0, 0],
				[1, 1, 0, 1, 0, 0, 0, 1, 0, 0],
				[0, 0, 1, 0, 1, 0, 0, 0, 1, 0],
				[1, 0, 0, 1, 0, 0, 0, 0, 0, 1],
				[1, 0, 0, 0, 0, 0, 0, 1, 1, 0],
				[0, 1, 0, 0, 0, 0, 0, 0, 1, 1],
				[0, 0, 1, 0, 0, 1, 0, 0, 0, 1],
				[0, 0, 0, 1, 0, 1, 1, 0, 0, 0],
				[0, 0, 0, 0, 1, 0, 1, 1, 0, 0]]
num_vertices = len(adj_matrix)

degree = np.zeros((num_vertices), dtype=int)

colouring = np.zeros((num_vertices), dtype=int)

for i in range(num_vertices):
	degree[i] = sum(adj_matrix[i])

order = np.array(sorted([i for i in range(num_vertices)], key=lambda x : -degree[x]))

for v in order:
	probs = np.zeros((num_vertices), dtype=int)
	# Colour vertex v based on prob_matrix[v]