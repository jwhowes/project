import numpy as np

adj_matrix = [	[0, 1, 0, 0, 1, 1, 0, 0, 0, 0],
				[1, 0, 1, 0, 0, 0, 1, 0, 0, 0],
				[0, 1, 0, 1, 0, 0, 0, 1, 0, 0],
				[0, 0, 1, 0, 1, 0, 0, 0, 1, 0],
				[1, 0, 0, 1, 0, 0, 0, 0, 0, 1],
				[1, 0, 0, 0, 0, 0, 0, 1, 1, 0],
				[0, 1, 0, 0, 0, 0, 0, 0, 1, 1],
				[0, 0, 1, 0, 0, 1, 0, 0, 0, 1],
				[0, 0, 0, 1, 0, 1, 1, 0, 0, 0],
				[0, 0, 0, 0, 1, 0, 1, 1, 0, 0]]
n = len(adj_matrix)
k = 3

def gpx(p1, p2):
	"""Produces offspring of p1 and p2 using Greedy Partition Crossover"""
	a = p1.copy()
	b = p2.copy()
	ret = [[] for i in range(k)]
	for i in range(k):
		if i % 2 == 0:
			parent = a
		else:
			parent = b
		m = 1
		for j in range(1, k):
			if len(parent[j]) > len(parent[m]):
				m = j
		ret[i] = parent[m]
		# Remove the vertices of parent[m] from a and b
	# Colour uncoloured nodes in ret randomly

def generate_member():  # Currently members are created randomly
	assignments = np.random.randint(0, k, (n))
	col = [[] for i in range(k)]
	for i in range(n):
		col[assignments[i]].append(i)
	return col

pop_size = 100

pop = np.array([generate_member() for i in range(pop_size)], dtype=object)
