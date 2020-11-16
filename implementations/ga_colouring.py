# My own formulation (loosely based on 1998 paper)
# TODO:
	# Implement parent selection (roullette approach)
	# Implement fitness function for partition variant (currently it's for assignment variant)
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

def f(x):
	"""Returns the number of conflicts in a given colouring x"""
	num = 0
	for i in range(len(adj_matrix)):
		for j in range(i):
			if adj_matrix[i][j] == 1 and x[i] == x[j]:
				num += 1
	return num

def gpx(p1, p2):
	"""Produces offspring of p1 and p2 using Greedy Partition Crossover"""
	a = np.array([x[:] for x in p1], dtype=object)
	b = np.array([x[:] for x in p2], dtype=object)
	ret = [[] for i in range(k)]
	uncoloured = np.array([i for i in range(n)])
	for i in range(k):
		if i % 2 == 0:
			parent = a; non_parent = b
		else:
			parent = b; non_parent = a
		m = 0
		for j in range(1, k):
			if len(parent[j]) > len(parent[m]):
				m = j
		ret[i] = parent[m]
		uncoloured = np.setdiff1d(uncoloured, parent[m])
		parent[m] = []
		for j in ret[i]:
			for l in range(k):
				if j in non_parent[l]:
					non_parent[l].remove(j)
					break
	for i in uncoloured:
		ret[np.random.randint(0, k)].append(i)
	return ret

def mutate(x):
	"""Randomly moves a vertex from one colour class to another (in place)"""
	c = np.random.choice([i for i in range(k) if len(x[i]) > 0])
	d = np.random.randint(0, k)
	v = np.random.choice(x[c])
	x[c].remove(v)
	x[d].append(v)

def generate_member():  # Currently members are created randomly
	assignments = np.random.randint(0, k, (n))
	col = [[] for i in range(k)]
	for i in range(n):
		col[assignments[i]].append(i)
	return col

pop_size = 50
num_iterations = 100
mutation_prob = 0.1

pop = np.array([generate_member() for i in range(pop_size)], dtype=object)
for t in range(num_iterations):
	newP = np.empty((pop_size,), dtype=object)
	for i in range(pop_size):
		newP[i] = [[] for j in range(k)]
	for i in range(pop_size):
		# Select p1 and p2 with roulette wheel approach
		p1 = pop[0]
		p2 = pop[1]
		newP[i] = gpx(p1, p2)
		if np.random.uniform(0, 1) < mutation_prob:
			mutate(newP[i])
	P = newP

P = P[np.array([f(p) for p in P]).argsort()]
print(P[0])