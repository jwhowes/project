import numpy as np

adj_matrix = [  [0, 1, 1, 1, 1, 1],
				[1, 0, 1, 1, 1, 1],
				[1, 1, 0, 1, 1, 1],
				[1, 1, 1, 0, 1, 1],
				[1, 1, 1, 1, 0, 1],
				[1, 1, 1, 1, 1, 0]]
n = len(adj_matrix)

def f(x):
	"""Returns the number of colours used by a colouring"""
	return x.max() + 1

def d(x, y):
	"""Returns the distance between x and cuckoo y"""
	ret = 0
	for i in range(n):
		if x[i] != y[i]:
			ret += 1
	return ret

def d_bar_sum(S1, S2):
	num = 0
	for i in S1:
		for j in S2:
			num += dist_matrix[i][j]
	return num

def tri_dist(i, j, S, dbss, clusters):
	"""Finds the triangular distance between cuckoo i and cluster j"""
	if len(clusters[j]) == 0:
		return -1
	return 2*d_bar_sum([i], clusters[j])/len(clusters[j]) - dbss[j]/(len(clusters[j]) * len(clusters[j]))

k = 3
def goal_point():
	populate_dist_matrix()
	S = np.random.randint(0, k, (n_pop))
	clusters = [[] for i in range(k)]
	for i in range(n_pop):
		clusters[S[i]].append(i)
	db_sum_self = np.zeros((k), dtype=int)  # Stores the sum of distance from every member to every other member (for calculation of d_bar(i, i))
	for i in range(k):
		for j in range(len(clusters[i])):
			for l in range(j):
				db_sum_self[i] += 2*dist_matrix[clusters[i][j]][clusters[i][l]]
	converged = False
	while not converged:
		converged = True
		for i in range(n_pop):
			cluster = -1
			tdc = -1
			for j in range(k):
				td = tri_dist(i, j, S, db_sum_self, clusters)
				if (cluster == -1 or td < tdc) and td != -1:
					cluster = j
					tdc = td
			if cluster != S[i]:
				clusters[S[i]].remove(i)
				clusters[cluster].append(i)
				converged = False
				for j in clusters[S[i]]:
					db_sum_self[S[i]] -= 2*dist_matrix[j][i]
				for j in clusters[cluster]:
					db_sum_self[cluster] += 2*dist_matrix[i][j]
				S[i] = cluster
	cluster_fitness = [[f(cuckoos[c]) for c in clusters[i]] for i in range(k)]
	best_cluster = -1
	best_cluster_mean_fitness = 0
	for i in range(k):
		if len(cluster_fitness[i]) > 0 and (best_cluster == -1 or np.mean(cluster_fitness[i]) < best_cluster_mean_fitness):
			best_cluster = i
			best_cluster_mean_fitness = np.mean(cluster_fitness[i])
	gp = -1
	for i in range(n_pop):
		if S[i] == best_cluster and (gp == -1 or f(cuckoos[i]) < f(cuckoos[gp])):
			gp = i
	return gp

def valid(v, c, col):
	"""Returns whether or not vertex v can be coloured c in col"""
	for i in range(n):
		if adj_matrix[v][i] == 1 and col[i] == c:
			return False
	return True

def generate_cuckoo():
	colouring = np.array([-1 for i in range(n)])
	order = np.array([i for i in range(n)])
	np.random.shuffle(order)
	for v in order:
		c = 0
		while True:
			if valid(v, c, colouring):
				colouring[v] = c
				break
			c += 1
	return colouring

def get_egg(x, elr):
	"""Returns a valid colouring within distance elr from x"""
	num = np.random.randint(0, elr + 1)
	# colour num vertices greedily
	col = x.copy()
	for i in range(num):
		v = np.random.randint(0, n)
		for c in range(n):
			if c != x[v] and valid(v, c, col):
				col[v] = c
				break
	return col

def minimise(col):
	"""Reduces the max number of a colouring (while keeping all colours consistent)"""
	colours = list(np.unique(col, axis=0))
	for i in range(n):
		col[i] = colours.index(col[i])

def migrate(x, y):
	"""Migrates x towards y (in place)"""
	# Generate random proportion of distance to travel
	r = np.random.uniform(0, 1)
	# Let I be the vertices on which x and y disagree
	I = np.array([i for i in range(n) if x[i] != y[i]])
	for i in range(int(r * len(I))):
		v = I[i]
		# Replace i with y's colouring for i
		x[v] = y[v]
		# Clean up col
		for j in range(i + 1, n):
			# For every conflicting edge (v,j), replace j's colour with the smallest legal colour different than y[j]
			if adj_matrix[j][v] == 1 and x[j] == x[v]:
				c = 0
				while True:
					if c != y[j] and valid(j, c, x):
						x[j] = c
						break
					c += 1
	#minimise(x)  # Not sure if necessary

def populate_dist_matrix():
	for i in range(n_pop):
		for j in range(i):
			dist_matrix[i][j] = d(cuckoos[i], cuckoos[j])
			dist_matrix[j][i] = dist_matrix[i][j]

alpha = 10
n_pop = 5
n_max = 50
num_iterations = 3000
p = 0.1

cuckoos = np.array([generate_cuckoo() for i in range(n_pop)])

dist_matrix = np.zeros((n_max, n_max), dtype=int)
populate_dist_matrix()

for t in range(num_iterations):
	num_eggs = np.random.randint(5, 21, (n_pop))
	tot_eggs = num_eggs.sum()
	eggs = np.zeros((tot_eggs, n), dtype=int)
	egg = 0
	for i in range(n_pop):
		elr = alpha * num_eggs[i]/tot_eggs * n  # I've just put n instead of (v_hi - v_lo)
		for j in range(num_eggs[i]):
			eggs[egg] = get_egg(cuckoos[i], elr)
			egg += 1
	eggs = eggs[np.array([f(e) for e in eggs]).argsort()]
	cuckoos = np.append(eggs[:int(tot_eggs * (1 - p))], cuckoos, axis=0)
	if len(cuckoos) > n_max:
		cuckoos = cuckoos[np.array([f(c) for c in cuckoos]).argsort()]
		cuckoos.resize((n_max, n))
	n_pop = len(cuckoos)
	gp = goal_point()
	for i in range(n_pop):
		migrate(cuckoos[i], cuckoos[gp])

cuckoos = cuckoos[np.array([f(c) for c in cuckoos]).argsort()]
print(cuckoos[0])