# Time for 6 nodes: 57.87
# Time for 10 nodes: 1:07.21 = 67.21
# Ratio of time: 1.16
# Ratio of nodes: 1.67 so it's looking sublinear (thank god)
# Also it coloured the petersen graph in 3 colours

# Time for 100 nodes: 5.06 minutes (with num_eggs = 1)
# Estimated time for 1000 nodes: > 20 hrs (with num_eggs = 1)

import numpy as np
import time

edge_list = [	[1, 4, 5],
				[0, 2, 6],
				[1, 3, 7],
				[2, 4, 8],
				[3, 1, 9],
				[0, 7, 8],
				[1, 8, 9],
				[2, 5, 8],
				[3, 5, 6],
				[4, 6, 7]]

def make_graph(num_vertices, edge_probability):
	global edge_list
	edge_list = np.empty((num_vertices), dtype=object)
	adj_matrix = np.zeros((num_vertices, num_vertices), dtype=int)
	for i in range(num_vertices):
		edge_list[i] = []
	for i in range(num_vertices):
		for j in range(i):
			if np.random.uniform(0, 1) < edge_probability:
				edge_list[i].append(j)
				edge_list[j].append(i)
make_graph(100, 0.5)
#print("100 vertices, edge probability 0.5")
n = len(edge_list)

def f(x):
	"""Returns the number of colours used by a colouring"""
	return x.max() + 1

#def f(x):  # The fitness function proposed by local search colouring survey
#	classes = np.zeros(x.max() + 1, dtype=int)
#	ret = 0
#	for i in range(n):
#		classes[x[i]] += 1
#	for c in classes:
#		ret -= c*c
#	return ret

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

def populate_dist_matrix():
	for i in range(n_pop):
		for j in range(i):
			dist_matrix[i][j] = d(cuckoos['cuckoo'][i], cuckoos['cuckoo'][j])
			dist_matrix[j][i] = dist_matrix[i][j]

k = 3
def goal_point():  # Could try precomputing d_bar to self for each cluster? (at each iteration)
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
	cluster_fitness = [[cuckoos['fitness'][c] for c in clusters[i]] for i in range(k)]
	best_cluster = -1
	best_cluster_mean_fitness = 0
	for i in range(k):
		if len(cluster_fitness[i]) > 0 and (best_cluster == -1 or np.mean(cluster_fitness[i]) < best_cluster_mean_fitness):
			best_cluster = i
			best_cluster_mean_fitness = np.mean(cluster_fitness[i])
	gp = -1
	for i in range(n_pop):
		if S[i] == best_cluster and (gp == -1 or cuckoos['fitness'][i] < cuckoos['fitness'][gp]):
			gp = i
	return gp

def valid(v, c, col):
	"""Returns whether or not vertex v can be coloured c in col"""
	for i in edge_list[v]:
		if col[i] == c:
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
	"""Migrates cuckoo x towards cuckoo y (in place)"""
	# Generate adaptive migration coefficient
	if f_max == f_min:
		F = F_max
	else:
		F = F_min + ((cuckoos['fitness'][x] - f_min)/(f_max - f_min))*(F_max - F_min)
	# Generate random proportion of distance to travel
	r = np.random.uniform(0, 1)
	# Let I be the vertices on which x and y disagree
	I = np.array([i for i in range(n) if cuckoos['cuckoo'][x][i] != cuckoos['cuckoo'][y][i]])
	for i in range(int(F * r * len(I))):
		v = I[i]
		# Replace i with y's colouring for i
		cuckoos['cuckoo'][x][v] = cuckoos['cuckoo'][y][v]
		# Clean up col
		for j in range(i + 1, n):
			# For every conflicting edge (v,j), replace j's colour with the smallest legal colour different than y[j]
			if cuckoos['cuckoo'][x][j] == cuckoos['cuckoo'][x][v] and j in edge_list[v]:
				c = 0
				while True:
					if c != cuckoos['cuckoo'][y][j] and valid(j, c, cuckoos['cuckoo'][x]):
						cuckoos['cuckoo'][x][j] = c
						break
					c += 1

alpha_max = 20
alpha_min = 0

F_max = 1  # In the continuous domain F_max and F_min can take any value (F_max > F_min)
F_min = 0  # In this domain it only makes sense that F_max, F_min come from (0, 1)

n_pop = 5
n_max = 50
num_iterations = 3000
p = 0.1

start = time.time()

cuckoos = np.array([(generate_cuckoo(), 0) for i in range(n_pop)], dtype=[('cuckoo', np.ndarray), ('fitness', int)])
cuckoos['fitness'] = np.array([f(c) for c in cuckoos['cuckoo']])

dist_matrix = np.zeros((n_max, n_max), dtype=int)

for t in range(num_iterations):
	num_eggs = np.random.randint(5, 21, (n_pop))
	tot_eggs = num_eggs.sum()
	eggs = np.array([(np.zeros(n, dtype=int), 0) for i in range(tot_eggs)], dtype=[('cuckoo', np.ndarray), ('fitness', int)])
	egg = 0
	alpha = alpha_max - ((alpha_max - alpha_min)/(num_iterations - t))  # I don't think I need the +1 in the denominator as t < num_iterations always
	for i in range(n_pop):
		elr = alpha * num_eggs[i]/tot_eggs * n  # I've just put n instead of (v_hi - v_lo)
		for j in range(num_eggs[i]):
			eggs['cuckoo'][egg] = get_egg(cuckoos['cuckoo'][i], elr)
			eggs['fitness'][egg] = f(eggs['cuckoo'][egg])
			egg += 1
	eggs = eggs[eggs['fitness'].argsort()]
	cuckoos = np.append(eggs[:int(tot_eggs * (1 - p))], cuckoos, axis=0)
	if len(cuckoos) > n_max:
		cuckoos = cuckoos[cuckoos['fitness'].argsort()]
		cuckoos.resize((n_max))
	n_pop = len(cuckoos)
	f_max = cuckoos['fitness'][cuckoos['fitness'].argmax()]
	f_min = cuckoos['fitness'][cuckoos['fitness'].argmin()]
	gp = goal_point()
	for i in range(n_pop):
		migrate(i, gp)

cuckoos = cuckoos[cuckoos['fitness'].argsort()]
print(cuckoos[0])
print("Time taken:", time.time() - start)