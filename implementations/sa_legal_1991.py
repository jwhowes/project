# This is slightly modified from the 1991 paper (I use colourings as assignments, rather than partitions)
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

def f(x):  # The fitness function proposed by local search colouring survey
	classes = np.zeros(x.max() + 1, dtype=int)
	ret = 0
	for i in range(n):
		classes[x[i]] += 1
	for c in classes:
		ret -= c*c
	return ret

def num_colours(x):
	return x.max() + 1

def valid(v, c, col):
	"""Returns whether or not vertex v can be coloured c in col"""
	for i in range(n):
		if adj_matrix[v][i] == 1 and col[i] == c:
			return False
	return True

def initial_solution():
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

def kempe_chain(c, d, v):
	"""Returns a (c,d)-kempe in s including a vertex v"""
	K = [v]
	i = 0
	while i < len(K):
		for j in range(n):
			if j not in K and adj_matrix[K[i]][j] == 1 and (s[j] == c or s[j] == d):
				K.append(j)
		i += 1
	return K


def get_neighbour():
	neighbour = s.copy()
	v = np.random.randint(0, n)
	c = s[v]
	d = np.random.choice(np.setdiff1d(list(range(num_colours(s))), [c]))
	K = kempe_chain(c, d, v)
	for i in K:
		if s[i] == c:
			neighbour[i] = d
		else:
			neighbour[i] = c
	return neighbour

freeze_lim = 5
size_factor = 16
cutoff = 0.1
min_percent = 0.2

T = 10000
beta = 1.0005

epsilon = 0.0001

s = initial_solution()
c = f(s)
best = s
best_c = c

freeze_count = 0
while freeze_count < freeze_lim and T > epsilon:  # I added the epsilon check (the algorithm wasn't terminating otherwise)
	changes = 0; trials = 0
	best_changed = False
	while trials < size_factor*num_colours(s)*n and changes < cutoff*num_colours(s)*n:
		trials += 1
		neighbour = get_neighbour()
		neighbour_fitness = f(neighbour)
		d = neighbour_fitness - c
		if c < 0:
			changes += 1
			s = neighbour
			c = neighbour_fitness
			if c < best_c:
				best_changed = True
				best_s = s
				best_c = c
		elif np.random.uniform(0, 1) <= np.exp(-d/T):
			changes += 1
			if c < best_c:
				best_changed = True
				best_s = s
				best_c = c
	T /= beta
	if best_changed:
		freeze_count = 0
	if changes/trials < min_percent:
		freeze_count += 1

print(best)