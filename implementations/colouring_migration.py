import numpy as np

adj_matrix = [  [0, 1, 1, 0, 1, 0, 0],
				[1, 0, 0, 0, 1, 0, 1],
				[1, 0, 0, 0, 0, 1, 1],
				[0, 0, 0, 0, 0, 0, 1],
				[1, 1, 0, 0, 0, 0, 0],
				[0, 0, 1, 0, 0, 0, 0],
				[0, 1, 1, 1, 0, 0, 0]]
n = len(adj_matrix)

def d(x, y):
	ret = 0
	for i in range(n):
		if x[i] != y[i]:
			ret += 1
	return ret

def valid(v, c, col):
	"""Returns whether or not vertex v can be coloured c in col"""
	for i in range(n):
		if adj_matrix[v][i] == 1 and col[i] == c:
			return False
	return True

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
			# For every conflicting edge (i,j), replace j's colour with the smallest legal colour different than y[j]
			if adj_matrix[j][v] == 1 and x[j] == x[v]:
				for c in range(n):
					if c != y[j] and valid(j, c, x):
						x[j] = c
						break
	minimise(x)
	return x

def get_egg(x, elr):
	"""Returns a valid colouring within distance elr from x"""
	num = np.random.randint(0, elr + 1)
	# colour num vertices greedily
	col = x.copy()
	for i in range(num):
		v = np.random.randint(0, n)
		for c in range(n):
			if c != col[v] and valid(v, c, col):
				col[v] = c
				break
	return col

x = np.array([1, 2, 2, 2, 0, 1, 1])
y = np.array([0, 1, 2, 1, 2, 0, 0])