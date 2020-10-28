import numpy as np
from math import sin, gamma, pi, exp

adj_matrix = [  [0, 1, 1, 0, 1, 0, 0],
				[1, 0, 0, 0, 1, 0, 1],
				[1, 0, 0, 0, 0, 1, 1],
				[0, 0, 0, 0, 0, 0, 1],
				[1, 1, 0, 0, 0, 0, 0],
				[0, 0, 1, 0, 0, 0, 0],
				[0, 1, 1, 1, 0, 0, 0]]
k = 3
n = len(adj_matrix)

def f(x):
	"""Returns the number of conflicts in a given colouring x"""
	num = 0
	for i in range(len(adj_matrix)):
		for j in range(i):
			if adj_matrix[i][j] == 1 and colour(x[i]) == colour(x[j]):
				num += 1
	return num

beta = 1.5
alpha = 1

sigma_p = 	((gamma(1+beta)*sin(pi*beta/2))/
			(gamma((1+beta)/2)*beta*2**((beta-1)/2)))**(2/beta)  # The variance of p

def colour(v):
	"""Converts from binary colour vector to single colour"""
	return np.where(v == 1)[0][0]

def levy():
	"""Generates a number with a levy distribution"""
	p = np.random.normal(0, sigma_p)
	q = np.random.normal(0, 1)
	return p/(abs(q)**(1/beta))

def sigmoid(x):
	"""Applies the sigmoid fucntion to a scalar x"""
	try:
		return 1/(1+exp(-x))
	except OverflowError:  # If |x| is too big exp(x) returns an overflow error
		if x < 0:
			return 0
		return 1

def generate_cuckoo():
	"""Generates a random (potentially invalid) colouring for the graph"""
	colouring = np.zeros((n, k), dtype=int)
	ncl = 0  # The current colour being applied
	uncoloured = np.array([i for i in range(n)])  # An array of all uncoloured vertices
	while ncl < k:
		# Colour as many as possible ncl before moving on
		Y = uncoloured
		while len(Y) > 0:
			# Choose a random vertex v and colour it ncl
			v = np.random.choice(Y)
			colouring[v][ncl] = 1
			# Remove v from uncoloured and remove v and v's neighbourhood from Y
			uncoloured = np.delete(uncoloured, np.where(uncoloured == v))
			Y = np.delete(Y, [x for x in range(len(Y)) if Y[x] == v or adj_matrix[v][Y[x]] == 1])
		ncl += 1
	# Colour remaining uncoloured nodes randomly
	for v in uncoloured:
		colouring[v][np.random.randint(0, k)] = 1
	return colouring

def any_neighbour_coloured(v, c, colouring):
	"""Returns whether or not any neighbour of v is coloured c in colouring"""
	for i in range(n):
		if adj_matrix[v][i] == 1 and colouring[i][c] == 1:
			return True
	return False

def levy_flight(colouring):
	"""Returns the result of performing a levy flight on colouring"""
	# Apply the sigmoid function to every entry of colouring
	S = np.vectorize(lambda x: sigmoid(x + alpha*levy()))(colouring)
	ret = np.zeros((n, k), dtype=int)
	for i in range(n):
		for j in range(k):
			# Vertex i is coloured colour j with probability S[i][j] and any neighbour of i is coloured j in colouring
			if S[i][j] > np.random.uniform(0, 1) and any_neighbour_coloured(i, j, colouring):  # This can colour the same vertex more than one colour (but for some reason it never returns one of those)
				ret[i][j] = 1
	# Colour any uncoloured nodes uniformly at random
	for i in range(n):  # Could possibly make more efficient
		if 1 not in ret[i]:
			ret[i][np.random.randint(0, k)] = 1
	return ret

pa = 0.25
num_nests = 50
num_iterations = 100

nests = np.array([generate_cuckoo() for i in range(num_nests)])

for t in range(num_iterations):
	for i in range(num_nests):
		j = levy_flight(nests[i])
		if f(j) < f(nests[i]):
			nests[i] = j
	nests = nests[np.array([f(n) for n in nests]).argsort()]
	for i in range(int(num_nests * pa)):
		nests[num_nests - i - 1] = generate_cuckoo()

nest = nests[np.array([f(n) for n in nests]).argmin()]
print(nest, f(nest))