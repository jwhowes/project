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
			if adj_matrix[i][j] == 1 and x[i] == x[j]:
				num += 1
	return num

beta = 1.5
alpha = 1

sigma_p = 	((gamma(1+beta)*sin(pi*beta/2))/
			(gamma((1+beta)/2)*beta*2**((beta-1)/2)))**(2/beta)  # The variance of p

def levy():
	"""Generates a number with a levy distribution"""
	p = np.random.normal(0, sigma_p)
	q = np.random.normal(0, 1)
	return p/(abs(q)**(1/beta))

def sigmoid(x):
	return 1/(1+exp(-x))

def generate_cuckoo():
	colouring = np.zeros((n, k), dtype=int)
	ncl = 0
	Y = np.array([i for i in range(n)])
	uncoloured = np.array([i for i in range(n)])
	while ncl < k:
		Y = uncoloured
		while len(Y) > 0:
			v = np.random.choice(Y)
			colouring[v][ncl] = 1
			uncoloured = np.delete(uncoloured, np.where(uncoloured == v))
			Y = np.delete(Y, [x for x in range(len(Y)) if Y[x] == v or adj_matrix[v][Y[x]] == 1])
		ncl += 1
	# Colour remaining nodes randomly
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
	S = np.vectorize(lambda x: sigmoid(x + alpha*levy()))(colouring)
	ret = np.zeros((n, k), dtype=int)
	for i in range(n):
		for j in range(k):
			if S[i][j] > np.random.uniform(0, 1) and any_neighbour_coloured(i, j, colouring):  # This can colour the same vertex more than one colour.
				ret[i][j] = 1
	# Colour any uncoloured nodes uniformly at random
	for i in range(n):  # Could possibly make more efficient
		if 1 not in ret[i]:
			ret[i][np.random.randint(0, k)] = 1
	return ret

pa = 0.25
num_nests = 50
num_iterations = 1000

nests = np.array([generate_cuckoo() for i in range(num_nests)])