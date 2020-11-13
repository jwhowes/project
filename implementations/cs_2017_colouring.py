# Running time on 100 vertices: 7.72 mins
import numpy as np
from math import sin, gamma, pi
import time

adj_matrix = [  [0, 1, 1, 0, 1, 0, 0],
				[1, 0, 0, 0, 1, 0, 1],
				[1, 0, 0, 0, 0, 1, 1],
				[0, 0, 0, 0, 0, 0, 1],
				[1, 1, 0, 0, 0, 0, 0],
				[0, 0, 1, 0, 0, 0, 0],
				[0, 1, 1, 1, 0, 0, 0]]
k = 20
n = len(adj_matrix)

def make_graph(num_vertices, edge_probability):
	global adj_matrix
	adj_matrix = np.zeros((num_vertices, num_vertices), dtype=int)
	for i in range(num_vertices):
		for j in range(i):
			if np.random.uniform(0, 1) < edge_probability:
				adj_matrix[i][j] = 1
				adj_matrix[j][i] = 1
n = len(adj_matrix)

def f(x):
	"""Returns the number of conflicts in a given colouring x"""
	num = 0
	for i in range(n):
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

def levy_flight(colouring):
	"""Returns the result of performing a levy flight on a colouring"""
	# Generate M
	M = int(alpha*levy()) + 1
	ret = colouring.copy()
	# Randomly change M colours in ret
	for i in range(M):
		j = np.random.randint(0, n)
		ret[j] = np.random.choice([i for i in range(k) if i != ret[j]])
	return ret

pa = 0.25
num_nests = 50
num_iterations = 3000
parasitism_comparison = True

#nests = np.random.randint(0, k, (num_nests, n))
nests = np.array([(np.random.randint(0, k, (n)), 0) for i in range(num_nests)], dtype=[('nest', np.ndarray), ('fitness', int)])
nests['fitness'] = np.array([f(n) for n in nests['nest']])

start = time.time()

for t in range(num_iterations):
	for i in range(num_nests):
		# Generate first candidate replacement and replace if better
		u_1 = levy_flight(nests['nest'][i])
		f_u_1 = f(u_1)
		if f_u_1 <= nests['fitness'][i]:
			nests['nest'][i] = u_1
			nests['fitness'][i] = f_u_1
		# With random chance pa produce a second candidate replacement and replace if better
		if np.random.uniform(0, 1) <= pa:
			u_2 = levy_flight(nests['nest'][i])
			f_u_2 = f(u_2)
			if not parasitism_comparison or f_u_2 <= nests['fitness'][i]:  # parasitism_comparison toggles the parasitism comparison
				nests['nest'][i] = u_2
				nests['fitness'][i] = f_u_2

best = nests[nests['fitness'].argmin()]
print(best)
print("Time taken:", time.time() - start)