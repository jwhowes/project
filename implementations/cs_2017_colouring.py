import numpy as np
from math import sin, gamma, pi

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
num_iterations = 1000
parasitism_comparison = True

nests = np.random.randint(0, k, (num_nests, n))

for t in range(num_iterations):  # This is pretty trash, I feel like it's not aggressive enough
	for i in range(num_nests):
		u_1 = levy_flight(nests[i])
		if f(u_1) <= f(nests[i]):
			nests[i] = u_1
		if np.random.uniform(0, 1) <= pa:
			u_2 = levy_flight(nests[i])
			if not parasitism_comparison or f(u_2) <= f(nests[i]):
				nests[i] = u_2

best = nests[np.array([f(n) for n in nests]).argmin()]
print(best, f(best))