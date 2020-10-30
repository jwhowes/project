import numpy as np

k = 3
n_pop = 30
n = 7

cuckoos = np.random.randint(0, k, (n_pop, n))

def d(x, y):
	ret = 0
	for i in range(n):
		if x[i] != y[i]:
			ret += 1
	return ret

def d_bar(S1, S2):
	num = 0
	for i in S1:
		for j in S2:
			num += d(i, j)
	return num/(len(S1)*len(S2))

def tri_dist(x, i):
	"""Finds the triangular distance between x and cluster i"""
	Si = cuckoos[np.where(S == i)]
	if len(Si) == 0:
		return -1
	return 2*d_bar([x], Si) - d_bar(Si, Si)

# Assign each cuckoo to a random cluster
S = np.random.randint(0, k, (n_pop))
tri_dist(cuckoos[0], 0)

print(S)

converged = False
while not converged:
	converged = True
	for i in range(n_pop):
		cluster = -1
		tdc = -1
		for j in range(k):
			td = tri_dist(cuckoos[i], j)
			if (cluster == -1 or td < tdc) and td != -1:
				cluster = j
				tdc = td
		if cluster != S[i]:
			converged = False
			S[i] = cluster

print(S)