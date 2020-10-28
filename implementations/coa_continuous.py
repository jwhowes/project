# TODO:
#	- Add clipping
import numpy as np
from math import sin, pi, sqrt

dimension = 2
v_lo = 0
v_hi = 5
def f(x):  # Min is at (2.20319, 1.57049)
	return -sin(x[0])*(sin(x[0]*x[0]/pi)**20) - sin(x[1])*(sin(2*x[1]*x[1]/pi)**20)

k = 5
def goal_point():  # I should try a different initialisation algorithm (convergence is too quick with random init)
	"""Clusters the cuckoos and returns the mean of the best cluster"""
	# I should also make the goal point the best point of the best cluster (rather than the mean of the best cluster)
	# Initialise means
	#	Using uniform random initialisation
	means = np.random.uniform(v_lo, v_hi, (k, dimension))
	converged = False
	# Initialise cluster array (will be populated later)
	clusters = np.zeros((n_pop), dtype=int)  # clusters[i] is the cluster cuckoo i belongs to
	while not converged:
		converged = True
		# Assign each cuckoo to the mean closest to it
		for i in range(n_pop):
			cluster = 0
			for j in range(1, k):
				if np.linalg.norm(cuckoos[i] - means[j]) < np.linalg.norm(cuckoos[i] - means[cluster]):
					cluster = j
			if cluster != clusters[i]:
				converged = False
				clusters[i] = cluster
		# Reset each mean value to the mean of its cluster
		_means = np.zeros((k, dimension))  # A placeholder for the new means (if a cluster is empty we don't want to overwrite its mean)
		num = np.zeros((k))  # num[i] is the number of cuckoos in cluster i
		for i in range(n_pop):
			_means[clusters[i]] += cuckoos[i]
			num[clusters[i]] += 1
		for i in range(k):
			if num[i] != 0:
				means[i] = _means[i]/num[i]
	best = 0
	for i in range(1, k):
		if f(means[i]) < f(means[best]):
			best = i
	return means[best]



alpha = 0.2
n_pop = 30
n_max = 30
num_iterations = 500
p = 0.1
F = 1

cuckoos = np.random.uniform(v_lo, v_hi, (n_pop, dimension))

for t in range(num_iterations):  # Performs surprisingly well considering I haven't implemented k-means
	# Assign each cuckoo a random number of eggs
	num_eggs = np.random.randint(5, 21, (n_pop))
	# Lay eggs
	tot_eggs = num_eggs.sum()
	eggs = np.zeros((tot_eggs, dimension))
	egg = 0
	for i in range(n_pop):
		ELR = alpha * num_eggs[i]/tot_eggs * (v_hi - v_lo)
		for j in range(num_eggs[i]):
			d = np.random.uniform(0, ELR)  # The distance the egg will be from the cuckoo
			v = np.random.uniform(0, 1, (dimension))
			v /= np.linalg.norm(v)  # v is a random unit vector
			eggs[egg] = cuckoos[i] + v * d  # Generates a random vector distance d from cuckoos[i] and appends it to eggs
			egg += 1
	# Kill p of the worst eggs
	#	Or just add the best n(1 - p) eggs to cuckoos
	eggs = eggs[np.array([f(e) for e in eggs]).argsort()]
	cuckoos = np.append(eggs[:int(tot_eggs * (1 - p))], cuckoos, axis=0)
	# Kill worst cuckoos until pop <= n_max (probably kill the eggs before adding them to pop)
	cuckoos = cuckoos[np.array([f(c) for c in cuckoos]).argsort()]
	if len(cuckoos) > n_max:
		cuckoos.resize((n_max, dimension))
	n_pop = len(cuckoos)
	# Cluster cuckoos and find goal point
	gp = goal_point()
	# Migrate cuckoos
	#	Using the linear random method (as it generalises to > 2 dimensions)
	for i in range(n_pop):
		r = np.random.uniform(0, 1)
		cuckoos[i] = cuckoos[i] + F * r * (gp - cuckoos[i])

cuckoos = cuckoos[np.array([f(c) for c in cuckoos]).argsort()]
print(cuckoos[0], f(cuckoos[0]))