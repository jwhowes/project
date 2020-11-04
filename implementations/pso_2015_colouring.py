import numpy as np

adj_matrix = [  [0, 1, 1, 0, 1, 0, 0],
				[1, 0, 0, 0, 1, 0, 1],
				[1, 0, 0, 0, 0, 1, 1],
				[0, 0, 0, 0, 0, 0, 1],
				[1, 1, 0, 0, 0, 0, 0],
				[0, 0, 1, 0, 0, 0, 0],
				[0, 1, 1, 1, 0, 0, 0]]
k = 3
n = len(adj_matrix)
m = np.sum(adj_matrix) // 2

def f(x):
	num = 0
	for i in range(len(adj_matrix)):
		for j in range(i):
			if adj_matrix[i][j] == 1 and x[i] == x[j]:
				num += 1
	return 1 - num/m

def d(x, y):
	num = 0
	for i in range(n):
		if x[i] != y[i]:
			num += 1
	return 1 - num/n

num_particles = 50
num_iterations = 1000

w = 0.05
c1 = 7
c2 = 0.03

particles = np.random.randint(0, k, (num_particles, n))
p = particles.copy()
g = particles[np.array([f(p) for p in particles]).argmax()]
particles_old = np.zeros((num_particles, n))

for t in range(num_iterations):
	for i in range(num_particles):
		r1 = np.random.uniform(0, 1)
		r2 = np.random.uniform(0, 1)
		v_rand = w * d(particles[i], particles_old[i])
		v_p = c1 * r1 * d(particles[i], p[i])
		v_g = c2 * r2 * d(particles[i], g)
		V = v_rand + v_p + v_g
		p_rand = v_rand/V
		p_p = v_p/V
		p_g = v_g/V
		particles_old[i] = particles[i]
		for j in range(n):
			r = np.random.uniform(0, 1)
			if r <= p_rand:
				particles[i][j] = np.random.randint(0, k)
			elif r <= p_rand + p_p:
				particles[i][j] = p[i][j]
			else:
				particles[i][j] = g[j]
		if f(particles[i]) > f(p[i]):
			p[i] = particles[i]
			if f(particles[i]) > f(g):
				g = particles[i]

print(g, (1 - f(g))*m)  # Print best solution and number of conflicts