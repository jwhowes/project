import numpy as np

adj_matrix = [  [0, 1, 1, 1, 1, 1],
				[1, 0, 1, 1, 1, 1],
				[1, 1, 0, 1, 1, 1],
				[1, 1, 1, 0, 1, 1],
				[1, 1, 1, 1, 0, 1],
				[1, 1, 1, 1, 1, 0]]
k = 6
n = len(adj_matrix)

def f(x):
	"""Returns the number of conflicts in a given colouring x"""
	num = 0
	for i in range(len(adj_matrix)):
		for j in range(i):
			if adj_matrix[i][j] == 1 and x[i] == x[j]:
				num += 1
	return num

class TabuList:
	def __init__(self, tenure):
		self.front = 0
		self.ls = [(0, 0)] * tenure
		self.tenure = tenure
		for i in range(tenure):
			self.ls[i] = (np.random.randint(0, n), np.random.randint(0, k))
	def contains(self, x):
		return x in self.ls
	def add(self, x):
		self.ls[self.front] = x
		self.front = (self.front + 1) % self.tenure

def A(z):
	if z in A_dict:
		return A_dict[z]
	return z - 1

def get_neighbour():
	best = None
	best_v = 0
	best_c = 0
	for v in range(n):
		for c in range(k):
			neighbour = s.copy()
			neighbour[v] = c
			if f(neighbour) <= A(f(s)):
				A_dict[f(s)] = f(neighbour) - 1
				if f(neighbour) < f(s):
					T.add((v, c))
					return neighbour
				best = neighbour
				best_v = v
				best_c = c
			elif not T.contains((v, c)):
				if f(neighbour) < f(s):
					T.add((v, c))
					return neighbour
				if best is None or f(neighbour) < f(best):
					best = neighbour
					best_v = v
					best_c = c
	T.add((best_v, best_c))
	return best

# Generate initial random tabu list (|T| = 7)
T = TabuList(7)

# Generate a random solution
s = np.random.randint(0, k, (n))

A_dict = {}

t = 0
num_iterations = 100

while f(s) > 0 and t < num_iterations:
	s = get_neighbour()  # Could probably do this inline (I thought there would be more after this)
	t += 1

print(s, f(s))