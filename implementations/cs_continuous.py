import numpy as np
from math import sin, gamma, pi

dimension = 2
v_lo = 0
v_hi = 5
def f(x):  # Max is at (2.20319, 1.57049)
	return -sin(x[0])*(sin(x[0]*x[0]/pi)**20) - sin(x[1])*(sin(2*x[1]*x[1]/pi)**20)


pa = 0.25
num_nests = 50
num_iterations = 1000

beta = 1.5
alpha = 0.2

sigma_p = 	((gamma(1+beta)*sin(pi*beta/2))/
			(gamma((1+beta)/2)*beta*2**((beta-1)/2)))**(2/beta)  # The variance of p

def levy():
	"""Generates a number with a levy distribution"""
	p = np.random.normal(0, sigma_p)
	q = np.random.normal(0, 1)
	return p/(abs(q)**(1/beta))

def levy_flight(nest):
	"""Returns the result of performing a levy flight on nest"""
	ret = nest.copy()
	for i in range(dimension):
		ret += alpha*levy()
	return ret

nests = np.random.uniform(v_lo, v_hi, (num_nests, dimension))  # Generate initial population uniformly at random

nests = nests[np.array([f(n) for n in nests]).argsort()]

for t in range(num_iterations):
	for i in range(num_nests):
		j = np.clip(levy_flight(nests[i]), v_lo, v_hi)
		if f(j) < f(nests[i]):
			nests[i] = j
	nests = nests[np.array([f(n) for n in nests]).argsort()]
	for i in range(int(num_nests * pa)):
		nests[num_nests - i - 1] = np.random.uniform(v_lo, v_hi, (dimension))

nest = nests[np.array([f(n) for n in nests]).argmin()]
print(nest, f(nest))