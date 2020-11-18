#define _SECURE_SCL 0

// TODO:
	// See if you can make the inner loop run for nC2 rather than n^2 (similar to how you did it with dist_matrix)
	// For some reason uni(seed) isn't terminating?
		// For now I've just done it with 0.5 replacing it each time
	// Consider some kind of local search at the end (to fix any conflicts)
		// Remember local optimum of this fitness function is a valid colouring
	// Potentially consider adding a minimise function (to make fitness calculating more efficient as it takes it as a partition)

#include <iostream>
#include <array>
#include <vector>
#include <math.h>
#include <chrono>
#include <algorithm>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>

using namespace std;
using namespace boost::random;

const int NUM_VERTICES = 10;

int adj_matrix[NUM_VERTICES][NUM_VERTICES] = {
	{0, 1, 0, 0, 1, 1, 0, 0, 0, 0},
	{1, 0, 1, 0, 0, 0, 1, 0, 0, 0},
	{0, 1, 0, 1, 0, 0, 0, 1, 0, 0},
	{0, 0, 1, 0, 1, 0, 0, 0, 1, 0},
	{1, 0, 0, 1, 0, 0, 0, 0, 0, 1},
	{1, 0, 0, 0, 0, 0, 0, 1, 1, 0},
	{0, 1, 0, 0, 0, 0, 0, 0, 1, 1},
	{0, 0, 1, 0, 0, 1, 0, 0, 0, 1},
	{0, 0, 0, 1, 0, 1, 1, 0, 0, 0},
	{0, 0, 0, 0, 1, 0, 1, 1, 0, 0}
};

float G = 1.0f;

const int num_particles = 50;
const int num_iterations = 3000;

int dist_matrix[NUM_VERTICES][NUM_VERTICES];

struct Particle {
	int colouring[NUM_VERTICES];
	int velocity[NUM_VERTICES];
	float mass;
};

Particle particles[num_particles];

mt19937 seed;
uniform_real_distribution<float> uni(0, 1);
uniform_int_distribution<int> random_vertex(0, NUM_VERTICES - 1);

int num_colours(int * x) {
	int num = 0;
	for (int i = 0; i < NUM_VERTICES; i++) {
		if (x[i] > num) {
			num = x[i];
		}
	}
	return num + 1;
}

int * colour_class_size;
int * edge_conflicts;
int f(int * x) {
	delete[] colour_class_size;
	delete[] edge_conflicts;
	int n = num_colours(x);
	colour_class_size = new int[n];
	edge_conflicts = new int[n];
	for (int i = 0; i < n; i++) {
		colour_class_size[i] = 0;
		edge_conflicts[i] = 0;
	}
	for (int i = 0; i < NUM_VERTICES; i++) {
		colour_class_size[x[i]]++;
		for (int j = 0; j < i; j++) {
			if (adj_matrix[i][j] == 1 && x[i] == x[j]) {
				edge_conflicts[x[i]]++;
			}
		}
	}
	int ret = 0;
	for (int i = 0; i < n; i++) {
		ret += 2 * colour_class_size[i] * edge_conflicts[i] - colour_class_size[i] * colour_class_size[i];
	}
	return ret;
}

bool valid(int v, int c, int * col) {  // Returns whether or not vertex v can be coloured colour c in colouring col (legally)
	for (int i = 0; i < NUM_VERTICES; i++) {
		if (adj_matrix[v][i] == 1 && col[i] == c) {  // If v is adjacent to some vertex i coloured c then this is not a valid assignment
			return false;
		}
	}
	return true;  // If no such i can be found then the assignment is valid
}

int order[NUM_VERTICES];
void generate_particle(int * particle) {
	random_shuffle(begin(order), end(order));
	for (int v : order) {
		int c = 0;
		while (true) {
			if (valid(v, c, particle)) {
				particle[v] = c;
				break;
			}
			c++;
		}
	}
}

int fitness[NUM_VERTICES];
void update_mass() {
	int worst = -1; int best = -1;
	for (int i = 0; i < num_particles; i++) {
		fitness[i] = f(particles[i].colouring);
		if (worst == -1 || fitness[i] > worst) {
			worst = fitness[i];
		}
		if (best == -1 || fitness[i] < best) {
			best = fitness[i];
		}
	}
	float m_sum = 0;
	for (int i = 0; i < num_particles; i++) {
		m_sum += (fitness[i] - worst) / (float)(best - worst);
	}
	for (int i = 0; i < num_particles; i++) {
		particles[i].mass = ((fitness[i] - worst) / (float)(best - worst)) / m_sum;
	}
}

int d(int * x, int * y) {
	int num = 0;
	for (int i = 0; i < NUM_VERTICES; i++) {
		num += abs(x[i] - y[i]);
	}
	return num;
}

void populate_dist_matrix() {
	for (int i = 0; i < NUM_VERTICES; i++) {
		for (int j = 0; j < i; j++) {
			dist_matrix[i][j] = d(particles[i].colouring, particles[j].colouring);
			dist_matrix[j][i] = dist_matrix[i][j];
		}
	}
}

void force(int p, int q, float * F) {  // Probably precompute distance matrix
	for (int i = 0; i < NUM_VERTICES; i++) {
		if (dist_matrix[p][q] != 0) {
			F[i] = G * ((particles[p].mass*particles[q].mass) / dist_matrix[p][q]) * (particles[q].colouring[i] - particles[p].colouring[i]);
		} else {
			F[i] = 0;
		}
	}
}

int main(){
	float F_pq[NUM_VERTICES];
	float F_p[NUM_VERTICES];
	for (int i = 0; i < NUM_VERTICES; i++) {
		order[i] = i;
	}
	for (int i = 0; i < num_particles; i++) {
		generate_particle(particles[i].colouring);
	}
	update_mass();
	for (int t = 0; t < num_iterations; t++) {
		populate_dist_matrix();
		for (int p = 0; p < num_particles; p++) {
			for (int i = 0; i < NUM_VERTICES; i++) {
				F_p[i] = 0;
			}
			for (int q = 0; q < num_particles; q++) {  // See of you can find some way to only do nC2 of these rather than n^2
				if (p != q) {
					// Find F_{pq} and put it into F
					force(p, q, F_pq);
					// Update F_p
					for (int i = 0; i < NUM_VERTICES; i++) {
						F_p[i] += 0.5 * F_pq[i];
					}
				}
			}
			// Update velocity and colouring
			for (int i = 0; i < NUM_VERTICES; i++) {
				particles[p].velocity[i] = 0.5 * particles[p].velocity[i] + (int)(F_p[i] / particles[i].mass);
			}
		}
		// Update colourings
		for (int p = 0; p < num_particles; p++) {
			for (int i = 0; i < NUM_VERTICES; i++) {
				particles[p].colouring[i] += particles[p].velocity[i];
				if (particles[p].colouring[i] < 0) {
					particles[p].colouring[i] = 0;
				}
			}
		}
		update_mass();
	}
	for (int i = 0; i < NUM_VERTICES; i++) {
		cout << particles[0].colouring[i] << " ";
	}
	cout << endl;
	return 0;
}