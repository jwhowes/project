#define _SECURE_SCL 0

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

const int NUM_VERTICES = 350;

int adj_matrix[NUM_VERTICES][NUM_VERTICES];/* = {
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
};*/

float G = 1.0f;

const int num_particles = 50;
const int num_iterations = 3000;

int dist_matrix[NUM_VERTICES][NUM_VERTICES];

struct Particle {
	int colouring[NUM_VERTICES];
	float mass;
};

Particle particles[num_particles];

mt19937 seed;
uniform_real_distribution<float> uni(0, 1);
uniform_int_distribution<int> random_vertex(0, NUM_VERTICES - 1);

void make_graph(float edge_probability) {  // Populates adj_matrix with a random graph
	for (int i = 0; i < NUM_VERTICES; i++) {
		for (int j = 0; j < i; j++) {
			if (uni(seed) < edge_probability) {
				adj_matrix[i][j] = 1;
				adj_matrix[j][i] = 1;
			}
		}
	}
}

int f(int * x) {
	int max = 0;
	for (int i = 0; i < NUM_VERTICES; i++) {
		if (x[i] > max) {
			max = x[i];
		}
	}
	return max + 1;
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
		if (x[i] != y[i]) {
			num++;
		}
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

int I[NUM_VERTICES];
void migrate(int * x, int * y, float dist) {  // Migrates x towards y
	// Populate I with all vertices on which x and y disagree
	int I_length = 0;
	for (int i = 0; i < NUM_VERTICES; i++) {
		if (x[i] != y[i]) {
			I[I_length] = i;
			I_length++;
		}
	}
	for (int i = 0; i < dist * I_length; i++) {  // For a random
		int v = I[i];
		// Assign x y's colour for *it
		x[v] = y[v];
		// Clean up x
		for (int j = i + 1; j < I_length; j++) {
			// For every vertex I[j] in I with I-index j > i that conflicts with v, assign the smallest legal colour different thatn y[I[j]]
			if (x[I[j]] == x[v] && adj_matrix[I[j]][v] == 1) {//find(edge_list[v].begin(), edge_list[v].end(), j) != edge_list[v].end()) {
				int c = 0;
				while (true) {
					if (c != y[I[j]] && valid(I[j], c, x)) {
						x[I[j]] = c;
						break;
					}
					c++;
				}
			}
		}
	}
}

bool compare_particles(Particle & p1, Particle & p2) {
	return f(p1.colouring) < f(p2.colouring);
}

float F[NUM_VERTICES];
int main(){
	make_graph(0.5);
	for (int i = 0; i < NUM_VERTICES; i++) {
		order[i] = i;
	}
	for (int i = 0; i < num_particles; i++) {
		generate_particle(particles[i].colouring);
	}
	for (int t = 0; t < num_iterations; t++) {
		update_mass();
		for (int p = 0; p < num_particles; p++) {
			float F_p = 0;
			for (int q = 0; q < num_particles; q++) {
				if (p != q) {
					F[q] = G * ((particles[q].mass) / (float)dist_matrix[p][q]);  // F[q] stores F_pq
					F_p += F[q];
				}
			}
			for (int q = 0; q < num_particles; q++) {
				if (p != q) {
					migrate(particles[p].colouring, particles[q].colouring, F[q]/F_p);
				}
			}
		}
	}
	sort(begin(particles), end(particles), compare_particles);
	for (int i = 0; i < NUM_VERTICES; i++) {
		cout << particles[0].colouring[i] << " ";
	}
	cout << endl << "Number of colours: " << f(particles[0].colouring);
	return 0;
}