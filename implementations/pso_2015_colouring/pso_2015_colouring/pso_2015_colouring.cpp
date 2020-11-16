#include <iostream>
#include <vector>
#include <array>
#include <chrono>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>

// Running time (population 50, 3000 iterations):
	// 500 vertices: 246.89 secs
		// Couldn't colour with 200 colours (my_coa_basic coloured 1000 with 143 colours)

using namespace std;
using namespace boost::random;

const int k = 200;
const int NUM_VERTICES = 500;
int m = 0;

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

const int num_particles = 50;
const int num_iterations = 3000;

const float w = 0.05;
const float c1 = 7;
const float c2 = 0.03;

int particles[num_particles][NUM_VERTICES];  // The population of particles
int particles_old[num_particles][NUM_VERTICES];  // The previous generation of particles
int p[num_particles][NUM_VERTICES];  // The personal best for each particle
int g[NUM_VERTICES];  // The global best particle found

mt19937 seed;
uniform_real_distribution<float> uni(0, 1);
uniform_int_distribution<int> random_colour(0, k - 1);
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

float f(int * x) {
	int num = 0;
	for (int i = 0; i < NUM_VERTICES; i++) {
		for (int j = 0; j < i; j++) {
			if (adj_matrix[i][j] == 1 && x[i] == x[j]) {
				num++;
			}
		}
	}
	return 1 - (num / (float)m);
}

float d(int * x, int * y) {
	int num = 0;
	for (int i = 0; i < NUM_VERTICES; i++) {
		if (x[i] != y[i]) {
			num++;
		}
	}
	return 1 - num / (float)NUM_VERTICES;
}

void generate_particle(int * particle) {
	for (int i = 0; i < NUM_VERTICES; i++) {
		particle[i] = random_colour(seed);
	}
}

int main(){  // Could probably precompute some fitness values (not as many as in coa but still could)
	make_graph(0.5);
	// Calculate number of edges
	for (int i = 0; i < NUM_VERTICES; i++) {
		for (int j = 0; j < i; j++) {
			m += adj_matrix[i][j];
		}
	}
	auto start = chrono::high_resolution_clock::now();
	// Generate initial population
	float g_fitness = 0;
	for (int i = 0; i < num_particles; i++) {
		generate_particle(particles[i]);
		copy(begin(particles[i]), end(particles[i]), begin(p[i]));
		if (f(particles[i]) > g_fitness) {
			copy(begin(particles[i]), end(particles[i]), begin(g));
			g_fitness = f(particles[i]);
		}
	}
	// Begin main loop
	for (int t = 0; t < num_iterations; t++) {
		for (int i = 0; i < num_particles; i++) {
			float r1 = uni(seed);
			float r2 = uni(seed);
			float v_rand = w * d(particles[i], particles_old[i]);
			float v_p = c1 * r1 * d(particles[i], p[i]);
			float v_g = c2 * r2 * d(particles[i], g);
			float V = v_rand + v_p + v_g;
			float p_rand = v_rand / V; float p_p = v_p / V;  // No need to actually calculate p_g
			copy(begin(particles[i]), end(particles[i]), begin(particles_old[i]));
			for (int j = 0; j < NUM_VERTICES; j++) {
				float r = uni(seed);
				if (r <= p_rand) {
					particles[i][j] = random_colour(seed);
				} else if (r <= p_rand + p_p) {
					particles[i][j] = p[i][j];
				} else {
					particles[i][j] = g[j];
				}
			}
			if (f(particles[i]) > f(p[i])) {
				copy(begin(particles[i]), end(particles[i]), begin(p[i]));
				if (f(particles[i]) > f(g)) {
					copy(begin(particles[i]), end(particles[i]), begin(g));
				}
			}
		}
		//cout << chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start).count() << endl;
	}
	for (int i = 0; i < NUM_VERTICES; i++) {
		cout << g[i] << " ";
	}
	cout << endl << "Number of conflicts: " << (1 - f(g)) * m << endl;
	cout << "Time taken (seconds): " << chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start).count() / (float)1000000 << endl;
	return 0;
}
