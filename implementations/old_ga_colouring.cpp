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

const int k = 3;
const int NUM_VERTICES = 10;

struct Partition {
	vector<int> partition[k];
	int fitness;
};

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

const int pop_size = 50;
const int num_iterations = 3000;
const float mutation_prob = 0.1f;

Partition population[pop_size];
Partition new_population[pop_size];

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

int f(vector<int> * x) {
	int num = 0;
	for (int c = 0; c < k; c++) {
		for (int i = 0; i < x[c].size(); i++) {
			for (int j = 0; j < i; j++) {
				if (adj_matrix[x[c][i]][x[c][j]] == 1) {
					num++;
				}
			}
		}
	}
	return num;
}

vector<int> parents[2][k];
void gpx(vector<int> * p1, vector<int> * p2, vector<int> * x) {
	vector<int> uncoloured(NUM_VERTICES);
	int parent;
	for (int i = 0; i < NUM_VERTICES; i++) {
		uncoloured[i] = i;
	}
	for (int i = 0; i < k; i++) {
		parents[0][i].assign(p1[i].begin(), p1[i].end());
		parents[1][i].assign(p2[i].begin(), p2[i].end());
		x[i].clear();
	}
	for (int i = 0; i < k; i++) {
		parent = i % 2;
		int m = 0;
		for (int j = 1; j < k; j++) {
			if (parents[parent][j].size() > parents[parent][m].size()) {
				m = j;
			}
		}
		x[i].assign(parents[parent][m].begin(), parents[parent][m].end());
		for (int i = 0; i < parents[parent][m].size(); i++) {
			uncoloured.erase(remove(uncoloured.begin(), uncoloured.end(), parents[parent][m][i]), uncoloured.end());
		}
		parents[parent][m].clear();
		for (int j = 0; j < x[i].size(); j++) {
			for (int l = 0; l < k; l++) {
				auto pos = find(parents[1 - parent][l].begin(), parents[1 - parent][l].end(), x[i][j]);
				if (pos != parents[1 - parent][l].end()) {
					parents[1 - parent][l].erase(pos);
					break;
				}
			}
		}
	}
	for (int i = 0; i < uncoloured.size(); i++) {
		x[random_colour(seed)].push_back(uncoloured[i]);
	}
}

void generate_member(vector<int> * x) {
	for (int i = 0; i < NUM_VERTICES; i++) {
		x[random_colour(seed)].push_back(i);
	}
}

void mutate(vector<int> * x) {  // Randomly moves a vertex from one colour class to another
	int c = random_colour(seed);
	while (x[c].size() == 0) {
		c = (c + 1) % k;
	}
	int d = random_colour(seed);
	int v = x[c][uniform_int_distribution<int>(0, x[c].size() - 1)(seed)];
	x[c].erase(remove(x[c].begin(), x[c].end(), v), x[c].end());
	x[d].push_back(v);
}

bool compare_partitions(Partition & x, Partition & y) {
	return x.fitness < y.fitness;
}

void get_parents(int * p1, int * p2, int F) {
	float r = uni(seed);
	for (int i = 0; i < pop_size; i++) {
		r -= (population[pop_size - 1].fitness - population[i].fitness) / (float)F;
		if (r <= 0.0f) {
			*p1 = i;
			break;
		}
	}
	r = uni(seed);
	for (int i = 0; i < pop_size; i++) {
		r -= (population[pop_size - 1].fitness - population[i].fitness) / (float)F;
		if (r <= 0.0f) {
			*p2 = i;
			break;
		}
	}
}

int main() {
	make_graph(0.5);
	int p1; int p2;
	int F;
	for (int i = 0; i < pop_size; i++) {
		generate_member(population[i].partition);
		population[i].fitness = f(population[i].partition);
	}
	auto start = chrono::high_resolution_clock::now();
	for (int t = 0; t < num_iterations; t++) {
		sort(begin(population), end(population), compare_partitions);  // You need to pass values through by reference (&), see if there's some way of doing this without defining a Partition struct
		F = 0;
		for (int i = 0; i < pop_size; i++) {
			F += (population[pop_size - 1].fitness - population[i].fitness);
		}
		for (int i = 0; i < pop_size; i++) {
			get_parents(&p1, &p2, F);
			gpx(population[p1].partition, population[p2].partition, new_population[i].partition);
			if (uni(seed) < mutation_prob) {
				mutate(new_population[i].partition);
			}
		}
		for (int i = 0; i < pop_size; i++) {
			copy(begin(new_population[i].partition), end(new_population[i].partition), begin(population[i].partition));
			population[i].fitness = f(population[i].partition);
		}
	}
	for (int i = 0; i < k; i++) {
		for (int j = 0; j < population[0].partition[i].size(); j++) {
			cout << population[0].partition[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl << "Number of conflicts: " << population[0].fitness;
	cout << endl << "Time taken (seconds): " << chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start).count() / (float)1000000 << endl;
	return 0;
}