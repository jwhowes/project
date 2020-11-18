#define _SECURE_SCL 0

// Running time (k = 0.2V):
	// 100 vertices: 76.3081 secs

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

const int k = 20;
const int NUM_VERTICES = 100;

struct Partition {
	int partition[k][NUM_VERTICES];
	int partition_length[k];
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

int f(Partition & x) {
	int num = 0;
	for (int c = 0; c < k; c++) {
		for (int i = 0; i < x.partition_length[c]; i++) {
			for (int j = 0; j < i; j++) {
				if (adj_matrix[x.partition[c][i]][x.partition[c][j]] == 1) {
					num++;
				}
			}
		}
	}
	return num;
}

Partition parents[2];
void gpx(Partition & p1, Partition & p2, Partition & x) {
	vector<int> uncoloured(NUM_VERTICES);
	int parent;
	for (int i = 0; i < NUM_VERTICES; i++) {
		uncoloured[i] = i;
	}
	for (int i = 0; i < k; i++) {
		copy(begin(p1.partition[i]), begin(p1.partition[i]) + p1.partition_length[i], begin(parents[0].partition[i]));
		parents[0].partition_length[i] = p1.partition_length[i];
		copy(begin(p2.partition[i]), begin(p2.partition[i]) + p2.partition_length[i], begin(parents[1].partition[i]));
		parents[1].partition_length[i] = p2.partition_length[i];
		x.partition_length[i] = 0;
	}
	for (int i = 0; i < k; i++) {
		parent = i % 2;
		int m = 0;
		for (int j = 1; j < k; j++) {
			if (parents[parent].partition_length[j] > parents[parent].partition_length[m]) {
				m = j;
			}
		}
		copy(begin(parents[parent].partition[m]), begin(parents[parent].partition[m]) + parents[parent].partition_length[m], begin(x.partition[i]));
		x.partition_length[i] = parents[parent].partition_length[m];
		for (int i = 0; i < parents[parent].partition_length[m]; i++) {
			uncoloured.erase(remove(uncoloured.begin(), uncoloured.end(), parents[parent].partition[m][i]), uncoloured.end());
		}
		parents[parent].partition_length[m] = 0;
		// Remove everything from x.partition[i] from parents[1-parent].partition
		for (int j = 0; j < x.partition_length[i]; j++) {
			for (int l = 0; l < k; l++) {
				auto pos = find(begin(parents[1 - parent].partition[l]), begin(parents[1 - parent].partition[l]) + parents[1 - parent].partition_length[l], x.partition[i][j]);
				//auto pos = find(parents[1-parent][l].begin(), parents[1-parent][l].end(), x[i][j]);
				if (pos != begin(parents[1 - parent].partition[l]) + parents[1 - parent].partition_length[l]) {
					copy(pos + 1, begin(parents[1 - parent].partition[l]) + parents[1 - parent].partition_length[l], pos);
					parents[1 - parent].partition_length[l]--;
					//parents[1-parent][l].erase(pos);
					break;
				}
			}
		}
	}
	for (int i = 0; i < uncoloured.size(); i++) {
		int c = random_colour(seed);
		x.partition[c][x.partition_length[c]] = uncoloured[i];
		x.partition_length[c]++;
	}
}

void generate_member(Partition & x) {
	for (int i = 0; i < NUM_VERTICES; i++) {
		int c = random_colour(seed);
		x.partition[c][x.partition_length[c]] = i;
		x.partition_length[c]++;
		//x[random_colour(seed)].push_back(i);
	}
}

void mutate(Partition & x) {  // Randomly moves a vertex from one colour class to another
	int c = random_colour(seed);
	while (x.partition_length[c] == 0) {
		c = (c + 1) % k;
	}
	int d = random_colour(seed);
	int v = x.partition[c][uniform_int_distribution<int>(0, x.partition_length[c] - 1)(seed)];
	auto pos = find(begin(x.partition[c]), begin(x.partition[c]) + x.partition_length[c], v);
	copy(pos + 1, begin(x.partition[c]) + x.partition_length[c], pos);
	x.partition_length[c]--;
	//x[c].erase(remove(x[c].begin(), x[c].end(), v), x[c].end());
	x.partition[d][x.partition_length[d]] = v;
	x.partition_length[d]++;
	//x[d].push_back(v);
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
		generate_member(population[i]);
		population[i].fitness = f(population[i]);
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
			gpx(population[p1], population[p2], new_population[i]);
			if (uni(seed) < mutation_prob) {
				mutate(new_population[i]);
			}
		}
		for (int i = 0; i < pop_size; i++) {
			for (int j = 0; j < k; j++) {
				copy(begin(new_population[i].partition[j]), begin(new_population[i].partition[j]) + new_population[i].partition_length[j], begin(population[i].partition[j]));
			}
			copy(begin(new_population[i].partition_length), end(new_population[i].partition_length), begin(population[i].partition_length));
			population[i].fitness = f(population[i]);
		}
		//cout << chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start).count() << endl;
	}
	sort(begin(population), end(population), compare_partitions);
	for (int i = 0; i < k; i++) {
		for (int j = 0; j < population[0].partition_length[i]; j++) {
			cout << population[0].partition[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl << "Number of conflicts: " << population[0].fitness;
	cout << endl << "Time taken (seconds): " << chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start).count() / (float)1000000 << endl;
	return 0;
}