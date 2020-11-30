#define _USE_MATH_DEFINES
#define _SECURE_SCL 0

// The bottleneck is the eta function

// Starting from highest degree vertex:
	// 100 vertices:
		// 20 colours: 0.8017 secs (0 conflicts)
		// ...
		// 17 colours: 2.03416 secs (2 conflicts)
	// 500 vertices:
		// 75 colours: 35.3846 secs (9 conflicts)
		// 70 colours: 33.517 secs (22 conflicts)

// Starting from random vertex:
	// 100 vertices:
		// 20 colours: 0.50653 secs (0 conflicts)
		// 19 colours: 1.37949 secs (0 conflicts)
		// 18 colours: 2.00710 secs (1 conflict)
	// 500 vertices:
		// 75 colours: 33.6179 secs (11 conflicts)
		// 70 colours: 32.5405 secs (26 conflicts)

// Starting from highest degree, setting M = num_vertices:
	// 100 vertices:
		// 20 colours: 0.24180 secs (0 conflicts)
		// 19 colours: 0.43910 secs (0 conflicts)
		// 18 colours: 0.63418 secs (0 conflicts)
		// 17 colours: It crashed?

// Starting from highest degree, alpha = 1, beta = 0.5
	// 100 vertices:
		// 20 colours: 0.14728 secs (0 conflicts)
		// 19 colours: 0.11779 secs (0 conflicts)
		// 18 colours: 0.40009 secs (0 conflicts)
		// 17 colours: 0.93647 secs (0 conflicts)
		// 16 colours: 87.5969 secs (2 conflicts)
	// 500 vertices:
		// 75 colours: 13.8668 secs (0 conflicts)
		// 70 colours: 28.6291 secs (0 conflicts)
		// 69 colours: 37.6985 secs (0 conflicts)
		// 68 colours: 34.3267 secs (0 conflicts)
		// 67 colours: 53.7195 secs (0 conflicts)
		// 66 colours: 74.6965 secs (0 conflicts), 75 iterations (estimated 49.7977 mins for 3000 iterations)

#include <iostream>
#include <array>
#include <vector>
#include <math.h>
#include <chrono>
#include <algorithm>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>

using namespace std;
using namespace boost::random;

const int num_vertices = 100;
int adj_matrix[num_vertices][num_vertices];/* = {
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
int adj_list[num_vertices][num_vertices];/* = {
	{1, 4, 5, 0, 0, 0, 0, 0, 0, 0},
	{0, 2, 6, 0, 0, 0, 0, 0, 0, 0},
	{1, 3, 7, 0, 0, 0, 0, 0, 0, 0},
	{2, 4, 8, 0, 0, 0, 0, 0, 0, 0},
	{0, 3, 9, 0, 0, 0, 0, 0, 0, 0},
	{0, 7, 8, 0, 0, 0, 0, 0, 0, 0},
	{1, 8, 9, 0, 0, 0, 0, 0, 0, 0},
	{2, 5, 9, 0, 0, 0, 0, 0, 0, 0},
	{3, 5, 6, 0, 0, 0, 0, 0, 0, 0},
	{4, 6, 7, 0, 0, 0, 0, 0, 0, 0}
};*/
int adj_list_length[num_vertices];// = { 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 };
const int k = 18;

const int num_iterations = 3000;

const float rho = 0.5;
const float t_pow = 1;
const float e_pow = 1;

const int num_nests = 50;
const float pa = 0.25;

const float beta = 0.5;
const float alpha = 1;

const float sigma_p = pow((tgamma(1 + beta)*sin(M_PI*beta / 2)) / (tgamma((1 + beta) / 2)*beta*pow(2, (beta - 1) / 2)), 2 / beta);

int nests[num_nests][num_vertices];
int fitness[num_nests];

mt19937 seed;
uniform_real_distribution<float> uni(0, 1);
uniform_int_distribution<int> random_vertex(0, num_vertices - 1);
uniform_int_distribution<int> random_colour(0, k - 1);
normal_distribution<float> normal_q(0, 1);
normal_distribution<float> normal_p(0, sigma_p);

float tau[num_vertices][num_vertices];
float d_tau[num_vertices][num_vertices];

void make_graph(float edge_probability) {  // Populates adj_matrix with a random graph
	for (int i = 0; i < num_vertices; i++) {
		adj_list_length[i] = 0;
		for (int j = 0; j < i; j++) {
			if (uni(seed) < edge_probability) {
				adj_list[i][adj_list_length[i]] = j;
				adj_list[j][adj_list_length[j]] = i;
				adj_list_length[i]++;
				adj_list_length[j]++;
				adj_matrix[i][j] = 1;
				adj_matrix[j][i] = 1;
			}
		}
	}
}

void get_cuckoo(int * nest) {
	for (int i = 0; i < num_vertices; i++) {
		nest[i] = random_colour(seed);
	}
}

void initialise_pheromones() {
	for (int i = 0; i < num_vertices; i++) {
		for (int j = 0; j < num_vertices; j++) {
			tau[i][j] = 1 - adj_matrix[i][j];
		}
	}
}

int colour_class_size[k];
int edge_conflicts[k];
int f(int * x) {
	for (int i = 0; i < k; i++) {
		colour_class_size[i] = 0;
		edge_conflicts[i] = 0;
	}
	for (int i = 0; i < num_vertices; i++) {
		colour_class_size[x[i]]++;
		for (int j = 0; j < adj_list_length[i]; j++) {
			if (x[i] == x[adj_list[i][j]]) {
				edge_conflicts[x[i]]++;
			}
		}
	}
	int ret = 0;
	for (int i = 0; i < k; i++) {
		ret += colour_class_size[i] * edge_conflicts[i] - colour_class_size[i] * colour_class_size[i];
	}
	return ret;
}

float levy() {
	float p = normal_p(seed);
	float q = normal_q(seed);
	float M = p / pow(abs(q), 1 / beta);
	return p / pow(abs(q), 1 / beta);
}

int num_conflicts(int * nest) {
	int num = 0;
	for (int i = 0; i < num_vertices; i++) {
		for (int j = 0; j < i; j++) {
			if (adj_matrix[i][j] == 1 && nest[i] == nest[j]) {
				num++;
			}
		}
	}
	return num;
}

bool found[k];
int neighbouring_colours[k];
float eta(int * nest, int v) {
	// Returns the heuristic value of v in nest
	int num_neighbouring = 0;
	for (int i = 0; i < k; i++) {
		found[i] = false;
	}
	for (int i = 0; i < adj_list_length[v]; i++) {
		if (!found[nest[adj_list[v][i]]]) {
			found[nest[adj_list[v][i]]] = true;
			num_neighbouring++;
		}
	}
	return num_neighbouring;
}

float weight[num_vertices];
int vertices[num_vertices];
int colour_counts[k];
bool tabu[num_vertices];
int e[num_vertices];
int levy_flight(int * nest, int start, int nest_fitness) {
	for (int i = 0; i < num_vertices; i++) {
		tabu[i] = false;
		e[i] = eta(nest, i);
	}
	int u = start;
	float M = abs(alpha * levy()) + 1;
	if (M > num_vertices) {
		M = num_vertices;
	}
	int v;
	for (int i = 0; i < M; i++) {
		vertices[i] = u;
		// Select a vertex v
		float weight_sum = 0;
		for (int w = 0; w < num_vertices; w++) {
			if (tabu[w]) {
				weight[w] = 0;
			}else {
				weight[w] = pow(tau[u][w], t_pow) * pow(eta(nest, w), e_pow);
				weight_sum += weight[w];
			}
		}
		float r = uni(seed);
		v = -1;
		do {
			v++;
			r -= weight[v] / weight_sum;
		} while (r > 0);
		// Recolour v to colour causing fewest conflicts
		int c = 0;
		for (int j = 0; j < k; j++) {
			colour_counts[j] = 0;
		}
		for (int j = 0; j < adj_list_length[v]; j++) {
			colour_counts[nest[adj_list[v][j]]]++;
		}
		for (int j = 1; j < k; j++) {
			if (colour_counts[j] < colour_counts[c]) {
				c = j;
			}
		}
		// Print the colourings somewhere around here to see what's going wrong
		nest_fitness += colour_counts[c] - colour_counts[nest[v]];
		nest[v] = c;
		tabu[v] = true;
		u = v;
	}
	for (int i = 0; i < M - 1; i++) {
		d_tau[vertices[i]][vertices[i + 1]] += 1 / (nest_fitness + 1);
	}
	return nest_fitness;
}

int main(){
	make_graph(0.5);
	for (int i = 0; i < num_nests; i++) {
		get_cuckoo(nests[i]);
		fitness[i] = num_conflicts(nests[i]);
	}
	int u = 0;
	for (int i = 1; i < num_vertices; i++) {
		if (adj_list_length[i] > adj_list_length[u]) {
			u = i;
		}
	}
	initialise_pheromones();
	int nest_temp[num_vertices];
	auto start = chrono::high_resolution_clock::now();
	for (int t = 0; t < num_iterations; t++) {
		// Reset d_tau
		for (int i = 0; i < num_vertices; i++) {
			for (int j = 0; j < num_vertices; j++) {
				d_tau[i][j] = 0;
			}
		}
		for (int c = 0; c < num_nests; c++) {
			copy(begin(nests[c]), end(nests[c]), begin(nest_temp));
			int l_f = levy_flight(nest_temp, u, fitness[c]);
			if (uni(seed) < pa || l_f < fitness[c]) {
				copy(begin(nest_temp), end(nest_temp), begin(nests[c]));
				fitness[c] = l_f;
				if (l_f == 0) {
					cout << "Iterations: " << t << endl;
					t = num_iterations;
					break;
				}
			}
		}
		for (int i = 0; i < num_vertices; i++) {
			for (int j = 0; j < num_vertices; j++) {
				tau[i][j] = (1 - rho) * tau[i][j] + d_tau[i][j];
			}
		}
		//cout << chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start).count() << endl;
	}
	int best = 0;
	int best_n = fitness[0];
	for (int i = 1; i < num_nests; i++) {
		int n = fitness[i];
		if (n < best_n) {
			best = i;
			best_n = n;
		}
	}
	for (int i = 0; i < num_vertices; i++) {
		cout << nests[best][i] << " ";
	}
	cout << endl << "Number of conflicts: " << best_n;
	cout << endl << "Time taken (seconds): " << chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start).count() / (float)1000000 << endl;
	return 0;
}