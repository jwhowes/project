#define _USE_MATH_DEFINES
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
#include <boost/random/normal_distribution.hpp>

using namespace std;
using namespace boost::random;

const int num_vertices = 10;
int adj_matrix[num_vertices][num_vertices] = {
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
const int k = 3;

const int num_iterations = 1000;

const float rho = 0.5;
const float t_pow = 1;
const float e_pow = 1;

const int num_nests = 50;
const float pa = 0.25;

const float beta = 1.5;
const float alpha = 1.0;

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
		for (int j = 0; j < i; j++) {
			if (uni(seed) < edge_probability) {
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

float levy() {
	float p = normal_p(seed);
	float q = normal_q(seed);
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

int neighbouring_colours[num_vertices];
float eta(int * nest, int v) {
	// Returns the heuristic value of v in nest
	int num_neighbouring = 0;
	for (int i = 0; i < num_vertices; i++) {
		if (adj_matrix[v][i] == 1 && nest[i] != -1 && find(begin(neighbouring_colours), begin(neighbouring_colours) + num_neighbouring, nest[i]) == begin(neighbouring_colours) + num_neighbouring) {
			neighbouring_colours[num_neighbouring] = nest[i];
			num_neighbouring++;
		}
	}
	return num_neighbouring;
}

int weight[num_vertices];
int vertices[num_vertices];
void levy_flight(int * nest) {
	int u = random_vertex(seed);
	nest[u] = random_colour(seed);
	int M = alpha * levy();
	if (M > num_vertices) {
		M = num_vertices;
	}
	int v;
	for (int i = 0; i < M; i++) {
		vertices[i] = u;
		// Select a vertex v
		float weight_sum = 0;
		for (int w = 0; w < num_vertices; w++) {
			weight[w] = pow(tau[u][w], t_pow) * pow(eta(nest, w), e_pow) + 1;
			weight_sum += weight[w];
		}
		float r = uni(seed);
		v = -1;
		do {
			v++;
			r -= weight[v] / weight_sum;
		} while (r > 0);
		// Recolour v to u's colour
		nest[v] = nest[u];
		u = v;
	}
	for (int i = 0; i < M - 1; i++) {
		d_tau[vertices[i]][vertices[i + 1]] += 1 / (num_conflicts(nest) + 1);
	}
}

int main(){
	//make_graph(0.5);
	for (int i = 0; i < num_nests; i++) {
		get_cuckoo(nests[i]);
	}
	initialise_pheromones();
	int nest_temp[num_vertices];
	for (int t = 0; t < num_iterations; t++) {
		// Reset d_tau
		for (int i = 0; i < num_vertices; i++) {
			for (int j = 0; j < num_vertices; j++) {
				d_tau[i][j] = 0;
			}
		}
		for (int c = 0; c < num_nests; c++) {
			copy(begin(nests[c]), end(nests[c]), begin(nest_temp));
			levy_flight(nest_temp);
			if (uni(seed) < pa || num_conflicts(nest_temp) < num_conflicts(nests[c])) {
				copy(begin(nest_temp), end(nest_temp), begin(nests[c]));
			}
		}
		for (int i = 0; i < num_vertices; i++) {
			for (int j = 0; j < num_vertices; j++) {
				tau[i][j] = (1 - rho) * tau[i][j] + d_tau[i][j];
			}
		}
	}
	int best = 0;
	int best_n = num_conflicts(nests[0]);
	for (int i = 1; i < num_nests; i++) {
		int n = num_conflicts(nests[i]);
		if (n < best_n) {
			best = i;
			best_n = n;
		}
	}
	for (int i = 0; i < num_vertices; i++) {
		cout << nests[best][i] << " ";
	}
	cout << endl << "Number of conflicts: " << best_n;
	return 0;
}