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

// Preliminary results (these are with running the program for 1 minute so take them with a grain of salt):
	// 100 vertices: 19 colours (0 conflicts)
	// 250 vertices: 40 colours (0 conflicts)
	// 500 vertices: 70 colours (0 conflicts)
	

using namespace std;
using namespace boost::random;

const int num_vertices = 250;
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
int k;

const int num_iterations = 100;
const auto duration = chrono::minutes{1};

const float rho = 0.5;
const float t_pow = 1;
const float e_pow = 1;
const int w = 5;

const int num_nests = 50;
const float pa = 0.1;

const float beta = 0.5;
const float alpha = 1;

const float sigma_p = pow((tgamma(1 + beta)*sin(M_PI*beta / 2)) / (tgamma((1 + beta) / 2)*beta*pow(2, (beta - 1) / 2)), 2 / beta);

int nests[num_nests][num_vertices];
int nest_path[num_nests][num_vertices];
int nest_path_length[num_nests];
int fitness[num_nests];

mt19937 seed;
uniform_real_distribution<float> uni(0, 1);
uniform_int_distribution<int> random_vertex(0, num_vertices - 1);
uniform_int_distribution<int> random_colour;
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

float recip_sigmoid(float x) {
	return 1 + exp(-0.01*x);
}

int num_colours(int * x) {
	/*vector<int> colours_used;
	for (int i = 0; i < num_vertices; i++) {
		if (find(colours_used.begin(), colours_used.end(), x[i]) == colours_used.end()) {
			colours_used.push_back(x[i]);
		}
	}
	return colours_used.size();*/
	int max = 0;
	for (int i = 0; i < num_vertices; i++) {
		if (x[i] > max) {
			max = x[i];
		}
	}
	return max + 1;
}

int colour_class_size[num_vertices];
int edge_conflicts[num_vertices];
int f(int * x) {
	int n = num_colours(x);
	for (int i = 0; i < n; i++) {
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
	for (int i = 0; i < n; i++) {
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

bool found[num_vertices];
int neighbouring_colours[num_vertices];
float eta(int * nest, int v) {
	// Returns the heuristic value of v in nest
	int n = num_colours(nest);
	int num_neighbouring = 0;
	for (int i = 0; i < n; i++) {
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
int colour_counts[num_vertices];
bool tabu[num_vertices];
int e[num_vertices];
float levy_flight(int * nest, int start, int index) {
	for (int i = 0; i < num_vertices; i++) {
		tabu[i] = false;
		e[i] = eta(nest, i);
	}
	int u = start;
	float M = abs(alpha * levy()) + 1;
	if (M > num_vertices) {
		M = num_vertices;
	}
	nest_path_length[index] = M;
	int v;
	for (int i = 0; i < M; i++) {
		nest_path[index][i] = u;
		// Select a vertex v
		float weight_sum = 0;
		for (int w = 0; w < num_vertices; w++) {
			if (tabu[w]) {
				weight[w] = 0;
			}
			else {
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
		nest[v] = c;
		tabu[v] = true;
		u = v;
	}
	int fitness = f(nest);
	for (int i = 0; i < M - 1; i++) {
		d_tau[nest_path[index][i]][nest_path[index][i + 1]] += 1 / recip_sigmoid(fitness);
	}
	int n = num_colours(nest);
	if (n < k && num_conflicts(nest) == 0) {
		k = n;
		random_colour = uniform_int_distribution<int>(0, k - 1);
	}
	return fitness;
}

bool valid(int v, int c, int * col) {  // Returns whether or not vertex v can be coloured colour c in colouring col (legally)
	for (int i = 0; i < adj_list_length[v]; i++) {
		if (col[adj_list[v][i]] == c) {  // If v is adjacent to some vertex i coloured c then this is not a valid assignment
			return false;
		}
	}
	return true;  // If no such i can be found then the assignment is valid
}

int chromatic_bound() {
	int colouring[num_vertices];
	int ret = 0;
	for (int i = 0; i < num_vertices; i++) {
		int c = 0;
		while (true) {
			if (valid(i, c, colouring)) {
				colouring[i] = c;
				if (c > ret) {
					ret = i;
				}
				break;
			}
			c++;
		}
	}
	return ret;
}

int main() {
	make_graph(0.5);
	k = chromatic_bound();
	random_colour = uniform_int_distribution<int>(0, k - 1);
	for (int i = 0; i < num_nests; i++) {
		get_cuckoo(nests[i]);
		fitness[i] = f(nests[i]);
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
	while(chrono::duration_cast<chrono::minutes>(chrono::high_resolution_clock::now() - start) < duration){
		// Reset d_tau
		for (int i = 0; i < num_vertices; i++) {
			for (int j = 0; j < num_vertices; j++) {
				d_tau[i][j] = 0;
			}
		}
		for (int c = 0; c < num_nests; c++) {
			copy(begin(nests[c]), end(nests[c]), begin(nest_temp));
			int l_f = levy_flight(nest_temp, u, c);
			if (uni(seed) < pa || l_f < fitness[c]) {
				copy(begin(nest_temp), end(nest_temp), begin(nests[c]));
				fitness[c] = l_f;
			}
		}
		int best = 0;
		for (int i = 1; i < num_nests; i++) {
			if (fitness[i] < fitness[best]) {
				best = i;
			}
		}
		for (int i = 0; i < num_vertices; i++) {
			for (int j = 0; j < num_vertices; j++) {
				tau[i][j] = (1 - rho) * tau[i][j] + d_tau[i][j];
			}
		}
		for (int i = 0; i < nest_path_length[best] - 1; i++) {
			tau[nest_path[best][i]][nest_path[best][i + 1]] += w * 1 / recip_sigmoid(fitness[best]);
		}
		//cout << chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start).count() << endl;
	}
	int best = 0;
	for (int i = 1; i < num_nests; i++) {
		if (fitness[i] < fitness[best]) {
			best = i;
		}
	}
	for (int i = 0; i < num_vertices; i++) {
		cout << nests[best][i] << " ";
	}
	cout << endl << "Number of colours: " << num_colours(nests[best]);
	cout << endl << "Number of conflicts: " << num_conflicts(nests[best]);
	cout << endl << "Time taken (seconds): " << chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start).count() / (float)1000000 << endl;
	return 0;
}