#define _USE_MATH_DEFINES
#define _SECURE_SCL 0

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <array>
#include <math.h>
#include <chrono>
#include <algorithm>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>

// Preliminary results (these are with running the program for 1 minute so take them with a grain of salt) (need to rerun for new eta precomputation method):
	// 100 vertices: 19 colours (0 conflicts)
	// 250 vertices: 39 colours (0 conflicts)
	// 500 vertices: 70 colours (0 conflicts)

// On DIMACS graphs (again, only running for 2 minutes):
	// latin_square (900 vertices, best known = 97): 137 colours

using namespace std;
using namespace boost::random;

const string graph_directory = "C:/Users/taydo/OneDrive/Documents/computer_science/year3/project/implementations/graphs/";

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
const auto duration = chrono::minutes{2};

const float rho = 0.5;
const float t_pow = 1;
const float e_pow = 1;
const int w = 5;

const int num_nests = 50;
const float pa = 0.1;
const float p = 0.1;

const float beta = 0.5;
const float alpha = 1;

const float sigma_p = pow((tgamma(1 + beta)*sin(M_PI*beta / 2)) / (tgamma((1 + beta) / 2)*beta*pow(2, (beta - 1) / 2)), 2 / beta);

struct Nest {
	int nest[num_vertices];
	int path[num_vertices];
	int eta[num_vertices];
	int path_length;
	int fitness;
};

Nest nests[num_nests];

/*int nests[num_nests][num_vertices];
int nest_path[num_nests][num_vertices];
int nest_path_length[num_nests];
int fitness[num_nests];*/

int colouring[num_vertices];

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

void read_graph(string filename) {
	string line;
	ifstream file;
	int u; int v;
	file.open(graph_directory + filename);
	if (file.is_open()) {
		while (getline(file, line)) {
			stringstream line_stream(line);
			line_stream >> line;
			if (line == "e") {
				line_stream >> u; line_stream >> v;
				u--; v--;
				adj_matrix[u][v] = 1; adj_matrix[v][u] = 1;
				adj_list[u][adj_list_length[u]] = v; adj_list[v][adj_list_length[v]] = u;
				adj_list_length[u]++; adj_list_length[v]++;
			}
		}
	}else {
		cout << "Couldn't open file." << endl;
		exit(1);
	}
	file.close();
}

bool found[num_vertices];
int num_colours(int * x) {
	int num = 0;
	for (int i = 0; i < num_vertices; i++) {
		found[i] = false;
	}
	for (int i = 0; i < num_vertices; i++) {
		if (!found[x[i]]) {
			found[x[i]] = true;
			num++;
		}
	}
	return num;
	/*vector<int> colours_used;
	for (int i = 0; i < num_vertices; i++) {
		if (find(colours_used.begin(), colours_used.end(), x[i]) == colours_used.end()) {
			colours_used.push_back(x[i]);
		}
	}
	return colours_used.size();
	int max = 0;
	for (int i = 0; i < num_vertices; i++) {
		if (x[i] > max) {
			max = x[i];
		}
	}
	return max + 1;*/
}

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
void get_cuckoo(int * nest, int * e) {
	for (int i = 0; i < num_vertices; i++) {
		nest[i] = random_colour(seed);
	}
	for (int i = 0; i < num_vertices; i++) {
		e[i] = eta(nest, i);
	}
}

void initialise_pheromones() {
	for (int i = 0; i < num_vertices; i++) {
		for (int j = 0; j < num_vertices; j++) {
			tau[i][j] = 1 - adj_matrix[i][j];
		}
	}
}

float sigmoid(float x) {
	return 1 / (1 + exp(-x));
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

bool is_legal(int * nest) {
	for (int i = 0; i < num_vertices; i++) {
		for (int j = 0; j < i; j++) {
			if (adj_matrix[i][j] == 1 && nest[i] == nest[j]) {
				return false;
			}
		}
	}
	return true;
}

bool tabu[num_vertices];
float levy_flight(int * nest, int * e, int start, int index) {
	for (int i = 0; i < num_vertices; i++) {
		tabu[i] = false;
	}
	int v = start;
	float M = abs(alpha * levy()) + 1;
	if (M > num_vertices) {
		M = num_vertices;
	}
	nests[index].path_length = M;
	int u;
	for (int i = 0; i < M; i++) {
		// Add v to the path and recolour it
		nests[index].path[i] = v;
		int c = 0;
		for (int j = 0; j < k; j++) {
			colour_counts[j] = 0;
		}
		for (int j = 0; j < adj_list_length[v]; j++) {
			colour_counts[nest[adj_list[v][j]]]++;
		}
		for (int j = 0; j < k; j++) {
			if (colour_counts[j] == 0) {
				c = j;
				break;
			}
			if (colour_counts[j] < colour_counts[c]) {
				c = j;
			}
		}
		if (nest[v] != c) {
			for (int j = 0; j < adj_list_length[v]; j++) {
				bool found_v = false; bool found_c = false;
				int neighbour = adj_list[v][j];
				for (int l = 0; l < adj_list_length[neighbour]; l++) {
					if (adj_list[neighbour][l] != v && nest[adj_list[neighbour][l]] == nest[v]) {
						found_v = true;
						break;
					}
				}
				for (int l = 0; l < adj_list_length[neighbour]; l++) {
					if (adj_list[neighbour][l] != v && nest[adj_list[neighbour][l]] == c) {
						found_c = true;
						break;
					}
				}
				e[neighbour] += found_v - found_c;
			}
			nest[v] = c;
		}
		tabu[v] = true;
		// Select a vertex u to colour next
		float weight_sum = 0;
		for (int w = 0; w < num_vertices; w++) {
			if (tabu[w]) {
				weight[w] = 0;
			}
			else {
				weight[w] = pow(tau[v][w], t_pow) * pow(e[w], e_pow);
				weight_sum += weight[w];
			}
		}
		float r = uni(seed);
		u = -1;
		do {
			u++;
			r -= weight[u] / weight_sum;
		} while (r > 0);
		// Recolour v to colour causing fewest conflicts
		v = u;
	}
	int fitness = f(nest);
	for (int i = 0; i < M - 1; i++) {
		d_tau[nests[index].path[i]][nests[index].path[i + 1]] += 1 - sigmoid(fitness);
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
					ret = c;
				}
				break;
			}
			c++;
		}
	}
	return ret + 1;
}

bool compare_nests(Nest & nest1, Nest & nest2) {
	return nest1.fitness < nest2.fitness;
}

int main() {
	read_graph("dsjc250.5.col");
	//make_graph(0.5);
	int u = 0;
	for (int i = 1; i < num_vertices; i++) {
		if (adj_list_length[i] > adj_list_length[u]) {
			u = i;
		}
	}
	k = chromatic_bound();
	random_colour = uniform_int_distribution<int>(0, k - 1);
	for (int i = 0; i < num_nests; i++) {
		get_cuckoo(nests[i].nest, nests[i].eta);
		nests[i].fitness = f(nests[i].nest);
	}
	initialise_pheromones();
	int nest_temp[num_vertices];
	int eta_temp[num_vertices];
	auto start = chrono::high_resolution_clock::now();
	int t = 0;
	while (chrono::duration_cast<chrono::minutes>(chrono::high_resolution_clock::now() - start) < duration) {
		t++;
		//auto start2 = chrono::high_resolution_clock::now();
		// Reset d_tau
		for (int i = 0; i < num_vertices; i++) {
			for (int j = 0; j < num_vertices; j++) {
				d_tau[i][j] = 0;
			}
		}
		for (int c = 0; c < num_nests; c++) {
			copy(begin(nests[c].nest), end(nests[c].nest), begin(nest_temp));
			copy(begin(nests[c].eta), end(nests[c].eta), begin(eta_temp));
			int l_f = levy_flight(nest_temp, eta_temp, u, c);
			if (uni(seed) < pa || l_f < nests[c].fitness) {
				int n = num_colours(nest_temp);
				if (n <= k && is_legal(nest_temp)) {
					copy(begin(nest_temp), end(nest_temp), begin(colouring));
					k = n - 1;
					random_colour = uniform_int_distribution<int>(0, k - 1);
				}
				copy(begin(nest_temp), end(nest_temp), begin(nests[c].nest));
				copy(begin(eta_temp), end(eta_temp), begin(nests[c].eta));
				nests[c].fitness = l_f;
			}
		}
		/*int best = 0;
		for (int i = 1; i < num_nests; i++) {
			if (nests[i].fitness < nests[best].fitness) {
				best = i;
			}
		}*/
		sort(begin(nests), end(nests), compare_nests);
		for (int i = 0; i < num_vertices; i++) {
			for (int j = 0; j < num_vertices; j++) {
				tau[i][j] = (1 - rho) * tau[i][j] + d_tau[i][j];
			}
		}
		for (int i = 0; i < nests[0].path_length - 1; i++) {
			tau[nests[0].path[i]][nests[0].path[i + 1]] += w * (1 - sigmoid(nests[0].fitness));
		}
		for (int i = 0; i < num_nests * p; i++) {
			get_cuckoo(nests[num_nests - i - 1].nest, nests[num_nests - i - 1].eta);
			nests[num_nests - i - 1].fitness = f(nests[num_nests - i - 1].nest);
		}
		//cout << chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start2).count() << endl;
	}
	for (int i = 0; i < num_vertices; i++) {
		cout << colouring[i] << " ";
	}
	cout << endl << "Number of colours: " << num_colours(colouring);
	cout << endl << "Number of conflicts: " << num_conflicts(colouring);
	cout << endl << "Number of iterations: " << t;
	return 0;
}