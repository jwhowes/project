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

using namespace std;
using namespace boost::random;

const string graph_directory = "C:/Users/taydo/OneDrive/Documents/computer_science/year3/project/implementations/graphs/";
const string results_directory = "C:/Users/taydo/OneDrive/Documents/computer_science/year3/project/implementations/results/";

const int num_vertices = 500;
int adj_matrix[num_vertices][num_vertices];
int adj_list[num_vertices][num_vertices];
int adj_list_length[num_vertices];

int k;

const int num_iterations = 3000;
const auto duration = chrono::minutes{60};

const int num_nests = 50;
const float pa = 0.1;
const float p = 0.1;

const float rho = 0.2;
const float t_pow = 1;
const float e_pow = 1;
const int w = 5;

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

struct Reference {
	int nest[num_vertices];
	int fitness;
	float diversity;
};

Nest nests[num_nests];

const int b = 5;
const int d = 5;
Reference ref_set[b + d];
int worst_best_set = 0;
int most_diverse_div_set = 0;

int best_colouring[num_vertices];

mt19937 seed;
uniform_real_distribution<float> uni(0, 1);
uniform_int_distribution<int> random_vertex(0, num_vertices - 1);
uniform_int_distribution<int> random_colour;
normal_distribution<float> normal_q(0, 1);
normal_distribution<float> normal_p(0, sigma_p);
uniform_int_distribution<int> random_ref_set(0, b + d - 1);

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
	}
	else {
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

bool valid(int v, int c, int * col) {  // Returns whether or not vertex v can be coloured colour c in colouring col (legally)
	for (int i = 0; i < adj_list_length[v]; i++) {
		if (col[adj_list[v][i]] == c) {  // If v is adjacent to some vertex i coloured c then this is not a valid assignment
			return false;
		}
	}
	return true;  // If no such i can be found then the assignment is valid
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

float weight[num_vertices];
int colour_counts[num_vertices];
void get_random_cuckoo(int * nest, int * e) {
	for (int i = 0; i < num_vertices; i++) {
		nest[i] = random_colour(seed);
	}
	for (int i = 0; i < num_vertices; i++) {
		e[i] = eta(nest, i);
	}
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

int pi[num_vertices][num_vertices];
int distance(int * x, int * y) {
	int k = num_colours(x);
	int d = -k;
	for (int i = 0; i < k; i++) {
		for (int j = 0; j < k; j++) {
			pi[i][j] = 0;
		}
	}
	for (int i = 0; i < num_vertices; i++) {
		if (pi[x[i]][y[i]] == 0) {
			d++;
		}
		pi[x[i]][y[i]]++;
	}
	return d;
}

/*int distance(int * col1, int * col2) {  // Hamming distance
	int num = 0;
	for (int i = 0; i < num_vertices; i++) {
		if (col1[i] != col2[i]) {
			num++;
		}
	}
	return num;
}*/

void gpx(int * nest, int * p1, int * p2) {
	int colour_class_size_1[num_vertices];
	int colour_class_size_2[num_vertices];
	for (int i = 0; i < num_vertices; i++) {
		nest[i] = -1;
		colour_class_size_1[i] = 0;
		colour_class_size_1[i] = 0;
	}
	for (int i = 0; i < num_vertices; i++) {
		colour_class_size_1[p1[i]]++;
		colour_class_size_2[p2[i]]++;
	}
	int m;
	for (int i = 0; i < k; i++) {
		if (i % 2 == 0) {
			m = 0;
			for (int j = 1; j < k; j++) {
				if (colour_class_size_1[j] > colour_class_size_1[m]) {
					m = j;
				}
			}
			for (int j = 0; j < num_vertices; j++) {
				if (p1[j] == m) {
					nest[j] = i;  // Could make it nest[j] = i instead? (To maintain diversity)
					colour_class_size_2[p2[j]]--;
					p2[j] = -1;
					p1[j] = -1;
				}
			}
			colour_class_size_1[m] = 0;
		} else {
			m = 0;
			for (int j = 1; j < k; j++) {
				if (colour_class_size_2[j] > colour_class_size_2[m]) {
					m = j;
				}
			}
			for (int j = 0; j < num_vertices; j++) {
				if (p2[j] == m) {
					nest[j] = i;  // Could make it nest[j] = i instead? (To maintain diversity)
					colour_class_size_1[p2[j]]--;
					p1[j] = -1;
					p2[j] = -1;
				}
			}
			colour_class_size_2[m] = 0;
		}
	}
	// Colour uncoloured vertices in nest randomly
	for (int i = 0; i < num_vertices; i++) {
		if (nest[i] == -1) {
			//nest[i] = random_colour(seed);
			int c = 0;
			for (int j = 0; j < k; j++) {
				colour_counts[j] = 0;
			}
			for (int j = 0; j < adj_list_length[i]; j++) {
				if (nest[adj_list[i][j]] != -1) {
					colour_counts[nest[adj_list[i][j]]]++;
				}
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
			nest[i] = c;
		}
	}
}

void get_cuckoo(int * nest, int * e) {
	// Generate new cuckoo from RefSet via GPX
	int p1 = random_ref_set(seed);
	int p2 = random_ref_set(seed);
	int p1_temp[num_vertices];
	int p2_temp[num_vertices];
	if (p1 == p2) {
		copy(begin(ref_set[p1].nest), end(ref_set[p1].nest), nest);
	} else {
		copy(begin(ref_set[p1].nest), end(ref_set[p1].nest), begin(p1_temp));
		copy(begin(ref_set[p2].nest), end(ref_set[p2].nest), begin(p2_temp));
		gpx(nest, p1_temp, p2_temp);
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

float levy() {
	float p = normal_p(seed);
	float q = normal_q(seed);
	float M = p / pow(abs(q), 1 / beta);
	return p / pow(abs(q), 1 / beta);
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
		} while (r > 0 && u < num_vertices - 1);
		// Recolour v to colour causing fewest conflicts
		v = u;
	}
	int fitness = f(nest);
	for (int i = 0; i < M - 1; i++) {
		d_tau[nests[index].path[i]][nests[index].path[i + 1]] += 1 - sigmoid(fitness);
	}
	return fitness;
}

/*int chromatic_bound() {
	int ret = 0;
	for (int i = 0; i < num_vertices; i++) {
		best_colouring[i] = -1;
	}
	for (int i = 0; i < num_vertices; i++) {
		int max = -1;
		for (int j = 0; j < num_vertices; j++) {
			if (best_colouring[j] == -1 && (max == -1 || eta(best_colouring, j) > eta(best_colouring, max))) {
				max = j;
			}
		}
		int c = 0;
		while (!valid(max, c, best_colouring)) {
			c++;
		}
		best_colouring[max] = c;
		if (c > ret) {
			ret = c;
		}
	}
	return ret + 1;
}*/

int chromatic_bound() {
	int ret = 0;
	for (int i = 0; i < num_vertices; i++) {
		int c = 0;
		while (true) {
			if (valid(i, c, best_colouring)) {
				best_colouring[i] = c;
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

int critical_vertices[num_vertices];
int num_critical;
int tabucol_iterations;// = 100;
int tabu_list[num_vertices][num_vertices];
int gamma[num_vertices][num_vertices];
uniform_int_distribution<int> random_L(0, 9);
const float lambda = 0.6;
const float move_p = 0.5;

void get_critical_vertices(int * colouring) {
	num_critical = 0;
	for (int i = 0; i < num_vertices; i++) {
		if (gamma[i][colouring[i]] > 0) {
			critical_vertices[num_critical] = i;
			num_critical++;
		}
	}
}

void inline update_gamma(int v, int c, int * colouring) {
	for (int i = 0; i < adj_list_length[v]; i++) {
		gamma[adj_list[v][i]][c]++;
		gamma[adj_list[v][i]][colouring[v]]--;
	}
}

void tabucol_make_move(int t, int * colouring) {
	bool initial = true;
	int best_v; int best_c; int best_d;
	get_critical_vertices(colouring);
	for (int i = 0; i < num_critical; i++) {
		if (uni(seed) <= move_p) {
			int v = critical_vertices[i];
			for (int c = 0; c < k; c++) {
				int d = gamma[v][c] - gamma[v][colouring[v]];
				if (d < 0) {
					update_gamma(v, c, colouring);
					colouring[v] = c;
					tabu_list[v][c] = t + random_L(seed) + lambda * num_critical;
					return;
				}
				else if (tabu_list[v][c] <= t && (initial || d < best_d)) {
					initial = false;
					best_d = d; best_v = v; best_c = c;
				}
			}
		}
	}
	if (!initial) {
		update_gamma(best_v, best_c, colouring);
		colouring[best_v] = best_c;
		tabu_list[best_v][best_c] = t + random_L(seed) + lambda * num_critical;
	}
}

int tabucol(int * colouring) {
	if (num_conflicts(colouring) == 0) {
		return f(colouring);
	}
	tabucol_iterations = k * num_vertices / (move_p * (k - 1) * num_conflicts(colouring));
	for (int i = 0; i < num_vertices; i++) {
		for (int j = 0; j < k; j++) {
			tabu_list[i][j] = 0;
			gamma[i][j] = 0;
		}
	}
	for (int i = 0; i < num_vertices; i++) {
		for (int j = 0; j < adj_list_length[i]; j++) {
			gamma[i][colouring[adj_list[i][j]]]++;
		}
	}
	for (int t = 0; t < tabucol_iterations; t++) {
		tabucol_make_move(t, colouring);
		if (num_conflicts(colouring) == 0) {
			return f(colouring);
		}
	}
	return f(colouring);
}

void update_ref_set() {
	most_diverse_div_set = 0;
	for (int i = 0; i < b + d; i++) {
		ref_set[i].diversity = 0;
		for (int j = 0; j < b + d; j++) {
			ref_set[i].diversity += distance(ref_set[i].nest, ref_set[j].nest);
		}
		ref_set[i].diversity /= b + d;
		if (ref_set[i].diversity > ref_set[most_diverse_div_set].diversity) {
			most_diverse_div_set = i;
		}
	}
	for (int i = 0; i < b; i++) {
		if (ref_set[i].fitness > ref_set[worst_best_set].fitness) {
			worst_best_set = i;
		}
	}
}

int main() {
	cout << "CSACO_refset\n";
	ofstream ofile;
	ofile.open(results_directory + "dsjc500.5_csaco_refset4.2.txt");
	read_graph("dsjc500.5.col");
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
		get_random_cuckoo(nests[i].nest, nests[i].eta);
		nests[i].fitness = f(nests[i].nest);
	}
	for (int i = 0; i < b; i++) {
		copy(begin(nests[i].nest), end(nests[i].nest), begin(ref_set[i].nest));
		ref_set[i].fitness = nests[i].fitness;
	}
	update_ref_set();
	initialise_pheromones();
	int nest_temp[num_vertices];
	int eta_temp[num_vertices];
	auto start = chrono::high_resolution_clock::now();
	//int t = 0;
	//while (chrono::duration_cast<chrono::minutes>(chrono::high_resolution_clock::now() - start) < duration) {
	for (int t = 0; t < num_iterations; t++) {
		//t++;
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
					copy(begin(nest_temp), end(nest_temp), begin(best_colouring));
					k = n - 1;
					random_colour = uniform_int_distribution<int>(0, k - 1);
				}
				copy(begin(nest_temp), end(nest_temp), begin(nests[c].nest));
				copy(begin(eta_temp), end(eta_temp), begin(nests[c].eta));
				nests[c].fitness = l_f;
			}
		}
		sort(begin(nests), end(nests), compare_nests);
		copy(begin(nests[0].nest), end(nests[0].nest), begin(nest_temp));
		int t_f = tabucol(nest_temp);
		if (t_f < nests[0].fitness) {
			copy(begin(nest_temp), end(nest_temp), begin(nests[0].nest));
			nests[0].fitness = t_f;
		}
		for (int i = 0; i < num_nests; i++) {
			if (nests[i].fitness < ref_set[worst_best_set].fitness) {
				copy(begin(nests[i].nest), end(nests[i].nest), begin(ref_set[worst_best_set].nest));
				ref_set[worst_best_set].fitness = nests[i].fitness;
				update_ref_set();
			} else {
				float div = 0;
				for (int j = 0; j < b + d; j++) {
					div += distance(nests[i].nest, ref_set[j].nest);
				}
				div /= b + d;
				if (div > ref_set[most_diverse_div_set].diversity) {
					copy(begin(nests[i].nest), end(nests[i].nest), begin(ref_set[most_diverse_div_set].nest));
					ref_set[most_diverse_div_set].fitness = nests[i].fitness;
					update_ref_set();
				}
			}
		}
		for (int i = 0; i < num_vertices; i++) {
			for (int j = 0; j < num_vertices; j++) {
				tau[i][j] = (1 - rho) * tau[i][j] + d_tau[i][j];
			}
		}
		// Apply additional pheromone to trail generated by the best nest
		for (int i = 0; i < nests[0].path_length - 1; i++) {
			tau[nests[0].path[i]][nests[0].path[i + 1]] += w * (1 - sigmoid(nests[0].fitness));
		}
		// Kill a fraction p of the worst nests
		for (int i = 0; i < num_nests * p; i++) {
			get_cuckoo(nests[num_nests - i - 1].nest, nests[num_nests - i - 1].eta);
			nests[num_nests - i - 1].fitness = f(nests[num_nests - i - 1].nest);
		}
		if (t % 10 == 0) {
			ofile << num_colours(best_colouring) << endl;
		}
	}
	ofile.close();
	for (int i = 0; i < num_vertices; i++) {
		cout << best_colouring[i] << " ";
	}
	cout << endl << "Number of colours: " << num_colours(best_colouring) << endl;
	cout << "Number of conflicts: " << num_conflicts(best_colouring) << endl;
	//cout << "Number of iterations: " << t << endl;
	cout << "Time taken (ms): " << chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - start).count() << endl;
	return 0;
}