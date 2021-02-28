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

const int num_vertices = 450;
int adj_matrix[num_vertices][num_vertices];
int adj_list[num_vertices][num_vertices];
int adj_list_length[num_vertices];

int k;

const int num_iterations = 3000;
const auto duration = chrono::minutes{ 60 };

const int num_nests = 10;
const int eggs_per_nest = 5;
const int num_cuckoos = 5;

const float p = 0.2;

struct Egg {
	int col[num_vertices];
	int fitness;
};

Egg nests[num_nests][eggs_per_nest];
int nest_sum_fitness[num_nests];
Egg cuckoos[num_cuckoos];

mt19937 seed;
uniform_real_distribution<float> uni(0, 1);
uniform_int_distribution<int> random_vertex(0, num_vertices - 1);
uniform_int_distribution<int> random_nest(0, num_nests - 1);
uniform_int_distribution<int> random_egg(0, eggs_per_nest - 1);
uniform_int_distribution<int> random_colour;

int best_colouring[num_vertices];

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

void generate_new_egg(int * col) {
	for (int i = 0; i < num_vertices; i++) {
		col[i] = random_colour(seed);
	}
}

void generate_new_nest(int n) {
	nest_sum_fitness[n] = 0;
	for (int i = 0; i < eggs_per_nest; i++) {
		generate_new_egg(nests[n][i].col);
		nests[n][i].fitness = f(nests[n][i].col);
		nest_sum_fitness[n] += nests[n][i].fitness;
	}
}

int colour_counts[num_vertices];
int nest_colour_class_size[eggs_per_nest][num_vertices];
void gpx_nest(int * egg, int nest) {
	for (int i = 0; i < num_vertices; i++) {
		egg[i] = -1;
		for (int j = 0; j < eggs_per_nest; j++) {
			nest_colour_class_size[j][i] = 0;
		}
	}
	for (int i = 0; i < num_vertices; i++) {
		for (int j = 0; j < eggs_per_nest; j++) {
			nest_colour_class_size[j][nests[nest][j].col[i]]++;
		}
	}
	int m;
	for (int i = 0; i < k; i++) {
		int e = i % eggs_per_nest;
		m = 0;
		for (int j = 1; j < k; j++) {
			if (nest_colour_class_size[e][j] > nest_colour_class_size[e][m]) {
				m = j;
			}
		}
		for (int j = 0; j < num_vertices; j++) {
			if (egg[j] == -1 && nests[nest][e].col[j] == m) {
				egg[j] = i;
				for (int l = 0; l < eggs_per_nest; l++) {
					nest_colour_class_size[l][j]--;
				}
			}
		}
		nest_colour_class_size[e][m] = 0;
	}
	// Colour uncoloured vertices in egg randomly
	for (int i = 0; i < num_vertices; i++) {
		if (egg[i] == -1) {
			//egg[i] = random_colour(seed);
			int c = 0;
			for (int j = 0; j < k; j++) {
				colour_counts[j] = 0;
			}
			for (int j = 0; j < adj_list_length[i]; j++) {
				if (egg[adj_list[i][j]] != -1) {
					colour_counts[egg[adj_list[i][j]]]++;
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
			egg[i] = c;
		}
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

int colour_class_size_1[num_vertices];
int colour_class_size_2[num_vertices];
void gpx2(int * egg, int * p1, int * p2) {
	for (int i = 0; i < num_vertices; i++) {
		egg[i] = -1;
		colour_class_size_1[i] = 0;
		colour_class_size_2[i] = 0;
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
				if (egg[j] == -1 && p1[j] == m) {
					egg[j] = i;
					colour_class_size_2[p2[j]]--;
				}
			}
			colour_class_size_1[m] = 0;
		}
		else {
			m = 0;
			for (int j = 1; j < k; j++) {
				if (colour_class_size_2[j] > colour_class_size_2[m]) {
					m = j;
				}
			}
			for (int j = 0; j < num_vertices; j++) {
				if (egg[j] == -1 && p2[j] == m) {
					egg[j] = i;  // Could make it egg[j] = i instead? (To maintain diversity)
					colour_class_size_1[p1[j]]--;
				}
			}
			colour_class_size_2[m] = 0;
		}
	}
	// Colour uncoloured vertices in egg randomly
	for (int i = 0; i < num_vertices; i++) {
		if (egg[i] == -1) {
			//egg[i] = random_colour(seed);
			int c = 0;
			for (int j = 0; j < k; j++) {
				colour_counts[j] = 0;
			}
			for (int j = 0; j < adj_list_length[i]; j++) {
				if (egg[adj_list[i][j]] != -1) {
					colour_counts[egg[adj_list[i][j]]]++;
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
			egg[i] = c;
		}
	}
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

int critical_vertices[num_vertices];
int num_critical;
int tabucol_iterations;// = 100;
int tabu_list[num_vertices][num_vertices];
int gamma[num_vertices][num_vertices];
uniform_int_distribution<int> random_L(0, 9);
const float move_p = 1;
const float lambda = 0.6;

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
	tabucol_iterations = 10 * k * num_vertices / (move_p * (k - 1) * num_conflicts(colouring));
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

int max_colour(int * x) {
	int m = 0;
	for (int i = 0; i < num_vertices; i++) {
		if (x[i] > m) {
			m = x[i];
		}
	}
	return m;
}

int D[num_vertices][num_vertices];
int px[num_vertices][num_vertices];
int px_length[num_vertices];
int py[num_vertices][num_vertices];
int py_length[num_vertices];
int min_D[num_vertices];
int mapping[num_vertices];
int colours[num_vertices];
bool taken[num_vertices];

bool compare_classes(int c1, int c2) {
	return min_D[c1] < min_D[c2];
}

int d(int * x, int * y) {
	int k = max_colour(x) + 1;
	int k_2 = max_colour(y) + 1;
	if (k_2 > k) {
		k = k_2;
	}
	for (int i = 0; i < k; i++) {
		min_D[i] = -1;
		colours[i] = i;
		px_length[i] = 0;
		py_length[i] = 0;
	}
	for (int i = 0; i < num_vertices; i++) {
		px_length[x[i]]++;
		py_length[y[i]]++;
	}
	for (int i = 0; i < k; i++) {
		taken[i] = false;
		for (int j = 0; j < k; j++) {
			D[j][i] = px_length[i] + py_length[j];
		}
	}
	for (int i = 0; i < num_vertices; i++) {
		D[y[i]][x[i]] -= 2;
		if (min_D[x[i]] == -1 || D[y[i]][x[i]] < min_D[x[i]]) {
			min_D[x[i]] = D[y[i]][x[i]];
		}
	}
	sort(begin(colours), begin(colours) + k, compare_classes);
	for (int i = 0; i < k; i++) {
		int c = colours[i];
		int min = -1;
		for (int j = 0; j < k; j++) {
			if (!taken[j] && (min == -1 || D[c][j] < D[c][min])) {
				min = j;
			}
		}
		mapping[c] = min;
		taken[min] = true;
	}
	int d = 0;
	for (int i = 0; i < num_vertices; i++) {
		if (mapping[x[i]] != y[i]) {
			d++;
		}
	}
	return d;
}

int main(){
	cout << "NEW\n";
	ofstream ofile;
	ofile.open(results_directory + "le450_5a_new.txt");
	read_graph("le450_5a.col");
	int egg_temp[num_vertices];
	int nest_egg[num_vertices];
	int cuckoo_sum_distance[num_vertices];
	k = chromatic_bound();
	random_colour = uniform_int_distribution<int>(0, k - 1);
	// Generate nest eggs
	for (int i = 0; i < num_nests; i++) {
		generate_new_nest(i);
	}
	// Generate cuckoo eggs
	for (int i = 0; i < num_cuckoos; i++) {
		generate_new_egg(cuckoos[i].col);
		cuckoos[i].fitness = f(cuckoos[i].col);
	}
	for (int t = 0; t < num_iterations; t++) {
		// Each cukoo creates and lays an egg
		for (int c = 0; c < num_cuckoos; c++) {
			int n = random_nest(seed);
			//int e = random_egg(seed);
			gpx_nest(nest_egg, n);
			gpx2(egg_temp, nest_egg, cuckoos[c].col);
			//gpx2(egg_temp, cuckoos[c].col, nests[n][e].col);
			int egg_fitness = tabucol(egg_temp);
			int max_fitness = 0;
			for (int i = 1; i < eggs_per_nest; i++) {
				if (nests[n][i].fitness > nests[n][max_fitness].fitness) {
					max_fitness = i;
				}
			}
			nest_sum_fitness[n] += egg_fitness - nests[n][max_fitness].fitness;
			copy(begin(egg_temp), end(egg_temp), begin(nests[n][max_fitness].col));
			nests[n][max_fitness].fitness = egg_fitness;
			copy(begin(egg_temp), end(egg_temp), begin(cuckoos[c].col));
			cuckoos[c].fitness = egg_fitness;
			int num = num_colours(cuckoos[c].col);
			if (num <= k && num_conflicts(cuckoos[c].col) == 0) {
				copy(begin(cuckoos[c].col), end(cuckoos[c].col), begin(best_colouring));
				k = num - 1;
				random_colour = uniform_int_distribution<int>(0, k - 1);
			}
		}
		// Worst nest is abandoned
		/*int worst_nest = 0;
		for (int i = 1; i < num_nests; i++) {
			if (nest_sum_fitness[i] > nest_sum_fitness[worst_nest]) {
				worst_nest = i;
			}
		}
		generate_new_nest(worst_nest);*/
		// Least diverse cuckoo is killed
		if (uni(seed) < p) {
			int worst_cuckoo = 0;
			int worst_cuckoo_sd = -1;
			for (int i = 0; i < num_cuckoos; i++) {
				int sum_distance = 0;
				for (int j = 0; j < num_cuckoos; j++) {
					sum_distance += d(cuckoos[i].col, cuckoos[j].col);
				}
				if (worst_cuckoo_sd = -1 || sum_distance < worst_cuckoo_sd) {
					worst_cuckoo_sd = sum_distance;
					worst_cuckoo = i;
				}
			}
			generate_new_egg(cuckoos[worst_cuckoo].col);
			cuckoos[worst_cuckoo].fitness = f(cuckoos[worst_cuckoo].col);
		}
		if (t % 10 == 0) {
			ofile << num_colours(best_colouring) << endl;
		}
	}
	for (int i = 0; i < num_vertices; i++) {
		cout << best_colouring[i] << " ";
	}
	cout << endl << "Number of colours: " << num_colours(best_colouring) << endl;
	cout << "Number of conflicts: " << num_conflicts(best_colouring) << endl;
}