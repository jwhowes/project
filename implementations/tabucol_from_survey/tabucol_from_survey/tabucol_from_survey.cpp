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

int k;
const int num_vertices = 300;

int adj_matrix[num_vertices][num_vertices];

int adj_list[num_vertices][num_vertices];/* = {
	{1, 4, 5, 0, 0, 0, 0, 0, 0, 0},
	{0, 2, 6, 0, 0, 0, 0, 0, 0, 0},
	{1, 3, 7, 0, 0, 0, 0, 0, 0, 0},
	{2, 4, 8, 0, 0, 0, 0, 0, 0, 0},
	{0, 3, 9, 0, 0, 0, 0, 0, 0, 0},
	{0, 7, 8, 0, 0, 0, 0, 0, 0, 0},
	{1, 8, 9, 0, 0, 0, 0, 0, 0, 0},
	{2, 5, 8, 0, 0, 0, 0, 0, 0, 0},
	{3, 5, 6, 0, 0, 0, 0, 0, 0, 0},
	{4, 6, 7, 0, 0, 0, 0, 0, 0, 0}
};*/
int adj_list_length[num_vertices];// = { 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 };

const float lambda = 0.6;
const int num_iterations = 3000;

chrono::time_point<chrono::steady_clock> start;
const auto duration = chrono::minutes{5};

int s[num_vertices];

int tabu_list[num_vertices][num_vertices];
int gamma[num_vertices][num_vertices];

mt19937 seed;
uniform_real_distribution<float> uni(0, 1);
uniform_int_distribution<int> random_vertex(0, num_vertices - 1);
uniform_int_distribution<int> random_L(0, 9);

void read_graph(string filename) {
	string line;
	ifstream file;
	int u; int v;
	file.open(graph_directory + filename);
	if (file.is_open()) {
		while (getline(file, line)) {
			char x;
			istringstream line_stream(line);
			line_stream >> x;
			if (x == 'e') {
				line_stream >> u >> v;
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

void populate_gamma() {
	for (int i = 0; i < num_vertices; i++) {
		for (int j = 0; j < k; j++) {
			gamma[i][j] = 0;
		}
	}
	for (int i = 0; i < num_vertices; i++) {
		for (int j = 0; j < adj_list_length[i]; j++) {
			gamma[i][s[adj_list[i][j]]]++;
		}
	}
}

int critical_vertices[num_vertices];
int num_critical;
void get_critical_vertices() {
	num_critical = 0;
	for (int i = 0; i < num_vertices; i++) {
		if (gamma[i][s[i]] > 0) {
			critical_vertices[num_critical] = i;
			num_critical++;
		}
	}
}

void update_gamma(int v, int c) {
	for (int i = 0; i < adj_list_length[v]; i++) {
		gamma[adj_list[v][i]][c]++;
		gamma[adj_list[v][i]][s[v]]--;
	}
}

int f(int * x) {
	int num = 0;
	for (int i = 0; i < num_vertices; i++) {
		for (int j = 0; j < i; j++) {
			if (adj_matrix[i][j] == 1 && x[i] == x[j]) {
				num++;
			}
		}
	}
	return num;
}

void make_move(int t) {
	bool initial = true;
	int best_v; int best_c; int best_d;
	get_critical_vertices();
	for (int i = 0; i < num_critical; i++) {
		int v = critical_vertices[i];
		for (int c = 0; c < k; c++) {
			int d = gamma[v][c] - gamma[v][s[v]];
			if (d < 0) {
				update_gamma(v, c);
				s[v] = c;
				tabu_list[v][c] = t + random_L(seed) + lambda * num_critical;
				return;
			} else if (tabu_list[v][c] <= t && (initial || d < best_d)) {
				initial = false;
				best_d = d; best_v = v; best_c = c;
			}
		}
	}
	if (!initial) {
		update_gamma(best_v, best_c);
		s[best_v] = best_c;
		tabu_list[best_v][best_c] = t + random_L(seed) + lambda * num_critical;
	}
}

/*void make_move(int t) {
	bool initial = true;
	int best_v; int best_c; int best_f;
	get_critical_vertices();
	for (int i = 0; i < num_vertices; i++) {
		int v = i;
		for (int c = 0; c < k; c++) {
			int old_c = s[v];
			int f_s = f(s);
			s[v] = c;
			if (f(s) < f_s) {
				cout << v << " " << c << endl;
				tabu_list[v][c] = t + random_L(seed) + lambda * num_critical;
				return;
			}
			else if (tabu_list[v][c] <= t && (initial || f_s < best_f)) {
				best_f = f_s;
				best_v = v; best_c = c;
				initial = false;
			}
			s[v] = old_c;
		}
	}
	if (!initial) {
		cout << best_v << " " << best_c << endl;
		s[best_v] = best_c;
		tabu_list[best_v][best_c] = t + random_L(seed) + lambda * num_critical;
	}
}*/

bool valid(int v, int c, int * col) {  // Returns whether or not vertex v can be coloured colour c in colouring col (legally)
	for (int i = 0; i < adj_list_length[v]; i++) {
		if (col[adj_list[v][i]] == c) {  // If v is adjacent to some vertex i coloured c then this is not a valid assignment
			return false;
		}
	}
	return true;  // If no such i can be found then the assignment is valid
}

int colouring[num_vertices];
int chromatic_bound() {
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

int colour_counts[num_vertices];
void generate_initial_solution() {
	for (int i = 0; i < num_vertices; i++) {
		s[i] = -1;
	}
	for (int i = 0; i < num_vertices; i++) {
		for (int j = 0; j < num_vertices; j++) {
			colour_counts[j] = 0;
		}
		for (int j = 0; j < adj_list_length[i]; j++) {
			if (s[adj_list[i][j]] >= 0) {
				colour_counts[s[adj_list[i][j]]]++;
			}
		}
		int c = 0;
		for (int j = 1; j < k; j++) {
			if (colour_counts[j] < colour_counts[c]) {
				c = j;
			}
		}
		s[i] = c;
	}
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

int t;
int global_t;
//ofstream ofile;
bool find_colouring() {
	generate_initial_solution();
	populate_gamma();
	while (chrono::duration_cast<chrono::minutes>(chrono::high_resolution_clock::now() - start) < duration) {
	//while(global_t < num_iterations){
		make_move(t);
		//if (global_t % 10 == 0) {
		//	ofile << num_colours(colouring) << endl;
		//}
		t++;
		global_t++;
		if (f(s) == 0) {
			copy(begin(s), end(s), begin(colouring));
			return true;
		}
	}
	return false;
}

int main(){
	cout << "TABUCOL\n";
	read_graph("flat300_26.col");
	//ofile.open(results_directory + "flat300_26_tabucol.txt");
	k = chromatic_bound() - 1;
	bool found_colouring = true;
	start = chrono::high_resolution_clock::now();
	global_t = 0;
	while(found_colouring){
		t = 0;
		found_colouring = find_colouring();
		k = num_colours(colouring) - 1;
	}
	//ofile.close();
	for (int i = 0; i < num_vertices; i++) {
		cout << colouring[i] << " ";
	}
	cout << endl << "Number of conflicts: " << f(colouring) << endl;
	cout << "Number of colours: " << num_colours(colouring) << endl;
	cout << "Number of iterations: " << global_t << endl;
	return 0;
}