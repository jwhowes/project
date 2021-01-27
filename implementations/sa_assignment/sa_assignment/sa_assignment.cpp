#define _SECURE_SCL 0

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <array>
#include <chrono>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>

using namespace std;
using namespace boost::random;

const string graph_directory = "C:/Users/taydo/OneDrive/Documents/computer_science/year3/project/implementations/graphs/";

const int num_vertices = 250;

mt19937 seed;
uniform_real_distribution<float> uni(0, 1);
uniform_int_distribution<int> random_vertex(0, num_vertices - 1);

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

void make_graph(float edge_probability) {
	for (int i = 0; i < num_vertices; i++) {
		for (int j = 0; j < i; j++) {
			if (uni(seed) < edge_probability) {
				adj_matrix[i][j] = 1;
				adj_matrix[j][i] = 1;
				//edge_list[i].push_back(j);
				//edge_list[j].push_back(i);
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
			}
		}
	}
	else {
		cout << "Couldn't open file." << endl;
		exit(1);
	}
	file.close();
}

const int freeze_lim = 5;
const int size_factor = 16;
const float cutoff = 0.1f;
const float min_percent = 0.2f;

float T = 10000;
const float beta = 1.0005f;

int s[num_vertices];
int c;

int neighbour[num_vertices];
int neighbour_c;

int best_s[num_vertices];
int best_c;

int num_colours(int * col) {
	int max = 0;
	for (int i = 0; i < num_vertices; i++) {
		if (col[i] > max) {
			max = col[i];
		}
	}
	return max + 1;
}

int classes[num_vertices];  // Are we completely confident this is bounded by num_vertices?
int num_classes;
int f(int * col) {
	num_classes = num_colours(col);
	int ret = 0;
	for (int i = 0; i < num_vertices; i++) {
		classes[col[i]]++;
	}
	for (int i = 0; i < num_classes; i++) {
		ret -= classes[i] * classes[i];
	}
	return ret;
}

bool valid(int v, int c, int * col) {
	for (int i = 0; i < num_vertices; i++) {
		if (adj_matrix[v][i] == 1 && col[i] == c) {
			return false;
		}
	}
	return true;
}

int order[num_vertices];
void generate_initial_solution() {
	random_shuffle(begin(order), end(order));
	for (int v : order) {
		int c = 0;
		while (true) {
			if (valid(v, c, s)) {
				s[v] = c;
				break;
			}
			c++;
		}
	}
}

vector<int> kempe_chain(int c, int d, int v) {
	vector<int> K = { v };
	int i = 0;
	while (i < K.size()) {
		for (int j = 0; j < num_vertices; j++) {
			if (adj_matrix[K[i]][j] == 1 && (s[j] == c || s[j] == d) && find(K.begin(), K.end(), j) == K.end()) {
				K.push_back(j);
			}
		}
		i++;
	}
	return K;
}

void get_neighbour() {  // Copy s into neighbour first
	int v = random_vertex(seed);
	int c = s[v];
	int n = num_colours(s) - 1;
	int d = uniform_int_distribution<int>(0, n)(seed);
	if (d == c) {
		d = (d + 1) % n;
	}
	vector<int> K = kempe_chain(c, d, v);
	for (int i = 0; i < K.size(); i++) {
		if (s[K[i]] == c) {
			neighbour[K[i]] = d;
		}else {
			neighbour[K[i]] = c;
		}
	}
}

int main(){
	cout << "SA\n";
	//make_graph(0.5);
	read_graph("dsjc250.5.col");
	for (int i = 0; i < num_vertices; i++) {
		order[i] = i;
	}
	auto start = chrono::high_resolution_clock::now();
	const auto duration = chrono::minutes{5};
	generate_initial_solution();
	c = f(s);
	copy(begin(s), end(s), begin(best_s));
	best_c = c;
	int freeze_count = 0;
	int t = 0;
	while (chrono::duration_cast<chrono::minutes>(chrono::high_resolution_clock::now() - start) < duration/*freeze_count < freeze_lim*/) {
		t++;
		int changes = 0; int trials = 0;
		bool best_changed = false;
		while (trials < size_factor * num_colours(s) * num_vertices && changes < cutoff * num_colours(s) * num_vertices) {  // Maybe make sure the float multiplication is returning a float
			trials++;
			copy(begin(s), end(s), begin(neighbour));
			get_neighbour();
			neighbour_c = f(neighbour);
			int d = neighbour_c - c;
			if (d < 0) {
				changes++;
				copy(begin(neighbour), end(neighbour), begin(s));
				c = neighbour_c;
				if (c < best_c) {
					best_changed = true;
					copy(begin(s), end(s), begin(best_s));
					best_c = c;
				} else if (uni(seed) <= exp(-d / T)) {
					changes++;
					copy(begin(neighbour), end(neighbour), begin(s));
					c = neighbour_c;
				}
			}
		}
		T /= beta;
		if (best_changed) {
			freeze_count = 0;
		}
		if (changes / (float)trials < min_percent) {
			freeze_count++;
		}
	}
	for (int i = 0; i < num_vertices; i++) {
		std::cout << best_s[i] << " ";
	}
	std::cout << endl << "Number of colours: " << num_colours(best_s) << endl;
	std::cout << "Number of iterations: " << chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start).count();
	return 0;
}
