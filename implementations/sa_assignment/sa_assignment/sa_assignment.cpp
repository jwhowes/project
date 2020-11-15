// sa_assignment.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <array>
#include <chrono>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>

using namespace std;
using namespace boost::random;

const int NUM_VERTICES = 100;

mt19937 seed;
uniform_real_distribution<float> uni(0, 1);
uniform_int_distribution<int> random_vertex(0, NUM_VERTICES - 1);

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

void make_graph(float edge_probability) {
	for (int i = 0; i < NUM_VERTICES; i++) {
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

const int freeze_lim = 5;
const int size_factor = 16;
const float cutoff = 0.1f;
const float min_percent = 0.2f;

float T = 10000;
const float beta = 1.0005f;

int s[NUM_VERTICES];
int c;

int neighbour[NUM_VERTICES];
int neighbour_c;

int best_s[NUM_VERTICES];
int best_c;

int num_colours(int * col) {
	int max = 0;
	for (int i = 0; i < NUM_VERTICES; i++) {
		if (col[i] > max) {
			max = col[i];
		}
	}
	return max + 1;
}

int classes[NUM_VERTICES];  // Are we completely confident this is bounded by NUM_VERTICES?
int num_classes;
int f(int * col) {
	num_classes = num_colours(col);
	int ret = 0;
	for (int i = 0; i < NUM_VERTICES; i++) {
		classes[col[i]]++;
	}
	for (int i = 0; i < num_classes; i++) {
		ret -= classes[i] * classes[i];
	}
	return ret;
}

bool valid(int v, int c, int * col) {
	for (int i = 0; i < NUM_VERTICES; i++) {
		if (adj_matrix[v][i] == 1 && col[i] == c) {
			return false;
		}
	}
	return true;
}

int order[NUM_VERTICES];
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
		for (int j = 0; j < NUM_VERTICES; j++) {
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
	make_graph(0.5);
	for (int i = 0; i < NUM_VERTICES; i++) {
		order[i] = i;
	}
	auto start = chrono::high_resolution_clock::now();
	generate_initial_solution();
	c = f(s);
	copy(begin(s), end(s), begin(best_s));
	best_c = c;
	int freeze_count = 0;
	while (freeze_count < freeze_lim) {
		int changes = 0; int trials = 0;
		bool best_changed = false;
		while (trials < size_factor * num_colours(s) * NUM_VERTICES && changes < cutoff * num_colours(s) * NUM_VERTICES) {  // Maybe make sure the float multiplication is returning a float
			trials++;
			copy(begin(s), end(s), begin(neighbour));
			get_neighbour();
			neighbour_c = f(neighbour);
			int d = neighbour_c - c;
			if (c < 0) {
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
	for (int i = 0; i < NUM_VERTICES; i++) {
		cout << best_s[i] << " ";
	}
	cout << endl << "Num colours: " << num_colours(best_s) << endl;
	cout << "Time taken: " << chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start).count();
	return 0;
}
