#define _SECURE_SCL 0

#include <iostream>
#include <vector>
#include <array>
#include <chrono>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>

using namespace std;
using namespace boost::random;

const int k = 3;
const int NUM_VERTICES = 10;

int adj_list[NUM_VERTICES][NUM_VERTICES] = {
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
};
int adj_list_length[NUM_VERTICES] = { 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 };

const float lambda = 0.6;
const int num_iterations = 3000;

int s[NUM_VERTICES];

int tabu_list[NUM_VERTICES][k];
int gamma[NUM_VERTICES][k];

mt19937 seed;
uniform_real_distribution<float> uni(0, 1);
uniform_int_distribution<int> random_colour(0, k - 1);
uniform_int_distribution<int> random_vertex(0, NUM_VERTICES - 1);
uniform_int_distribution<int> random_L(0, 9);

void populate_gamma() {
	for (int i = 0; i < NUM_VERTICES; i++) {
		for (int j = 0; j < adj_list_length[i]; j++) {
			gamma[i][s[j]]++;
		}
	}
}

vector<int> critical_vertices;
void get_critical_vertices() {
	critical_vertices.clear();
	for (int i = 0; i < NUM_VERTICES; i++) {
		if (gamma[i][s[i]] > 0) {
			critical_vertices.push_back(i);
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
	for (int i = 0; i < NUM_VERTICES; i++) {
		for (int j = i; j < adj_list_length[i]; j++) {
			if (s[adj_list[i][j]] == s[i]) {
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
	for (int v : critical_vertices) {
		for (int c = 0; c < k; c++) {
			int d = gamma[v][c] - gamma[v][s[v]];
			if (d < 0) {
				update_gamma(v, c);
				cout << f(s) << " ";
				s[v] = c;
				cout << f(s) << endl;
				tabu_list[v][c] = t + random_L(seed) + lambda * critical_vertices.size();
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
		tabu_list[best_v][best_c] = t + random_L(seed) + lambda * critical_vertices.size();
	}
}

void generate_initial_solution() {
	for (int i = 0; i < NUM_VERTICES; i++) {
		s[i] = random_colour(seed);
	}
}

int main(){
	generate_initial_solution();
	populate_gamma();
	for (int t = 0; t < num_iterations; t++) {
		make_move(t);
		if (f(s) == 0) {
			break;
		}
	}
	for (int i = 0; i < NUM_VERTICES; i++) {
		cout << s[i] << " ";
	}
	cout << endl << "Number of conflicts: " << f(s) << endl;
	return 0;
}