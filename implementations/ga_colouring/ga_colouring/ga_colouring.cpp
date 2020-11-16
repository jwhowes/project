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

using namespace std;
using namespace boost::random;

const int k = 3;
const int NUM_VERTICES = 10;

int adj_matrix[NUM_VERTICES][NUM_VERTICES] = {
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

const int pop_size = 50;
const int num_iterations = 100;
const float mutation_prob = 0.1;

vector<int> population[NUM_VERTICES][k];

mt19937 seed;
uniform_real_distribution<float> uni(0, 1);
uniform_int_distribution<int> random_colour(0, k - 1);
uniform_int_distribution<int> random_vertex(0, NUM_VERTICES - 1);

int f(vector<int> * x) {
	int num = 0;
	for(int c = 0; c < k; c++){
		for (int i = 0; i < x[c].size(); i++) {
			for (int j = 0; j < i; j++) {
				if (adj_matrix[x[c][i]][x[c][j]] == 1) {
					num++;
				}
			}
		}
	}
	return num;
}

vector<int> parents[2][k];
void gpx(vector<int> * p1, vector<int> * p2, vector<int> * x) {
	vector<int> uncoloured(NUM_VERTICES);
	int parent;
	for (int i = 0; i < NUM_VERTICES; i++) {
		uncoloured[i] = i;
	}
	for (int i = 0; i < k; i++) {
		parents[0][i].assign(p1[i].begin(), p1[i].end());
		parents[1][i].assign(p2[i].begin(), p2[i].end());
		x[i].clear();
	}
	for (int i = 0; i < k; i++) {
		parent = i % 2;
		int m = 0;
		for (int j = 1; j < k; j++) {
			if (parents[parent][j].size() > parents[parent][m].size()) {
				m = j;
			}
		}
		x[i].assign(parents[parent][m].begin(), parents[parent][m].end());
		for (int i = 0; i < parents[parent][m].size(); i++) {
			uncoloured.erase(remove(uncoloured.begin(), uncoloured.end(), parents[parent][m][i]), uncoloured.end());
		}
		parents[parent][m].clear();
		for (int j = 0; j < x[i].size(); j++) {
			for (int l = 0; l < k; l++) {
				auto pos = find(parents[1-parent][l].begin(), parents[1-parent][l].end(), x[i][j]);
				if (pos != parents[1-parent][l].end()) {
					parents[1-parent][l].erase(pos);
					break;
				}
			}
		}
	}
	for (int i = 0; i < uncoloured.size(); i++) {
		x[random_colour(seed)].push_back(uncoloured[i]);
	}
}

void generate_member(vector<int> * x) {
	for (int i = 0; i < NUM_VERTICES; i++) {
		x[random_colour(seed)].push_back(i);
	}
}

int main(){
	vector<int> p1[k];
	vector<int> p2[k];
	vector<int> x[k];
	generate_member(p1);
	generate_member(p2);
	for (int c = 0; c < k; c++) {
		for (int v = 0; v < p1[c].size(); v++) {
			cout << p1[c][v] << " ";
		}
		cout << endl;
	}
	cout << endl;
	for (int c = 0; c < k; c++) {
		for (int v = 0; v < p2[c].size(); v++) {
			cout << p2[c][v] << " ";
		}
		cout << endl;
	}
	cout << endl;
	gpx(p1, p2, x);
	for (int c = 0; c < k; c++) {
		for (int v = 0; v < x[c].size(); v++) {
			cout << x[c][v] << " ";
		}
		cout << endl;
	}
	cout << endl;
	return 0;
}