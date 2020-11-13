#define _USE_MATH_DEFINES
#define _SECURE_SCL 0
#define NUM_VERTICES 10
#define N_MAX 50

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

typedef struct Cuckoo {
	int cuckoo[NUM_VERTICES];
	int fitness;
};

vector<vector<int>> edge_list = {
	{1, 4, 5},
	{0, 2, 7},
	{1, 3, 7},
	{2, 4, 8},
	{3, 1, 9},
	{0, 7, 8},
	{1, 8, 9},
	{2, 5, 8},
	{3, 5, 6},
	{4, 6, 7}
};

int n_pop = 5;

const int alpha = 10;
const int num_iterations = 1;
const float p = 0.1;
const int min_eggs = 5;
const int max_eggs = 20;

Cuckoo cuckoos[N_MAX];

int dist_matrix[N_MAX][N_MAX];

int num_eggs[N_MAX];
Cuckoo eggs[N_MAX * max_eggs];

mt19937 seed;
uniform_real_distribution<float> uni(0, 1);
uniform_int_distribution<int> random_vertex(0, NUM_VERTICES - 1);
uniform_int_distribution<int> random_egg_num(min_eggs, max_eggs);

int f(int * x) {
	int max = 0;
	for (int i = 0; i < NUM_VERTICES; i++) {
		if (x[i] > max) {
			max = x[i];
		}
	}
	return max + 1;
}

int d(int * x, int * y) {
	int ret = 0;
	for (int i = 0; i < NUM_VERTICES; i++) {
		if (x[i] != y[i]) {
			ret++;
		}
	}
	return ret;
}

int d_bar_sum(vector<int> S1, vector<int> S2) {
	int ret = 0;
	for (int i = 0; i < S1.size(); i++) {
		for (int j = 0; j < S2.size(); j++) {
			ret += dist_matrix[i][j];
		}
	}
	return ret;
}

float tri_dist(int i, int j, int * dbss, vector<int> * clusters) {
	if (clusters[j].size() == 0) {
		return -1;
	}
	return 2 * d_bar_sum({ i }, clusters[j]) / clusters[j].size() - dbss[j] / (clusters[j].size() * clusters[j].size());
}

void populate_dist_matrix() {
	for (int i = 0; i < n_pop; i++) {
		for (int j = 0; j < i; j++) {
			dist_matrix[i][j] = d(cuckoos[i].cuckoo, cuckoos[j].cuckoo);
			dist_matrix[j][i] = dist_matrix[i][j];
		}
	}
}

const int k = 3;
int S[N_MAX];
vector<int> clusters[k];
int db_sum_self[k];
uniform_int_distribution<int> random_cluster(0, k - 1);

int goal_point() {
	populate_dist_matrix();
	// Assign each cuckoo to a random cluster
	for (int i = 0; i < n_pop; i++) {
		S[i] = random_cluster(seed);
		clusters[S[i]].push_back(i);
	}
	// Initialise db_sum_self
	for (int i = 0; i < k; i++) {
		db_sum_self[i] = 0;
		for (int j = 0; j < clusters[i].size(); j++) {
			for (int l = 0; l < j; l++) {
				db_sum_self[i] += 2 * dist_matrix[clusters[i][j]][clusters[i][l]];
			}
		}
	}
	bool converged = false;
	while (!converged) {
		converged = true;
		for (int i = 0; i < n_pop; i++) {
			int cluster = -1;
			float tdc = -1;
			for (int j = 0; j < k; j++) {
				float td = tri_dist(i, j, db_sum_self, clusters);
				if ((cluster == -1 || td < tdc) && td != -1) {
					cluster = j;
					tdc = td;
				}
			}
			if (cluster != S[i]) {
				// Remove i from clusters[S[i]]
				for (auto it = clusters[S[i]].begin(); it != clusters[S[i]].end(); ++it) {
					if (*it == i) {
						clusters[S[i]].erase(it);
						break;
					}
				}
				clusters[cluster].push_back(i);
				converged = false;
				for (int j = 0; j < clusters[S[i]].size(); j++) {
					db_sum_self[S[i]] -= 2 * dist_matrix[i][j];
				}
				for (int j = 0; j < clusters[cluster].size(); j++) {
					db_sum_self[cluster] += 2 * dist_matrix[i][j];
				}
				S[i] = cluster;
			}
		}
	}
	int best_cluster = -1;
	float best_cluster_mean_fitness = 0;
	for (int i = 0; i < k; i++) {
		float mean_cluster_fitness = 0;
		for (int j = 0; j < clusters[i].size(); j++) {
			mean_cluster_fitness += cuckoos[clusters[i][j]].fitness;
		}
		mean_cluster_fitness /= clusters[i].size();
		if (clusters[i].size() > 0 && (best_cluster == -1 || mean_cluster_fitness < best_cluster_mean_fitness)) {
			best_cluster = i;
			best_cluster_mean_fitness = mean_cluster_fitness;
		}
	}
	int gp = clusters[best_cluster][0];
	for (auto it = clusters[best_cluster].begin(); it != clusters[best_cluster].end(); ++it) {
		if (cuckoos[*it].fitness < cuckoos[gp].fitness) {
			gp = *it;
		}
	}
	return gp;
}

bool valid(int v, int c, int * col) {
	for (int i = 0; i < edge_list[v].size(); i++) {
		if (col[i] == c) {
			return false;
		}
	}
	return true;
}

int order[NUM_VERTICES];

void generate_cuckoo(int * cuckoo) {
	// Choose a random ordering of vertices
	random_shuffle(begin(order), end(order));
	for (int v : order) {
		int c = 0;
		while (true) {
			if (valid(v, c, cuckoo)) {
				cuckoo[v] = c;
				break;
			}
			c++;
		}
	}
}

void get_egg(int * cuckoo, float elr) {
	int num = uniform_int_distribution<int>(0, elr)(seed);
	for (int i = 0; i < num; i++) {
		int v = random_vertex(seed);
		int c = 0;
		while (true) {
			if (c != cuckoo[v] && valid(v, c, cuckoo)) {
				cuckoo[v] = c;
				break;
			}
			c++;
		}
	}
}

bool compare_cuckoos(Cuckoo c1, Cuckoo c2) {
	return c1.fitness < c2.fitness;
}

bool reverse_compare_cuckoos(Cuckoo c1, Cuckoo c2) {
	return c1.fitness > c2.fitness;
}

int main(){
	// Populate order array for generating cuckoos
	for (int i = 0; i < NUM_VERTICES; i++) {
		order[i] = i;
	}
	// Populate cuckoo array
	for (int i = 0; i < n_pop; i++) {
		generate_cuckoo(cuckoos[i].cuckoo);
		cuckoos[i].fitness = f(cuckoos[i].cuckoo);
	}
	for (int t = 0; t < num_iterations; t++) {
		// Lay eggs
		int tot_eggs = 0;
		int egg = 0;
		for (int i = 0; i < n_pop; i++) {
			num_eggs[i] = random_egg_num(seed);
			tot_eggs += num_eggs[i];
		}
		for (int i = 0; i < n_pop; i++) {
			float elr = alpha * num_eggs[i] / tot_eggs * NUM_VERTICES;
			for (int j = 0; j < num_eggs[i]; j++) {
				copy(begin(cuckoos[i].cuckoo), end(cuckoos[i].cuckoo), begin(eggs[egg].cuckoo));
				get_egg(eggs[egg].cuckoo, elr);
				eggs[egg].fitness = f(eggs[egg].cuckoo);
				egg++;
			}
		}
		// Kill a fraction p of the worst eggs
		sort(begin(eggs), begin(eggs) + tot_eggs, compare_cuckoos);
		tot_eggs *= (1 - p);
		if (n_pop + tot_eggs > N_MAX) {  // See if we can get away without using copy here
			// Add eggs until n_pop = N_MAX
			copy(begin(eggs), begin(eggs) + N_MAX - n_pop - 1, begin(cuckoos) + n_pop);
			// Let best remaining eggs replace worst cuckoos
			sort(begin(cuckoos), begin(cuckoos) + n_pop, reverse_compare_cuckoos);
			for (int i = N_MAX - n_pop; i < tot_eggs; i++) {
				if (eggs[i].fitness < cuckoos[i - N_MAX + n_pop].fitness) {
					copy(begin(eggs[i].cuckoo), end(eggs[i].cuckoo), begin(cuckoos[i - N_MAX + n_pop].cuckoo));
					cuckoos[i - N_MAX + n_pop].fitness = eggs[i].fitness;
				}
			}
			n_pop = N_MAX;
		}else{
			for (int i = 0; i < tot_eggs; i++) {
				n_pop++;
				cuckoos[n_pop] = eggs[i];
			}
		}
		int gp = goal_point();
		cout << gp;
		// Migrate cuckoos towards goal point (should be easy as the python code is alread inplace)
	}
	return 0;
}