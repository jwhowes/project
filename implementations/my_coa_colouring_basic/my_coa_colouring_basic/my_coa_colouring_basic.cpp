#define _SECURE_SCL 0
#define NUM_VERTICES 10
#define N_MAX 50

// Time taken (full parameters):
	// 100 vertices: 24.379093 secs
	// 500 vertices: 4.5324301 mins (why did it say 68 secs before?)
	// 1000 vertices (estimated): 27.8457 mins

// TODO:
	// Fix whatever's causing (0, 0, ..., 0) to be returned (it's probably a memory leak)
		// It's caused by migration (at least I've narrowed it down (thank god))
	// Try out using iterators more
	// Try adj matrix
	// Try implementing edge lists with arrays
		// edge_lists[NUM_VERTICES][NUM_VERTICES]
		// edge_list_length[NUM_VERTICES]
	// Consider laying eggs in place
	// Look at OG github code to see what's different

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

/*vector<int> edge_list[NUM_VERTICES] = {
	{1, 4, 5},
	{0, 2, 6},
	{1, 3, 7},
	{2, 4, 8},
	{3, 0, 9},
	{0, 7, 8},
	{1, 8, 9},
	{2, 5, 9},
	{3, 5, 6},
	{4, 6, 7}
};*/

int adj_matrix[NUM_VERTICES][NUM_VERTICES];

int n_pop = 5;

const int alpha = 1;
const int num_iterations = 3000;
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
			ret += dist_matrix[S1[i]][S2[j]];
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
const int max_clustering_iterations = 10;
int S[N_MAX];
vector<int> clusters[k];
int db_sum_self[k];
uniform_int_distribution<int> random_cluster(0, k - 1);

int goal_point() {
	populate_dist_matrix();
	// Assign each cuckoo to a random cluster
	for (int i = 0; i < k; i++) {
		clusters[i].clear();
	}
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
	int t = 0;
	while (!converged && t < max_clustering_iterations) {
		t++;
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
				// Remove i from clusters[S[i]] and add it to clusters[cluster]
				for (auto it = clusters[S[i]].begin(); it != clusters[S[i]].end(); ++it) {
					if (*it == i) {
						clusters[S[i]].erase(it);
						break;
					}
				}
				clusters[cluster].push_back(i);
				converged = false;
				for (int j = 0; j < clusters[S[i]].size(); j++) {
					db_sum_self[S[i]] -= 2 * dist_matrix[i][clusters[S[i]][j]];
				}
				for (int j = 0; j < clusters[cluster].size(); j++) {
					db_sum_self[cluster] += 2 * dist_matrix[i][clusters[cluster][j]];
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
	for (int i = 1; i < clusters[best_cluster].size(); i++) {
		if (cuckoos[clusters[best_cluster][i]].fitness < cuckoos[gp].fitness) {
			gp = clusters[best_cluster][i];
		}
	}
	return gp;
}

bool valid(int v, int c, int * col) {
	for (int i = 0; i < NUM_VERTICES; i++) {
		if (adj_matrix[v][i] == 1 && col[i] == c) {
			return false;
		}
	}
	/*for (int i = 0; i < edge_list[v].size(); i++) {
		if (col[edge_list[v][i]] == c) {
			return false;
		}
	}*/
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

int I[NUM_VERTICES];
void migrate(int * x, int * y) {
	float r = uni(seed);
	// Populate a vector I with all vertices on which x and y disagree
	int I_length = 0;
	for (int i = 0; i < NUM_VERTICES; i++) {
		if (x[i] != y[i]) {
			I[I_length] = i;
			I_length++;
		}
	}
	for (int i = 0; i < r * I_length; i++) {
		int v = I[i];
		// Assign x y's colour for *it
		x[v] = y[v];
		// Clean up x
		for (int j = i + 1; j < I_length; j++) {
			// For every conflicting edge (v, *jit), replace *jit's colour with the smallest legal colour different than y[*jit]
			if (x[I[j]] == x[v] && adj_matrix[I[j]][v] == 1) {//find(edge_list[v].begin(), edge_list[v].end(), j) != edge_list[v].end()) {
				int c = 0;
				while (true) {
					if (c != y[I[j]] && valid(I[j], c, x)) {
						x[I[j]] = c;
						break;
					}
					c++;
				}
			}
		}
	}
}

int main() {
	make_graph(0.5);
	// Populate order array for generating cuckoos
	for (int i = 0; i < NUM_VERTICES; i++) {
		order[i] = i;
	}
	// Populate cuckoo array
	for (int i = 0; i < n_pop; i++) {
		generate_cuckoo(cuckoos[i].cuckoo);
		cuckoos[i].fitness = f(cuckoos[i].cuckoo);
	}
	auto start = chrono::high_resolution_clock::now();
	for (int t = 0; t < num_iterations; t++) {
		// Lay eggs
		int tot_eggs = 0;
		int egg = 0;
		for (int i = 0; i < n_pop; i++) {
			num_eggs[i] = random_egg_num(seed);
			tot_eggs += num_eggs[i];
		}
		for (int i = 0; i < n_pop; i++) {
			float elr = alpha * (num_eggs[i] / (float)tot_eggs) * NUM_VERTICES;
			for (int j = 0; j < num_eggs[i]; j++) {
				copy(begin(cuckoos[i].cuckoo), end(cuckoos[i].cuckoo), begin(eggs[egg].cuckoo));
				get_egg(eggs[egg].cuckoo, elr);
				eggs[egg].fitness = f(eggs[egg].cuckoo);
				egg++;
			}
		}
		// Kill a fraction p of the worst eggs
		tot_eggs *= (1 - p);
		if (n_pop + tot_eggs > N_MAX) {  // See if we can get away without using copy here
			for (int i = 0; i < N_MAX - n_pop; i++) {
				copy(begin(eggs[i].cuckoo), end(eggs[i].cuckoo), begin(cuckoos[n_pop + i].cuckoo));
				cuckoos[n_pop + i].fitness = eggs[i].fitness;
			}
			sort(begin(cuckoos), begin(cuckoos) + n_pop, reverse_compare_cuckoos);
			for (int i = N_MAX - n_pop; i < tot_eggs; i++) {
				if (eggs[i].fitness < cuckoos[i - N_MAX + n_pop].fitness) {
					copy(begin(eggs[i].cuckoo), end(eggs[i].cuckoo), begin(cuckoos[i - N_MAX + n_pop].cuckoo));
					cuckoos[i - N_MAX + n_pop].fitness = eggs[i].fitness;
				}
			}
			n_pop = N_MAX;
		}
		else {
			for (int i = 0; i < tot_eggs; i++) {
				copy(begin(eggs[i].cuckoo), end(eggs[i].cuckoo), begin(cuckoos[n_pop].cuckoo));
				cuckoos[n_pop].fitness = eggs[i].fitness;
				n_pop++;
			}
		}
		int gp = goal_point();
		for (int i = 0; i < n_pop; i++) {
			migrate(cuckoos[i].cuckoo, cuckoos[gp].cuckoo);
			cuckoos[i].fitness = f(cuckoos[i].cuckoo);
		}
		//cout << chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start).count() << endl;
	}
	for (int i = 0; i < NUM_VERTICES; i++) {
		cout << cuckoos[0].cuckoo[i] << " ";
	}
	cout << endl << "Number of colours: " << cuckoos[0].fitness << endl;
	cout << "Time taken: " << chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start).count();
	return 0;
}