#define _USE_MATH_DEFINES
#define _SECURE_SCL 0
#define NUM_VERTICES 100
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

// Running times (eggs 5 to 20, pop 5 to 50, iterations 3000)
	// 10 vertices: 17.94 secs
	

using namespace std;
using namespace boost::random;

typedef struct Cuckoo {
	int cuckoo[NUM_VERTICES];
	int fitness;
};

vector<int> edge_list[NUM_VERTICES];/* = {
	{1, 4, 5},
	{0, 2, 6},
	{1, 3, 7},
	{2, 4, 8},
	{3, 1, 9},
	{0, 7, 8},
	{1, 8, 9},
	{2, 5, 8},
	{3, 5, 6},
	{4, 6, 7}
};*/

int n_pop = 5;

const float alpha = 1;
const int num_iterations = 3000;
const float p = 0.1;
const int min_eggs = 1;
const int max_eggs = 1;

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
				edge_list[i].push_back(j);
				edge_list[j].push_back(i);
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

int d(int x, int y) {
	int ret = 0;
	for (int i = 0; i < NUM_VERTICES; i++) {
		if (cuckoos[x].cuckoo[i] != cuckoos[y].cuckoo[i]) {
			ret++;
		}
	}
	return ret;
}

int d_bar_sum(vector<int> S1, vector<int> S2) {
	int ret = 0;
	for (auto it = S1.begin(); it != S1.end(); ++it) {
		for (auto jit = S2.begin(); jit != S2.end(); ++jit) {
			ret += dist_matrix[*it][*jit];
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
			dist_matrix[i][j] = d(i, j);
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
		for (auto it = clusters[i].begin(); it != clusters[i].end(); ++it) {
			for (auto jit = clusters[i].begin(); jit < it; ++jit) {
				db_sum_self[i] += 2 * dist_matrix[*it][*jit];
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
				for (auto it = clusters[S[i]].begin(); it != clusters[S[i]].end(); ++it){
					db_sum_self[S[i]] -= 2 * dist_matrix[i][*it];
				}
				for (auto it = clusters[cluster].begin(); it != clusters[cluster].end(); ++it) {
					db_sum_self[cluster] += 2 * dist_matrix[i][*it];
				}
				S[i] = cluster;
			}
		}
	}
	int best_cluster = -1;
	float best_cluster_mean_fitness = 0;
	for (int i = 0; i < k; i++) {
		float mean_cluster_fitness = 0;
		for (auto it = clusters[i].begin(); it != clusters[i].end(); ++it){
			mean_cluster_fitness += cuckoos[*it].fitness;
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

bool valid(int v, int c, int col) {
	for (auto it = edge_list[v].begin(); it != edge_list[v].end(); ++it){
		if (cuckoos[col].cuckoo[*it] == c) {
			return false;
		}
	}
	return true;
}

int order[NUM_VERTICES];

void generate_cuckoo(int cuckoo) {
	// Choose a random ordering of vertices
	random_shuffle(begin(order), end(order));
	for (int v : order) {
		int c = 0;
		while (true) {
			if (valid(v, c, cuckoo)) {
				cuckoos[cuckoo].cuckoo[v] = c;
				break;
			}
			c++;
		}
	}
}

void get_egg(int cuckoo, float elr) {
	int num = uniform_int_distribution<int>(0, elr)(seed);
	for (int i = 0; i < num; i++) {
		int v = random_vertex(seed);
		int c = 0;
		while (true) {
			if (c != cuckoos[cuckoo].cuckoo[v] && valid(v, c, cuckoo)) {
				cuckoos[cuckoo].cuckoo[v] = c;
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
void migrate(int x, int y) {
	float r = uni(seed);
	// Populate a vector I with all vertices on which x and y disagree
	int I_length = 0;
	for (int i = 0; i < NUM_VERTICES; i++) {
		if (cuckoos[x].cuckoo[i] != cuckoos[y].cuckoo[i]) {
			I[I_length] = i;
			I_length++;
		}
	}
	for (int i = 0; i < r * I_length; i++){
		int v = I[i];
		// Assign x y's colour for *it
		cuckoos[x].cuckoo[v] = cuckoos[y].cuckoo[v];
		// Clean up x
		for (int j = i + 1; j < NUM_VERTICES; j++) {
			// For every conflicting edge (v, *jit), replace *jit's colour with the smallest legal colour different than y[*jit]
			if (cuckoos[x].cuckoo[j] == cuckoos[x].cuckoo[v] && find(edge_list[v].begin(), edge_list[v].end(), j) != edge_list[v].end()) {
				int c = 0;
				while (true) {
					if (c != cuckoos[y].cuckoo[v] && valid(j, c, x)) {
						cuckoos[x].cuckoo[j] = c;
						break;
					}
					c++;
				}
			}
		}
	}
}

int main(){  // I'm pretty sure the bottleneck is clustering
	make_graph(0.5);
	// Populate order array for generating cuckoos
	for (int i = 0; i < NUM_VERTICES; i++) {
		order[i] = i;
	}
	
	// Populate cuckoo array
	for (int i = 0; i < n_pop; i++) {
		generate_cuckoo(i);
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
				get_egg(egg, elr);
				eggs[egg].fitness = f(eggs[egg].cuckoo);
				egg++;
			}
		}
		// Kill a fraction p of the worst eggs
		sort(begin(eggs), begin(eggs) + tot_eggs, compare_cuckoos);
		tot_eggs *= (1 - p);
		if (n_pop + tot_eggs > N_MAX) {  // See if we can get away without using copy here
			// Add eggs until n_pop = N_MAX
			for (int i = 0; i < N_MAX - n_pop; i++) {
				copy(begin(eggs[i].cuckoo), end(eggs[i].cuckoo), begin(cuckoos[n_pop + i].cuckoo));
				cuckoos[n_pop + i].fitness = eggs[i].fitness;
			}
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
				copy(begin(eggs[i].cuckoo), end(eggs[i].cuckoo), begin(cuckoos[n_pop].cuckoo));
				cuckoos[n_pop].fitness = eggs[i].fitness;
			}
		}
		int gp = goal_point();
		//int gp = 0;
		for (int i = 0; i < n_pop; i++) {
			migrate(i, gp);
			cuckoos[i].fitness = f(cuckoos[i].cuckoo);
		}
		//cout << chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start).count() << endl;
	}
	for (int i = 0; i < NUM_VERTICES; i++) {
		cout << cuckoos[0].cuckoo[i] << " ";
	}
	cout << endl << "Time taken: " << chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start).count();
	return 0;
}