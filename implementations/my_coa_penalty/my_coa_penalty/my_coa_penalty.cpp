#define _SECURE_SCL 0
#define NUM_VERTICES 1000
#define N_MAX 50

// Results (alpha = 1):
	// 100 vertices: 36.7039 secs (18 colours)
	// 500 vertices: 8.51443 mins (68 colours)
	// 1000 vertices: 31.2307 mins (123 colours)

// TODO:
	// See if doing a partition variant would improve run time (the fitness function would be much quicker)
		// Use the array based partition approach from ga_colouring

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

struct Cuckoo {
	int cuckoo[NUM_VERTICES];
	int fitness;
};

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

int adj_list[NUM_VERTICES][NUM_VERTICES];/* = {
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
int adj_list_length[NUM_VERTICES];// = { 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 };

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

int chromatic_bound = 0;

mt19937 seed;
uniform_real_distribution<float> uni(0, 1);
uniform_int_distribution<int> random_vertex(0, NUM_VERTICES - 1);
uniform_int_distribution<int> random_colour;
uniform_int_distribution<int> random_egg_num(min_eggs, max_eggs);

void make_graph(float edge_probability) {  // Populates adj_matrix with a random graph
	for (int i = 0; i < NUM_VERTICES; i++) {
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

int num_colours(int * x) {
	vector<int> colours_used;
	for (int i = 0; i < NUM_VERTICES; i++) {
		if (find(colours_used.begin(), colours_used.end(), x[i]) == colours_used.end()) {
			colours_used.push_back(x[i]);
		}
	}
	return colours_used.size();
}

int num_conflicts(int * x) {
	int ret = 0;
	for (int i = 0; i < NUM_VERTICES; i++) {
		for (int j = 0; j < i; j++) {
			if (adj_matrix[i][j] == 1 && x[i] == x[j]) {
				ret += 1;
			}
		}
	}
	return ret;
}

int max_colour(int * x) {
	int max = 0;
	for (int i = 0; i < NUM_VERTICES; i++) {
		if (x[i] > max) {
			max = x[i];
		}
	}
	return max + 1;
}

int * colour_class_size;
int * edge_conflicts;
int f(int * x) {
	delete[] colour_class_size;
	delete[] edge_conflicts;
	int n = max_colour(x);
	if (n < chromatic_bound) {
		chromatic_bound = n;
	}
	colour_class_size = new int[n];
	edge_conflicts = new int[n];
	for (int i = 0; i < n; i++) {
		colour_class_size[i] = 0;
		edge_conflicts[i] = 0;
	}
	for (int i = 0; i < NUM_VERTICES; i++) {
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

int d(int * x, int * y) {  // Returns the hamming distance between two colourings
	int ret = 0;
	for (int i = 0; i < NUM_VERTICES; i++) {
		if (x[i] != y[i]) {
			ret++;
		}
	}
	return ret;
}

int d_bar_sum(vector<int> S1, vector<int> & S2) {  // Calculates d_bar between two clusters S1 and S2
	// Doesn't divide by clusters size to make changing clusters more efficient
	int ret = 0;
	for (int i = 0; i < S1.size(); i++) {
		for (int j = 0; j < S2.size(); j++) {
			ret += dist_matrix[S1[i]][S2[j]];
		}
	}
	return ret;
}

float tri_dist(int i, int j, int * dbss, vector<int> * clusters) {  // Calculates tri dist between cuckoo i and cluster j
	if (clusters[j].size() == 0) {
		return -1;
	}
	return 2 * d_bar_sum({ i }, clusters[j]) / clusters[j].size() - dbss[j] / (clusters[j].size() * clusters[j].size());
}

void populate_dist_matrix() {  // Populates the distance matrix s.t. d[i][j] = d(cuckoos[i], cuckoos[j])
	for (int i = 0; i < n_pop; i++) {
		for (int j = 0; j < i; j++) {
			dist_matrix[i][j] = d(cuckoos[i].cuckoo, cuckoos[j].cuckoo);
			dist_matrix[j][i] = dist_matrix[i][j];
		}
	}
}

const int k = 3;  // The number of clusters
const int max_clustering_iterations = 10;
int S[N_MAX];
vector<int> clusters[k];
int db_sum_self[k];
uniform_int_distribution<int> random_cluster(0, k - 1);

int goal_point() {
	populate_dist_matrix();
	// Assign each cuckoo to a random cluster
	for (int i = 0; i < k; i++) {
		clusters[i].clear();  // Clear clusters from previous iteration
	}
	for (int i = 0; i < n_pop; i++) {
		S[i] = random_cluster(seed);  // S[i] stores the cluster of i
		clusters[S[i]].push_back(i);  // clusters[i] stores cluster i
	}
	// Initialise db_sum_self
	// db_sum_self stores the db sum of each cluster to itself, thus d_bar for cluster i is just db_sum_self[i]/clusters[i].size()
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
	while (!converged && t < max_clustering_iterations) {  // Not entirely sure if max_clustering_iterations is strictly necessary
		t++;
		converged = true;
		// Find the nearest cluster for each cuckoo
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
				// Update db_sum_self for affected clusters
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
	// Find the best cluster (the cluster with best mean fitness value)
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
	// Migration goal point is the best cuckoo of the best cluster
	int gp = clusters[best_cluster][0];
	for (int i = 1; i < clusters[best_cluster].size(); i++) {
		if (cuckoos[clusters[best_cluster][i]].fitness < cuckoos[gp].fitness) {
			gp = clusters[best_cluster][i];
		}
	}
	return gp;
}

bool valid(int v, int c, int * col) {  // Returns whether or not vertex v can be coloured colour c in colouring col (legally)
	for (int i = 0; i < adj_list_length[v]; i++) {
		if (col[adj_list[v][i]] == c) {  // If v is adjacent to some vertex i coloured c then this is not a valid assignment
			return false;
		}
	}
	return true;  // If no such i can be found then the assignment is valid
}

int order[NUM_VERTICES];
void generate_cuckoo(int * cuckoo) {  // Populates cuckoo with a random valid colouring
	// Choose a random ordering of vertices
	random_shuffle(begin(order), end(order));
	for (int v : order) {
		int c = 0;
		while (true) {
			if (valid(v, c, cuckoo)) {
				cuckoo[v] = c;
				if (c > chromatic_bound) {
					chromatic_bound = c;
				}
				break;
			}
			c++;
		}
	}
}

void get_egg(int * cuckoo, float elr) {
	int num = uniform_int_distribution<int>(0, (int)elr)(seed);
	random_colour = uniform_int_distribution<int>(0, chromatic_bound);
	for (int i = 0; i < num; i++) {
		cuckoo[random_vertex(seed)] = random_colour(seed);
	}
}

// Comparison operators for sorting
bool compare_cuckoos(Cuckoo & c1, Cuckoo & c2) {
	return c1.fitness < c2.fitness;
}
bool reverse_compare_cuckoos(Cuckoo & c1, Cuckoo & c2) {
	return c1.fitness > c2.fitness;
}

int I[NUM_VERTICES];
void migrate(int * x, int * y) {
	float r = uni(seed);
	int I_length = 0;
	for (int i = 0; i < NUM_VERTICES; i++) {
		if (x[i] != y[i]) {
			I[I_length] = i;
			I_length++;
		}
	}
	for (int i = 0; i < r * I_length; i++) {
		x[I[i]] = y[I[i]];
	}
}

int main() {
	make_graph(0.5);
	// Populate order array for generating cuckoos
	for (int i = 0; i < NUM_VERTICES; i++) {
		order[i] = i;
	}
	// Generate initial cuckoo population
	for (int i = 0; i < n_pop; i++) {
		generate_cuckoo(cuckoos[i].cuckoo);
		cuckoos[i].fitness = f(cuckoos[i].cuckoo);
	}
	auto start = chrono::high_resolution_clock::now();
	for (int t = 0; t < num_iterations; t++) {
		// Lay eggs
		int tot_eggs = 0;
		int egg = 0;
		// Generate number of eggs for each cuckoo
		for (int i = 0; i < n_pop; i++) {
			num_eggs[i] = random_egg_num(seed);
			tot_eggs += num_eggs[i];
		}
		for (int i = 0; i < n_pop; i++) {
			float elr = alpha * (num_eggs[i] / (float)tot_eggs) * NUM_VERTICES;  // Find cuckoo i's elr
			for (int j = 0; j < num_eggs[i]; j++) {  // Generate num_eggs[i] eggs for cuckoo i and append them to the eggs array
				copy(begin(cuckoos[i].cuckoo), end(cuckoos[i].cuckoo), begin(eggs[egg].cuckoo));
				get_egg(eggs[egg].cuckoo, elr);
				eggs[egg].fitness = f(eggs[egg].cuckoo);
				egg++;
			}
		}
		// Kill a fraction p of the worst eggs
		sort(begin(eggs), begin(eggs) + tot_eggs, compare_cuckoos);
		tot_eggs *= (1 - p);
		if (n_pop + tot_eggs > N_MAX) {  // If appending the eggs to cuckoos will exceed N_MAX then we append the eggs and kill the worst cuckoos until pop < N_MAX
			// Here is a way of achieving the above without exceeding the bounds of the cuckoos array
			for (int i = 0; i < N_MAX - n_pop; i++) {  // First append as many eggs as we can
				copy(begin(eggs[i].cuckoo), end(eggs[i].cuckoo), begin(cuckoos[n_pop + i].cuckoo));
				cuckoos[n_pop + i].fitness = eggs[i].fitness;
			}
			sort(begin(cuckoos), begin(cuckoos) + n_pop, reverse_compare_cuckoos);
			for (int i = N_MAX - n_pop; i < tot_eggs; i++) {  // Then iteratively replace the worst cuckoo with the best egg (if the egg is indeed better)
				if (eggs[i].fitness < cuckoos[i - N_MAX + n_pop].fitness) {
					copy(begin(eggs[i].cuckoo), end(eggs[i].cuckoo), begin(cuckoos[i - N_MAX + n_pop].cuckoo));
					cuckoos[i - N_MAX + n_pop].fitness = eggs[i].fitness;
				}
				else {
					// Eggs are increasing in fitness value while cuckoos are decreasing
					// Hence, if we reach a point where a cuckoo has better fitness than an egg, no eggs will ever have a better fitness than any cuckoo from that point so we can exit the loop
					break;

				}
			}
			n_pop = N_MAX;
		} else {  // If we can append all eggs without exceeding N_MAX then we can simply append all the eggs without worry
			for (int i = 0; i < tot_eggs; i++) {
				copy(begin(eggs[i].cuckoo), end(eggs[i].cuckoo), begin(cuckoos[n_pop].cuckoo));
				cuckoos[n_pop].fitness = eggs[i].fitness;
				n_pop++;
			}
		}
		// Cluster cuckoos to find goal point
		int gp = goal_point();
		// Migrate all cuckoos towards goal point
		for (int i = 0; i < n_pop; i++) {
			migrate(cuckoos[i].cuckoo, cuckoos[gp].cuckoo);
			cuckoos[i].fitness = f(cuckoos[i].cuckoo);
		}
		//cout << chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start).count() << endl;
	}
	sort(begin(cuckoos), end(cuckoos), compare_cuckoos);
	for (int i = 0; i < NUM_VERTICES; i++) {
		cout << cuckoos[0].cuckoo[i] << " ";
	}
	cout << endl << "Number of colours: " << num_colours(cuckoos[0].cuckoo) << endl;
	cout << "Number of conflicts: " << num_conflicts(cuckoos[0].cuckoo) << endl;
	cout << "Time taken (seconds): " << chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start).count() / (float)1000000 << endl;
	return 0;
}