#define _SECURE_SCL 0
#define NUM_VERTICES 250
#define N_MAX 50

// DON'T FORGET TO ADD A MINIMISE FUNCTION (I think this one needs it)

// In this version, if we can add all eggs without exceeding N_MAX, then we do it.
	// Else an egg replaces the parent (or a previous egg) if it is better or if it passes the SA check
// Basically, it's pretty rubbish atm but I think it could potentially be improved
	// Maybe rather than replacing the parent, replace the worst cuckoo (then move on to replacing the next cuckoo when that fails)
	// This means we'll probably have to have two separate loops (as we did before). One for adding up to N_MAX then another for replacing

// Results:
	// 100 vertices: 24 colours, 17.4263 seconds

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
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

const string graph_directory = "C:/Users/taydo/OneDrive/Documents/computer_science/year3/project/implementations/graphs/";

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
const int num_iterations = 10;
chrono::time_point<chrono::steady_clock> start;
const auto duration = chrono::minutes{5};
const float p = 0.1;
const int min_eggs = 5;
const int max_eggs = 20;

float T = 10000;
const float beta = 1.0005f;

Cuckoo cuckoos[N_MAX];

int dist_matrix[N_MAX][N_MAX];

int num_eggs[N_MAX];

mt19937 seed;
uniform_real_distribution<float> uni(0, 1);
uniform_int_distribution<int> random_vertex(0, NUM_VERTICES - 1);
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

int num_colours(int * x) {  // Returns the fitness of a cuckoo (number of colours used)
	int max = 0;
	for (int i = 0; i < NUM_VERTICES; i++) {
		if (x[i] > max) {
			max = x[i];
		}
	}
	return max + 1;
}

int * classes;
int f(int * x) {
	int num_classes = num_colours(x);
	delete[] classes;
	classes = new int[num_classes];
	int ret = 0;
	for (int i = 0; i < NUM_VERTICES; i++) {
		classes[x[i]]++;
	}
	for (int i = 0; i < num_classes; i++) {
		ret -= classes[i] * classes[i];
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
				break;
			}
			c++;
		}
	}
}

void get_egg(int * cuckoo, float elr) {  // Populates cuckoo with a random colouring within distance elr of it
	int num = uniform_int_distribution<int>(0, elr)(seed);
	for (int i = 0; i < num; i++) {
		int v = random_vertex(seed);
		int c = 0;
		while (true) {
			if (c != cuckoo[v] && valid(v, c, cuckoo)) {  // Changes num vertices to the smallest valid colour (different to their current one)
				cuckoo[v] = c;
				break;
			}
			c++;
		}
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
void migrate(int * x, int * y) {  // Migrates x towards y
	float r = uni(seed);
	// Populate I with all vertices on which x and y disagree
	int I_length = 0;
	for (int i = 0; i < NUM_VERTICES; i++) {
		if (x[i] != y[i]) {
			I[I_length] = i;
			I_length++;
		}
	}
	for (int i = 0; i < r * I_length; i++) {  // For a random
		int v = I[i];
		// Assign x y's colour for *it
		x[v] = y[v];
		// Clean up x
		for (int j = i + 1; j < I_length; j++) {
			// For every vertex I[j] in I with I-index j > i that conflicts with v, assign the smallest legal colour different thatn y[I[j]]
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

int egg_temp[NUM_VERTICES];
int parent_temp[NUM_VERTICES];
int main() {
	cout << "COA Full SA\n";
	//make_graph(0.5);
	read_graph("r250.5.col");
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
	int t = 0;
	while(chrono::duration_cast<chrono::minutes>(chrono::high_resolution_clock::now() - start) < duration) {
	//for(int t = 0; t < num_iterations; t++){
		t++;
		// Lay eggs
		int tot_eggs = 0;
		int egg = 0;
		// Generate number of eggs for each cuckoo
		for (int i = 0; i < n_pop; i++) {
			num_eggs[i] = random_egg_num(seed);
			tot_eggs += num_eggs[i];
		}
		tot_eggs *= (1 - p);
		int egg_fitness;
		egg = 0;
		for (int i = 0; i < n_pop; i++) {
			copy(begin(cuckoos[i].cuckoo), end(cuckoos[i].cuckoo), begin(parent_temp));
			float elr = alpha * (num_eggs[i] / (float)tot_eggs) * NUM_VERTICES;
			for (int j = 0; j < num_eggs[i]; j++) {
				copy(begin(parent_temp), end(parent_temp), begin(egg_temp));
				get_egg(egg_temp, elr);
				if (n_pop + egg < N_MAX) {
					copy(begin(egg_temp), end(egg_temp), begin(cuckoos[n_pop + egg].cuckoo));
					cuckoos[n_pop + egg].fitness = f(cuckoos[n_pop + egg].cuckoo);
					egg++;
				} else {
					egg_fitness = f(egg_temp);
					int d = egg_fitness - cuckoos[i].fitness;
					if (d < 0 || uni(seed) <= exp(-d / T)) {
						copy(begin(egg_temp), end(egg_temp), begin(cuckoos[i].cuckoo));
						cuckoos[i].fitness = egg_fitness;
					}
				}
			}
		}
		if (n_pop + tot_eggs < N_MAX) {
			n_pop += tot_eggs;
		}
		else {
			n_pop = N_MAX;
		}
		T /= beta;
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
	cout << "Time taken (ms): " << chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - start).count() << endl;
	cout << "Number of iterations: " << t << endl;
	return 0;
}