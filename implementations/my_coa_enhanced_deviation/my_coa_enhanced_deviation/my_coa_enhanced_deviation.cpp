#define _SECURE_SCL 0
#define _USE_MATH_DEFINES
#define N_MAX 50

// Time taken (full parameters):
	// 100 vertices: 21.1037 secs (17 colours)
	// 500 vertices: 2.015 mins (80 colours)
		// With alpha = 5, it took 6.9852 mins but coloured it in 73 colours
	// 1000 vertices:
		// With alpha = 5, estimated is 34.2772 minutes (this seems more like it)

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

const int num_vertices = 450;

struct Cuckoo {
	int cuckoo[num_vertices];
	int fitness;
};

int adj_matrix[num_vertices][num_vertices];

int adj_list[num_vertices][num_vertices];
int adj_list_length[num_vertices];

int n_pop = 5;

const int alpha = 2;
const int num_iterations = 3000;
chrono::time_point<chrono::steady_clock> start;
const auto duration = chrono::minutes{ 5 };
const float p = 0.1;
const int min_eggs = 5;
const int max_eggs = 20;

float T = 10000;
const float beta = 1.0005f;

Cuckoo cuckoos[N_MAX];

int dist_matrix[N_MAX][N_MAX];

int num_eggs[N_MAX];
Cuckoo eggs[N_MAX * max_eggs];

mt19937 seed;
uniform_real_distribution<float> uni(0, 1);
uniform_real_distribution<float> angle(0, M_PI / 6);
uniform_int_distribution<int> random_vertex(0, num_vertices - 1);
uniform_int_distribution<int> random_egg_num(min_eggs, max_eggs);

void make_graph(float edge_probability) {  // Populates adj_matrix with a random graph
	for (int i = 0; i < num_vertices; i++) {
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

int num_conflicts(int * x) {
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

int max_colour(int * x) {
	int max = 0;
	for (int i = 0; i < num_vertices; i++) {
		if (x[i] > max) {
			max = x[i];
		}
	}
	return max + 1;
}

int colour_class_size[num_vertices];
int f(int * x) {
	int n = max_colour(x);
	for (int i = 0; i < n; i++) {
		colour_class_size[i] = 0;
	}
	for (int i = 0; i < num_vertices; i++) {
		colour_class_size[x[i]]++;
	}
	int ret = 0;
	for (int i = 0; i < n; i++) {
		ret -= colour_class_size[i] * colour_class_size[i];
	}
	return ret;
}

int D[num_vertices][num_vertices];
int px[num_vertices][num_vertices];
int px_length[num_vertices];
int py[num_vertices][num_vertices];
int py_length[num_vertices];
int min_D[num_vertices];
int mapping[num_vertices];
int colours[num_vertices];
bool taken[num_vertices];

bool compare_classes(int c1, int c2) {
	return min_D[c1] < min_D[c2];
}

void create_mapping(int * x, int * y) {
	int k = max_colour(x) + 1;
	int k_2 = max_colour(y) + 1;
	if (k_2 > k) {
		k = k_2;
	}
	for (int i = 0; i < k; i++) {
		min_D[i] = -1;
		colours[i] = i;
		px_length[i] = 0;
		py_length[i] = 0;
	}
	for (int i = 0; i < num_vertices; i++) {
		px_length[x[i]]++;
		py_length[y[i]]++;
	}
	for (int i = 0; i < k; i++) {
		taken[i] = false;
		for (int j = 0; j < k; j++) {
			D[j][i] = px_length[i] + py_length[j];
		}
	}
	for (int i = 0; i < num_vertices; i++) {
		D[y[i]][x[i]] -= 2;
		if (min_D[x[i]] == -1 || D[y[i]][x[i]] < min_D[x[i]]) {
			min_D[x[i]] = D[y[i]][x[i]];
		}
	}
	sort(begin(colours), begin(colours) + k, compare_classes);
	for (int i = 0; i < k; i++) {
		int c = colours[i];
		int min = -1;
		for (int j = 0; j < k; j++) {
			if (!taken[j] && (min == -1 || D[c][j] < D[c][min])) {
				min = j;
			}
		}
		mapping[c] = min;
		taken[min] = true;
	}
}

int d(int * x, int * y) {
	create_mapping(x, y);
	int d = 0;
	for (int i = 0; i < num_vertices; i++) {
		if (mapping[x[i]] != y[i]) {
			d++;
		}
	}
	return d;
}

void impose(int * x, int * y) {
	create_mapping(x, y);
	for (int i = 0; i < num_vertices; i++) {
		x[i] = mapping[x[i]];
	}
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

int order[num_vertices];
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

int I[num_vertices];
void migrate(int * x, int * y) {  // Migrates x towards y
	float r = uni(seed);
	float phi = angle(seed);
	int distance = d(x, y) * sqrt(1 + r * r - 2 * r * cos(phi));
	// Populate I with all vertices on which x and y disagree
	int I_length = 0;
	for (int i = 0; i < num_vertices; i++) {
		if (x[i] != y[i]) {
			I[I_length] = i;
			I_length++;
		}
	}
	random_shuffle(begin(I), begin(I) + I_length);
	for (int i = 0; i < distance; i++) {
		int c = 0;
		while (!(valid(I[i], c, y) && c != y[I[i]] && c != x[I[i]])) {
			c++;
		}
		y[I[i]] = c;
	}
	impose(x, y);
	random_shuffle(begin(I), begin(I) + I_length);
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

int main() {
	cout << "COA_deviation\n";
	ofstream ofile;
	ofile.open(results_directory + "le450_5a_coa_deviation.txt");
	//make_graph(0.5);
	read_graph("le450_5a.col");
	// Populate order array for generating cuckoos
	for (int i = 0; i < num_vertices; i++) {
		order[i] = i;
	}
	// Generate initial cuckoo population
	for (int i = 0; i < n_pop; i++) {
		generate_cuckoo(cuckoos[i].cuckoo);
		cuckoos[i].fitness = f(cuckoos[i].cuckoo);
	}
	start = chrono::high_resolution_clock::now();
	//int t = 0;
	//while(chrono::duration_cast<chrono::minutes>(chrono::high_resolution_clock::now() - start) < duration) {
	for (int t = 0; t < num_iterations; t++) {
		//t++;
		// Lay eggs
		int tot_eggs = 0;
		int egg = 0;
		// Generate number of eggs for each cuckoo
		for (int i = 0; i < n_pop; i++) {
			num_eggs[i] = random_egg_num(seed);
			tot_eggs += num_eggs[i];
		}
		for (int i = 0; i < n_pop; i++) {
			float elr = alpha * (num_eggs[i] / (float)tot_eggs) * num_vertices;  // Find cuckoo i's elr
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
				int d = eggs[i].fitness - cuckoos[i - N_MAX + n_pop].fitness;
				if (d < 0 || uni(seed) <= exp(-d / T)) {
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
		}
		else {  // If we can append all eggs without exceeding N_MAX then we can simply append all the eggs without worry
			for (int i = 0; i < tot_eggs; i++) {
				copy(begin(eggs[i].cuckoo), end(eggs[i].cuckoo), begin(cuckoos[n_pop].cuckoo));
				cuckoos[n_pop].fitness = eggs[i].fitness;
				n_pop++;
			}
		}
		T /= beta;
		// Cluster cuckoos to find goal point
		int gp = goal_point();
		int mp[num_vertices];
		// Migrate all cuckoos towards goal point
		for (int i = 0; i < n_pop; i++) {
			copy(begin(cuckoos[gp].cuckoo), end(cuckoos[gp].cuckoo), begin(mp));
			migrate(cuckoos[i].cuckoo, mp);
			cuckoos[i].fitness = f(cuckoos[i].cuckoo);
		}
		//cout << chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start).count() << endl;
		if (t % 10 == 0) {
			ofile << num_colours(cuckoos[0].cuckoo) << endl;
		}
	}
	ofile.close();
	sort(begin(cuckoos), end(cuckoos), compare_cuckoos);
	for (int i = 0; i < num_vertices; i++) {
		cout << cuckoos[0].cuckoo[i] << " ";
	}
	cout << endl << "Number of colours: " << num_colours(cuckoos[0].cuckoo) << endl;
	cout << "Number of conflicts: " << num_conflicts(cuckoos[0].cuckoo) << endl;
	//cout << "Number of iterations: " << t << endl;
	cout << "Time taken (ms): " << chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - start).count() << endl;
	return 0;
}