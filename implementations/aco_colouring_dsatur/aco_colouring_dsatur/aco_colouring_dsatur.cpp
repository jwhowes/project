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

// They use small ant populations (< 50) on very small graphs (|V| <= 100) so I don't think this is a particularly efficient algorithm

// TODO:
	// Figure out why it's returning invalid colourings
		// weight_sum is sometimes 0 (if I use tau + 1 in the transition rule then this doesn't happen)

using namespace std;
using namespace boost::random;

const int num_vertices = 10;

int adj_matrix[num_vertices][num_vertices] = {
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

float tau[num_vertices][num_vertices];

const int num_ants = 50;
const int num_iterations = 3000;

int ant_solution[num_ants][num_vertices];
int ant_path[num_ants][num_vertices];  // Stores the path each ant takes through the construction graph (for Daemon actions)
int ant_pos[num_ants];

const float rho = 0.5;
const float q = 0.5;
const float alpha = 1;
const float beta = 1;

mt19937 seed;
uniform_real_distribution<float> uni(0, 1);
uniform_int_distribution<int> random_vertex(0, num_vertices - 1);

int f(int * x) {
	int max = 0;
	for (int i = 0; i < num_vertices; i++) {
		if (x[i] > max) {
			max = x[i];
		}
	}
	return max + 1;
}

void initialise_pheromones() {
	for (int i = 0; i < num_vertices; i++) {
		for (int j = 0; j < num_vertices; j++) {
			tau[i][j] = 1 - adj_matrix[i][j];
		}
	}
}

void initialise_ants() {
	for (int i = 0; i < num_ants; i++) {
		for (int j = 0; j < num_vertices; j++) {
			ant_solution[i][j] = -1;
		}
	}
}

bool valid(int v, int c, int * col) {  // Returns whether or not vertex v can be coloured colour c in colouring col (legally)
	for (int i = 0; i < num_vertices; i++) {
		if (adj_matrix[v][i] == 1 && col[i] == c) {  // If v is adjacent to some vertex i coloured c then this is not a valid assignment
			return false;
		}
	}
	return true;  // If no such i can be found then the assignment is valid
}

float eta(int ant, int v) {
	// Returns the heuristic value of v in colouring ant
	vector<int> neighbouring_colours;
	for (int i = 0; i < num_vertices; i++) {
		if (adj_matrix[v][i] == 1 && ant_solution[ant][i] != -1 && find(neighbouring_colours.begin(), neighbouring_colours.end(), ant_solution[ant][i]) == neighbouring_colours.end()) {
			neighbouring_colours.push_back(ant_solution[ant][i]);
		}
	}
	return neighbouring_colours.size();
}

int main(){
	int best[num_vertices];
	int f_best;
	float weight[num_vertices];
	initialise_pheromones();
	int v;
	for (int t = 0; t < num_iterations; t++) {
		// Initialise ants (to a random vertex)
		initialise_ants();
		for (int i = 0; i < num_ants; i++) {
			ant_pos[i] = random_vertex(seed);
			int c = 0;
			while (true) {
				if (valid(ant_pos[i], c, ant_solution[i])) {
					ant_solution[i][ant_pos[i]] = c;
					break;
				}
				c++;
			}
		}
		// Build ant solutions
		for (int i = 0; i < num_vertices; i++) {
			for (int a = 0; a < num_ants; a++) {
				// Select next vertex with transition rule
				float weight_sum = 0;
				for (int u = 0; u < num_vertices; u++) {
					if (ant_solution[a][u] != -1) {
						weight[u] = 0;
					} else {
						weight[u] = pow(tau[ant_pos[a]][u], alpha) * pow(eta(a, u) + 1, beta);
						weight_sum += weight[u];
					}
				}
				if (uni(seed) < q) {
					v = 0;
					for (int u = 1; u < num_vertices; u++) {
						if (weight[u] > weight[v]) {
							v = u;
						}
					}
				} else {
					float r = uni(seed);
					v = -1;
					do {
						v++;
						r -= weight[v] / weight_sum;
					} while (r > 0);
				}
				// Colour vertex lowest colour possible
				int c = 0;
				while (true) {
					if (valid(v, c, ant_solution[a])) {
						ant_solution[a][v] = c;
						break;
					}
					c++;
				}
				// Apply online stage by stage pheromone update
				tau[ant_pos[a]][v] = (1 - rho) * tau[ant_pos[a]][v] + rho;
				ant_path[a][i] = ant_pos[a];
				ant_pos[a] = v;
			}
		}
		// Find best solution (and update best if it's better
		int b = 0;
		int f_b = f(ant_solution[0]);
		for(int i = 1; i < num_ants; i++){
			int f_i = f(ant_solution[i]);
			if (f_i < f_b) {
				b = i;
				f_b = f_i;
			}
		}
		if (t == 0 || f_b < f_best) {
			copy(begin(ant_solution[b]), end(ant_solution[b]), begin(best));
			f_best = f_b;
		}
		// Apply Daemon actions
		for (int i = 0; i < num_vertices - 1; i++) {
			tau[ant_path[b][i]][ant_path[b][i + 1]] += rho * f_b;
		}
		// Apply pheromone evaporation
		for (int i = 0; i < num_vertices; i++) {
			for (int j = 0; j < num_vertices; j++) {
				tau[i][j] *= 1 - rho;
			}
		}
	}
	// Return best solution
	for (int i = 0; i < num_vertices; i++) {
		cout << best[i] << " ";
	}
	cout << endl << "Number of colours: " << f_best;
	return 0;
}