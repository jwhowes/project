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

int tau[num_vertices][num_vertices];

const int num_ants = 50;
const int num_iterations = 100;

int ant_solution[num_ants][num_vertices];
int ant_pos[num_ants];

int best[num_vertices];

mt19937 seed;
uniform_real_distribution<float> uni(0, 1);
uniform_int_distribution<int> random_vertex(0, num_vertices - 1);

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

int main(){
	initialise_pheromones();
	for (int t = 0; t < num_iterations; t++) {
		// Initialise ants (to a random vertex)
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
				// Colour vertex lowest colour possible
				// Apply online stage by stage pheromone update
			}
		}
		// Select best solution (update best if it's better)
		// Apply Daemon actions
		// Apply pheromone evaporation
	}
	// Return best solution
	return 0;
}