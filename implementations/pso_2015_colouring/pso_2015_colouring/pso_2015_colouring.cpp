#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <array>
#include <chrono>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>

// Running time (population 50, 3000 iterations):
	// 500 vertices: 246.89 secs
		// Couldn't colour with 200 colours (my_coa_basic coloured 1000 with 143 colours)

using namespace std;
using namespace boost::random;

const string graph_directory = "C:/Users/taydo/OneDrive/Documents/computer_science/year3/project/implementations/graphs/";

int k;
const int num_vertices = 250;
int m = 0;

int adj_matrix[num_vertices][num_vertices];/* = {
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

const int num_particles = 50;
const int num_iterations = 3000;
chrono::time_point<chrono::steady_clock> start;
const auto duration = chrono::minutes{2};

const float w = 0.05;
const float c1 = 7;
const float c2 = 0.03;

int particles[num_particles][num_vertices];  // The population of particles
int particles_old[num_particles][num_vertices];  // The previous generation of particles
int p[num_particles][num_vertices];  // The personal best for each particle
int g[num_vertices];  // The global best particle found

mt19937 seed;
uniform_real_distribution<float> uni(0, 1);
uniform_int_distribution<int> random_colour;//(0, k - 1);
uniform_int_distribution<int> random_vertex(0, num_vertices - 1);

void make_graph(float edge_probability) {  // Populates adj_matrix with a random graph
	for (int i = 0; i < num_vertices; i++) {
		for (int j = 0; j < i; j++) {
			if (uni(seed) < edge_probability) {
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
			}
		}
	}
	else {
		cout << "Couldn't open file." << endl;
		exit(1);
	}
	file.close();
}

float f(int * x) {
	int num = 0;
	for (int i = 0; i < num_vertices; i++) {
		for (int j = 0; j < i; j++) {
			if (adj_matrix[i][j] == 1 && x[i] == x[j]) {
				num++;
			}
		}
	}
	return 1 - (num / (float)m);
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

float d(int * x, int * y) {
	int num = 0;
	for (int i = 0; i < num_vertices; i++) {
		if (x[i] != y[i]) {
			num++;
		}
	}
	return 1 - num / (float)num_vertices;
}

void generate_particle(int * particle) {
	for (int i = 0; i < num_vertices; i++) {
		particle[i] = random_colour(seed);
	}
}

bool valid(int v, int c, int * col) {  // Returns whether or not vertex v can be coloured colour c in colouring col (legally)
	for (int i = 0; i < num_vertices; i++) {
		if (adj_matrix[i][v] == 1 && col[i] == c) {  // If v is adjacent to some vertex i coloured c then this is not a valid assignment
			return false;
		}
	}
	return true;  // If no such i can be found then the assignment is valid
}

int colouring[num_vertices];
int chromatic_bound() {
	int ret = 0;
	for (int i = 0; i < num_vertices; i++) {
		int c = 0;
		while (true) {
			if (valid(i, c, colouring)) {
				colouring[i] = c;
				if (c > ret) {
					ret = c;
				}
				break;
			}
			c++;
		}
	}
	return ret + 1;
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

bool find_colouring() {
	// Generate initial population
	float g_fitness = 0;
	for (int i = 0; i < num_particles; i++) {
		generate_particle(particles[i]);
		copy(begin(particles[i]), end(particles[i]), begin(p[i]));
		if (f(particles[i]) > g_fitness) {
			copy(begin(particles[i]), end(particles[i]), begin(g));
			g_fitness = f(particles[i]);
		}
	}
	// Begin main loop
	while(chrono::duration_cast<chrono::minutes>(chrono::high_resolution_clock::now() - start) < duration){
		for (int i = 0; i < num_particles; i++) {
			float r1 = uni(seed);
			float r2 = uni(seed);
			float v_rand = w * d(particles[i], particles_old[i]);
			float v_p = c1 * r1 * d(particles[i], p[i]);
			float v_g = c2 * r2 * d(particles[i], g);
			float V = v_rand + v_p + v_g;
			float p_rand = v_rand / V; float p_p = v_p / V;  // No need to actually calculate p_g
			copy(begin(particles[i]), end(particles[i]), begin(particles_old[i]));
			for (int j = 0; j < num_vertices; j++) {
				float r = uni(seed);
				if (r <= p_rand) {
					particles[i][j] = random_colour(seed);
				}
				else if (r <= p_rand + p_p) {
					particles[i][j] = p[i][j];
				}
				else {
					particles[i][j] = g[j];
				}
			}
			if (f(particles[i]) > f(p[i])) {
				copy(begin(particles[i]), end(particles[i]), begin(p[i]));
				if (f(particles[i]) > f(g)) {
					copy(begin(particles[i]), end(particles[i]), begin(g));
				}
			}
		}

		if (num_conflicts(g) == 0) {
			copy(begin(g), end(g), begin(colouring));
			return true;
		}
	}
	return false;
}

int main(){  // Could probably precompute some fitness values (not as many as in coa but still could)
	//make_graph(0.5);
	read_graph("dsjc250.5.col");
	// Calculate number of edges
	for (int i = 0; i < num_vertices; i++) {
		for (int j = 0; j < i; j++) {
			m += adj_matrix[i][j];
		}
	}
	k = chromatic_bound();
	random_colour = uniform_int_distribution<int>(0, k - 1);
	bool found_colouring = true;
	start = chrono::high_resolution_clock::now();
	while (found_colouring) {
		cout << k + 1 << endl;
		found_colouring = find_colouring();
		k = num_colours(colouring) - 1;
		random_colour = uniform_int_distribution<int>(0, k - 1);
	}
	for (int i = 0; i < num_vertices; i++) {
		cout << colouring[i] << " ";
	}
	cout << endl << "Number of colours: " << num_colours(colouring) << endl;
	cout << "Number of conflicts: " << num_conflicts(colouring) << endl;
	//cout << "Time taken (seconds): " << chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start).count() / (float)1000000 << endl;
	return 0;
}
