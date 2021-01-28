#define _SECURE_SCL 0

// Running time (k = 0.2V):
	// 100 vertices: 25.7828 secs

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

int k;
const int num_vertices = 100;

struct Partition {
	int partition[num_vertices][num_vertices];
	int partition_length[num_vertices];
	int fitness;
};

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

const int pop_size = 50;
const int num_iterations = 10;
chrono::time_point<chrono::steady_clock> start;
const auto duration = chrono::minutes{ 5 };
const float mutation_prob = 0.1f;

Partition population[pop_size];
Partition new_population[pop_size];

int colouring[num_vertices];

mt19937 seed;
uniform_real_distribution<float> uni(0, 1);
uniform_int_distribution<int> random_colour;
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

int f(Partition & x) {
	int num = 0;
	for (int c = 0; c < k; c++) {
		for (int i = 0; i < x.partition_length[c]; i++) {
			for (int j = 0; j < i; j++) {
				if (adj_matrix[x.partition[c][i]][x.partition[c][j]] == 1) {
					num++;
				}
			}
		}
	}
	return num;
}

Partition parents[2];
bool uncoloured[num_vertices];
void gpx(Partition & p1, Partition & p2, Partition & x) {
	for (int i = 0; i < num_vertices; i++) {
		uncoloured[i] = true;
	}
	int parent;
	for (int i = 0; i < num_vertices; i++) {
		uncoloured[i] = i;
	}
	for (int i = 0; i < k; i++) {
		copy(begin(p1.partition[i]), begin(p1.partition[i]) + p1.partition_length[i], begin(parents[0].partition[i]));
		parents[0].partition_length[i] = p1.partition_length[i];
		copy(begin(p2.partition[i]), begin(p2.partition[i]) + p2.partition_length[i], begin(parents[1].partition[i]));
		parents[1].partition_length[i] = p2.partition_length[i];
		x.partition_length[i] = 0;
	}
	for (int i = 0; i < k; i++) {
		parent = i % 2;
		int m = 0;
		for (int j = 1; j < k; j++) {
			if (parents[parent].partition_length[j] > parents[parent].partition_length[m]) {
				m = j;
			}
		}
		copy(begin(parents[parent].partition[m]), begin(parents[parent].partition[m]) + parents[parent].partition_length[m], begin(x.partition[i]));
		x.partition_length[i] = parents[parent].partition_length[m];
		for (int i = 0; i < parents[parent].partition_length[m]; i++) {
			uncoloured[parents[parent].partition[m][i]] = false;
			//uncoloured.erase(remove(uncoloured.begin(), uncoloured.end(), parents[parent].partition[m][i]), uncoloured.end());
		}
		parents[parent].partition_length[m] = 0;
		// Remove everything from x.partition[i] from parents[1-parent].partition
		for (int j = 0; j < x.partition_length[i]; j++) {
			for (int l = 0; l < k; l++) {
				auto pos = find(begin(parents[1 - parent].partition[l]), begin(parents[1 - parent].partition[l]) + parents[1 - parent].partition_length[l], x.partition[i][j]);
				//auto pos = find(parents[1-parent][l].begin(), parents[1-parent][l].end(), x[i][j]);
				if (pos != begin(parents[1 - parent].partition[l]) + parents[1 - parent].partition_length[l]) {
					copy(pos + 1, begin(parents[1 - parent].partition[l]) + parents[1 - parent].partition_length[l], pos);
					parents[1 - parent].partition_length[l]--;
					//parents[1-parent][l].erase(pos);
					break;
				}
			}
		}
	}
	for (int i = 0; i < num_vertices; i++) {
		if (uncoloured[i]) {
			int c = random_colour(seed);
			x.partition[c][x.partition_length[c]] = i;
			x.partition_length[c]++;
		}
	}
	/*for (int i = 0; i < uncoloured.size(); i++) {
		int c = random_colour(seed);
		x.partition[c][x.partition_length[c]] = uncoloured[i];
		x.partition_length[c]++;
	}*/
}

void generate_member(Partition & x) {
	for (int i = 0; i < num_vertices; i++) {
		int c = random_colour(seed);
		x.partition[c][x.partition_length[c]] = i;
		x.partition_length[c]++;
		//x[random_colour(seed)].push_back(i);
	}
}

void mutate(Partition & x) {  // Randomly moves a vertex from one colour class to another
	int c = random_colour(seed);
	while (x.partition_length[c] == 0) {
		c = (c + 1) % k;
	}
	int d = random_colour(seed);
	int v = x.partition[c][uniform_int_distribution<int>(0, x.partition_length[c] - 1)(seed)];
	auto pos = find(begin(x.partition[c]), begin(x.partition[c]) + x.partition_length[c], v);
	copy(pos + 1, begin(x.partition[c]) + x.partition_length[c], pos);
	x.partition_length[c]--;
	//x[c].erase(remove(x[c].begin(), x[c].end(), v), x[c].end());
	x.partition[d][x.partition_length[d]] = v;
	x.partition_length[d]++;
	//x[d].push_back(v);
}

bool compare_partitions(Partition & x, Partition & y) {
	return x.fitness < y.fitness;
}

void get_parents(int * p1, int * p2, int F) {
	float r = uni(seed);
	for (int i = 0; i < pop_size; i++) {
		r -= (population[pop_size - 1].fitness - population[i].fitness) / (float)F;
		if (r <= 0.0f) {
			*p1 = i;
			break;
		}
	}
	r = uni(seed);
	for (int i = 0; i < pop_size; i++) {
		r -= (population[pop_size - 1].fitness - population[i].fitness) / (float)F;
		if (r <= 0.0f) {
			*p2 = i;
			break;
		}
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

void partition_to_colouring(Partition & p) {
	for (int c = 0; c < k; c++) {
		for (int i = 0; i < p.partition_length[c]; i++) {
			colouring[p.partition[c][i]] = c;
		}
	}
}

int t;
bool find_colouring() {
	int p1; int p2;
	int F;
	for (int i = 0; i < pop_size; i++) {
		generate_member(population[i]);
		population[i].fitness = f(population[i]);
	}
	while(chrono::duration_cast<chrono::minutes>(chrono::high_resolution_clock::now() - start) < duration) {
		t++;
		sort(begin(population), end(population), compare_partitions);
		F = 0;
		for (int i = 0; i < pop_size; i++) {
			F += (population[pop_size - 1].fitness - population[i].fitness);
		}
		for (int i = 0; i < pop_size; i++) {
			get_parents(&p1, &p2, F);
			gpx(population[p1], population[p2], new_population[i]);
			if (uni(seed) < mutation_prob) {
				mutate(new_population[i]);
			}
		}
		for (int i = 0; i < pop_size; i++) {
			for (int j = 0; j < k; j++) {
				copy(begin(new_population[i].partition[j]), begin(new_population[i].partition[j]) + new_population[i].partition_length[j], begin(population[i].partition[j]));
			}
			copy(begin(new_population[i].partition_length), end(new_population[i].partition_length), begin(population[i].partition_length));
			population[i].fitness = f(population[i]);
			if (population[i].fitness == 0) {
				partition_to_colouring(population[i]);
				return true;
			}
		}
		//cout << chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start).count() << endl;
	}
	return false;
}

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

int main() {
	cout << "GA\n";
	//make_graph(0.5);
	read_graph("dsjc250.5.col");
	k = chromatic_bound() - 1;
	random_colour = uniform_int_distribution<int>(0, k - 1);
	bool found_colouring = true;
	start = chrono::high_resolution_clock::now();
	t = 0;
	while (found_colouring) {
		found_colouring = find_colouring();
		k = num_colours(colouring) - 1;
		random_colour = uniform_int_distribution<int>(0, k - 1);
	}
	for (int i = 0; i < num_vertices; i++) {
		cout << colouring[i] << " ";
	}
	cout << endl << "Number of colours: " << num_colours(colouring);
	cout << endl << "Number of iterations: " << t;
	return 0;
}