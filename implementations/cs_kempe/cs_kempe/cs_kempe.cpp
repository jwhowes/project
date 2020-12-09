#define _USE_MATH_DEFINES
#define _SECURE_SCL 0

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

const int num_vertices = 250;
int adj_matrix[num_vertices][num_vertices];
int adj_list[num_vertices][num_vertices];
int adj_list_length[num_vertices];

const auto duration = chrono::minutes{2};

const int num_nests = 50;
const float pa = 0.1;
const float p = 0.1;

const float beta = 0.5;
const float alpha = 1;

const float sigma_p = pow((tgamma(1 + beta)*sin(M_PI*beta / 2)) / (tgamma((1 + beta) / 2)*beta*pow(2, (beta - 1) / 2)), 2 / beta);

struct Nest {
	int nest[num_vertices];
	int fitness;
};

Nest nests[num_nests];

mt19937 seed;
uniform_real_distribution<float> uni(0, 1);
uniform_int_distribution<int> random_vertex(0, num_vertices - 1);
uniform_int_distribution<int> random_colour;
normal_distribution<float> normal_q(0, 1);
normal_distribution<float> normal_p(0, sigma_p);

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


int classes[num_vertices];  // Are we completely confident this is bounded by NUM_VERTICES?
int num_classes;
int f(int * col) {
	num_classes = num_colours(col);
	int ret = 0;
	for (int i = 0; i < num_vertices; i++) {
		classes[col[i]]++;
	}
	for (int i = 0; i < num_classes; i++) {
		ret -= classes[i] * classes[i];
	}
	return ret;
}

int K[num_vertices];
int length_k;
void kempe_chain(int c, int d, int v, int * nest) {
	K[0] = v;
	length_k = 1;
	int i = 0;
	while (i < length_k) {
		v = K[i];
		for (int j = 0; j < adj_list_length[v]; j++) {
			if ((nest[adj_list[v][j]] == c || nest[adj_list[v][j]] == d) && find(begin(K), begin(K) + length_k, adj_list[v][j]) == begin(K) + length_k) {
				K[length_k] = adj_list[v][j];
				length_k++;
			}
		}
		i++;
	}
}

float levy() {
	float p = normal_p(seed);
	float q = normal_q(seed);
	return p / pow(abs(q), 1 / beta);
}

int levy_flight(int * nest) {
	int n = num_colours(nest);
	uniform_int_distribution<int> random_colour(0, n);
	float M = abs(alpha * levy()) + 1;
	if (M > num_vertices) {
		M = num_vertices;
	}
	for (int i = 0; i < M; i++) {
		// Perform a kempe chain interchange
		int v = random_vertex(seed);
		int c = nest[v];
		int d = random_colour(seed);
		if (d == c) {
			d = (d + 1) % n;
		}
		kempe_chain(c, d, v, nest);
		for (int j = 0; j < length_k; j++) {
			if (nest[K[j]] == c) {
				nest[K[j]] = d;
			}else {
				nest[K[j]] = c;
			}
		}
	}
	return f(nest);
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
void get_cuckoo(int * nest) {
	for (int i = 0; i < num_vertices; i++) {
		nest[i] = -1;
	}
	random_shuffle(begin(order), end(order));
	for (int v : order) {
		int c = 0;
		while (!valid(v, c, nest)) {
			c++;
		}
		nest[v] = c;
	}
}

bool compare_nests(Nest & nest1, Nest & nest2) {
	return nest1.fitness < nest2.fitness;
}

int num_conflicts(int * nest) {
	int num = 0;
	for (int i = 0; i < num_vertices; i++) {
		for (int j = 0; j < i; j++) {
			if (adj_matrix[i][j] == 1 && nest[i] == nest[j]) {
				num++;
			}
		}
	}
	return num;
}

// I realise now that the tabucol (as it currently is) is useless. It's trying to reduce conflicts but all colourings are legal
// I should modify it to minimise f rather than num_conflicts

int critical_vertices[num_vertices];
int num_critical;
const int tabucol_iterations = 50;
int tabu_list[num_vertices][num_vertices];
int gamma[num_vertices][num_vertices];
uniform_int_distribution<int> random_L(0, 9);
int k;
const float lambda = 0.6;

void get_critical_vertices(int * colouring) {
	num_critical = 0;
	for (int i = 0; i < num_vertices; i++) {
		if (gamma[i][colouring[i]] > 0) {
			critical_vertices[num_critical] = i;
			num_critical++;
		}
	}
}

void inline update_gamma(int v, int c, int * colouring) {
	for (int i = 0; i < adj_list_length[v]; i++) {
		gamma[adj_list[v][i]][c]++;
		gamma[adj_list[v][i]][colouring[v]]--;
	}
}

void tabucol_make_move(int t, int * colouring) {
	bool initial = true;
	int best_v; int best_c; int best_d;
	get_critical_vertices(colouring);
	for (int i = 0; i < num_critical; i++) {
		int v = critical_vertices[i];
		for (int c = 0; c < k; c++) {
			int d = gamma[v][c] - gamma[v][colouring[v]];
			if (d < 0) {
				update_gamma(v, c, colouring);
				colouring[v] = c;
				tabu_list[v][c] = t + random_L(seed) + lambda * num_critical;
				return;
			}
			else if (tabu_list[v][c] <= t && (initial || d < best_d)) {
				initial = false;
				best_d = d; best_v = v; best_c = c;
			}
		}
	}
	if (!initial) {
		update_gamma(best_v, best_c, colouring);
		colouring[best_v] = best_c;
		tabu_list[best_v][best_c] = t + random_L(seed) + lambda * num_critical;
	}
}

int tabucol(int * colouring) {
	k = num_colours(colouring);
	for (int i = 0; i < num_vertices; i++) {
		for (int j = 0; j < k; j++) {
			tabu_list[i][j] = 0;
			gamma[i][j] = 0;
		}
	}
	for (int i = 0; i < num_vertices; i++) {
		for (int j = 0; j < adj_list_length[i]; j++) {
			gamma[i][colouring[adj_list[i][j]]]++;
		}
	}
	for (int t = 0; t < tabucol_iterations; t++) {
		tabucol_make_move(t, colouring);
		if (num_conflicts(colouring) == 0) {
			return f(colouring);
		}
	}
	return f(colouring);
}

int main(){
	read_graph("dsjc250.5.col");
	for (int i = 0; i < num_vertices; i++) {
		order[i] = i;
	}
	for (int i = 0; i < num_nests; i++) {
		get_cuckoo(nests[i].nest);
		nests[i].fitness = f(nests[i].nest);
	}
	int nest_temp[num_vertices];
	auto start = chrono::high_resolution_clock::now();
	int t = 0;
	while (chrono::duration_cast<chrono::minutes>(chrono::high_resolution_clock::now() - start) < duration) {
		t++;
		for (int n = 0; n < num_nests; n++) {
			copy(begin(nests[n].nest), end(nests[n].nest), begin(nest_temp));
			int l_f = levy_flight(nest_temp);
			if (uni(seed) < p || l_f < nests[n].fitness) {
				copy(begin(nest_temp), end(nest_temp), begin(nests[n].nest));
				nests[n].fitness = l_f;
			}
		}
		for (int n = 0; n < num_nests; n++) {
			copy(begin(nests[n].nest), end(nests[n].nest), begin(nest_temp));
			int t_f = tabucol(nest_temp);
			if (t_f < nests[n].fitness) {
				copy(begin(nest_temp), end(nest_temp), begin(nests[n].nest));
				nests[n].fitness = t_f;
			}
		}
		sort(begin(nests), end(nests), compare_nests);
		for (int i = 0; i < num_nests * pa; i++) {
			get_cuckoo(nests[num_nests - i - 1].nest);
			nests[num_nests - i - 1].fitness = f(nests[num_nests - i - 1].nest);
		}
	}
	sort(begin(nests), end(nests), compare_nests);
	for (int i = 0; i < num_vertices; i++) {
		cout << nests[0].nest[i] << " ";
	}
	cout << endl << "Number of colours: " << num_colours(nests[0].nest) << endl;
	cout << "Number of conflicts: " << num_conflicts(nests[0].nest) << endl;
	cout << "Number of iterations: " << t << endl;
	return 0;
}