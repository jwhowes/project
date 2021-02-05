#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>

using namespace std;
using namespace boost::random;

const string graph_directory = "C:/Users/taydo/OneDrive/Documents/computer_science/year3/project/implementations/graphs/";
const string results_directory = "C:/Users/taydo/OneDrive/Documents/computer_science/year3/project/implementations/results/";

const int num_vertices = 300;
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

int adj_list[num_vertices][num_vertices];/* = {
	{1, 4, 5, 0, 0, 0, 0, 0, 0, 0},
	{0, 2, 6, 0, 0, 0, 0, 0, 0, 0},
	{1, 3, 7, 0, 0, 0, 0, 0, 0, 0},
	{2, 4, 8, 0, 0, 0, 0, 0, 0, 0},
	{0, 3, 9, 0, 0, 0, 0, 0, 0, 0},
	{0, 7, 8, 0, 0, 0, 0, 0, 0, 0},
	{1, 8, 9, 0, 0, 0, 0, 0, 0, 0},
	{2, 5, 9, 0, 0, 0, 0, 0, 0, 0},
	{3, 5, 6, 0, 0, 0, 0, 0, 0, 0},
	{4, 6, 7, 0, 0, 0, 0, 0, 0, 0}
};*/

int adj_list_length[num_vertices];// = { 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 };

int k;

const int num_iterations = 3000;
chrono::time_point<chrono::steady_clock> start;
const auto duration = chrono::minutes{5};

mt19937 seed;
uniform_real_distribution<float> uni(0, 1);
uniform_int_distribution<int> random_colour;//(0, k - 1);

int colouring[num_vertices];
int best_colouring[num_vertices];

int orientation[num_vertices][num_vertices];
int d_plus[num_vertices];
int d_minus[num_vertices];

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

int f1(int * col) {
	int num = 0;
	for (int i = 0; i < num_vertices; i++) {
		for (int j = 0; j < i; j++) {
			if (adj_matrix[i][j] == 1 && col[i] == col[j]) {
				num++;
			}
		}
	}
	return num;
}

int f2(int * col) {
	int num = 0;
	for (int i = 0; i < num_vertices; i++) {
		if (col[i] < 0) {
			num++;
		}
	}
	return num;
}

void t12(int * col) {
	for (int i = 0; i < num_vertices; i++) {
		for (int j = 0; j < i; j++) {
			if (adj_matrix[i][j] == 1 && col[i] == col[j] && col[i] != k && col[j] != k) {
				if (uni(seed) < 0.5) {
					col[i] = k;
				} else {
					col[j] = k;
				}
			}
		}
	}
}

int order[num_vertices];
int colour_counts[num_vertices];
void t21(int * col) {
	random_shuffle(begin(order), end(order));
	for (int v : order) {
		if (col[v] == k) {
			// Assign to col[v] the colour causing the fewest conflicts
			int c = 0;
			for (int j = 0; j < k; j++) {
				colour_counts[k] = 0;
			}
			for (int j = 0; j < num_vertices; j++) {
				if (adj_matrix[v][j] == 1 && col[j] != k) {
					colour_counts[col[j]]++;
				}
			}
			for (int j = 0; j < k; j++) {
				if (colour_counts[j] == 0) {
					c = j;
					break;
				}
				if (colour_counts[j] < colour_counts[c]) {
					c = j;
				}
			}
			col[v] = c;
		}
	}
}

void t13(int * col) {  // Transforms a colouring (partial or otherwise), into an oriented graph
	for (int i = 0; i < num_vertices; i++) {
		for (int j = i; j < num_vertices; j++) {
			if (adj_matrix[i][j] == 1) {
				if (col[i] < col[j]) {
					orientation[i][j] = 1;
				} else {
					orientation[i][j] = -1;
				}
				orientation[j][i] = -orientation[i][j];
			} else {
				orientation[j][i] = orientation[i][j] = 0;
			}
		}
	}
}

int class_size[num_vertices];

bool cmp_class_size(int i, int j) {
	return class_size[i] > class_size[j];
}

void t32(int * col) {
	//copy(&col, &col + num_vertices, begin(d_minus));
	for (int i = 0; i < num_vertices; i++) {
		class_size[i] = 0;
	}
	for (int i = 0; i < num_vertices; i++) {
		class_size[d_minus[i]]++;
	}
	sort(begin(order), end(order), cmp_class_size);
	int curr = 0;
	col[order[0]] = 0;
	for (int i = 1; i < num_vertices; i++) {
		if (d_minus[order[i]] != d_minus[order[i - 1]] && curr < k) {
			curr++;
		}
		col[order[i]] = curr;
	}
}

int critical_vertices[num_vertices];
int num_critical;
int tabu_list[num_vertices][num_vertices];
int gamma[num_vertices][num_vertices];
uniform_int_distribution<int> random_L(0, 9);
const float tabu_tenure = 0.6;
const int I_t = 10;
const int tabucol_iterations = 100;
const int t_3_iterations = 1;

void get_critical_vertices(int * colouring) {
	num_critical = 0;
	for (int i = 0; i < num_vertices; i++) {
		if (gamma[i][colouring[i]] > 0) {
			critical_vertices[num_critical] = i;
			num_critical++;
		}
	}
}

int uncoloured[num_vertices];
int num_uncoloured;
void get_uncoloured(int * col) {
	num_uncoloured = 0;
	for (int i = 0; i < num_vertices; i++) {
		if (col[i] == k) {
			uncoloured[num_uncoloured] = i;
			num_uncoloured++;
		}
	}
}

void inline update_gamma(int v, int c, int * colouring) {
	for (int i = 0; i < num_vertices; i++) {
		if (adj_matrix[i][v] == 1) {
			gamma[i][c]++;
			gamma[i][colouring[v]]--;
		}
	}
}

bool make_move_1(int t, int * colouring) {
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
				tabu_list[v][c] = t + random_L(seed) + tabu_tenure * num_critical;
				return true;
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
		tabu_list[best_v][best_c] = t + random_L(seed) + tabu_tenure * num_critical;
	}
	return false;
}

bool make_move_2(int t, int * colouring) {
	int F = f2(colouring);
	int temp_col[num_vertices];
	bool initial = true;
	int best_v; int best_c; int best_f;
	get_uncoloured(colouring);
	for (int i = 0; i < num_uncoloured; i++) {
		int v = uncoloured[i];
		for (int c = 0; c < k; c++) {
			int f = F - 1;
			// Count vertices adjacent to v coloured c
			for (int j = 0; j < num_vertices; j++) {
				if (adj_matrix[v][j] == 1 && colouring[j] == c) {
					f++;
				}
			}
			if (f < F) {
				// I don't think we need to worry about looping over N(v) here
				colouring[v] = c;
				tabu_list[v][c] = tabu_tenure * num_uncoloured + random_L(seed);
				return true;
			}
			else if (tabu_list[v][c] <= t && (initial || f < best_f)) {
				initial = false;
				best_f = f; best_v = v; best_c = c;
			}
		}
	}
	if (!initial) {
		tabu_list[best_v][best_c] = tabu_tenure * num_uncoloured + random_L(seed);
		colouring[best_v] = best_c;
		for (int j = 0; j < num_vertices; j++) {
			if (adj_matrix[best_v][j] == 1 && colouring[j] == best_c) {
				colouring[j] = k;
			}
		}
		return false;
	}
}

void tabucol_1() {
	for (int i = 0; i < num_vertices; i++) {
		for (int j = 0; j < k; j++) {
			tabu_list[i][j] = 0;
			gamma[i][j] = 0;
		}
	}
	for (int i = 0; i < num_vertices; i++) {
		for (int j = 0; j < num_vertices; j++) {
			if (adj_matrix[i][j] == 1) {
				gamma[i][colouring[j]]++;
			}
		}
	}
	int t = 0;
	int num_without_increase = 0;
	while (t < tabucol_iterations) {//(num_without_increase < I_t) {
		if (make_move_1(t, colouring)) {
			num_without_increase++;
		} else {
			num_without_increase = 0;
		}
		if (f1(colouring) == 0) {
			return;
		}
		t++;
	}
}

void tabucol_2() {
	for (int i = 0; i < num_vertices; i++) {
		for (int j = 0; j < k; j++) {
			tabu_list[i][j] = 0;
		}
	}
	int t = 0;
	int num_without_increase = 0;
	while (t < tabucol_iterations) {//(num_without_increase < I_t) {
		if (make_move_2(t, colouring)) {
			num_without_increase++;
		} else {
			num_without_increase = 0;
		}
		if (f2(colouring) == 0) {
			return;
		}
		t++;
	}
}

int populate_distances(int orientation[num_vertices][num_vertices], int * d_plus, int * d_minus) {
	int indegree[num_vertices];
	int outdegree[num_vertices];
	int queue[num_vertices];
	int queue_length;
	int lambda = 1;
	for (int i = 0; i < num_vertices; i++) {
		d_minus[i] = 1;
		d_plus[i] = 1;
		indegree[i] = 0;
		outdegree[i] = 0;
	}
	queue_length = 0;
	for (int i = 0; i < num_vertices; i++) {
		for (int j = 0; j < adj_list_length[i]; j++) {
			if (orientation[i][adj_list[i][j]] == -1) {
				indegree[i]++;
			}
			else {
				outdegree[i]++;
			}
		}
	}
	for (int i = 0; i < num_vertices; i++) {
		if (outdegree[i] == 0) {
			queue[queue_length] = i;
			queue_length++;
		}
	}
	int num_checked = 0;
	while (num_checked < queue_length) {
		int i = queue[num_checked];
		for (int j = 0; j < num_vertices; j++) {
			if (orientation[i][j] == -1) {
				outdegree[j] -= 1;
				if (d_plus[i] + 1 > d_plus[j]) {
					d_plus[j] = d_plus[i] + 1;
					if (d_plus[j] > lambda) {
						lambda = d_plus[j];
					}
				}
				if (outdegree[j] == 0) {
					queue[queue_length] = j;
					queue_length++;
				}
			}
		}
		num_checked++;
	}
	num_checked = 0;
	queue_length = 0;
	for (int i = 0; i < num_vertices; i++) {
		if (indegree[i] == 0) {
			queue[queue_length] = i;
			queue_length++;
		}
	}
	while (num_checked < queue_length) {
		int i = queue[num_checked];
		for (int j = 0; j < num_vertices; j++) {
			if (orientation[i][j] == 1) {
				indegree[j] -= 1;
				if (d_minus[i] + 1 > d_minus[j]) {
					d_minus[j] = d_minus[i] + 1;
				}
				if (indegree[j] == 0) {
					queue[queue_length] = j;
					queue_length++;
				}
			}
		}
		num_checked++;
	}
	return lambda;
}

int on_longest_path[num_vertices];
int on_longest_length;

void make_move_3(int orientation[num_vertices][num_vertices], int v, int j) {

}

void tabucol_3() {
	int orientation_temp[num_vertices][num_vertices];
	int d_p_temp[num_vertices];
	int d_m_temp[num_vertices];
	bool found[num_vertices];
	int lambda = populate_distances(orientation, d_plus, d_minus);
	for (int i = 0; i < num_vertices; i++) {
		tabu_list[i][0] = 0;
		tabu_list[i][1] = 0;
	}
	for (int t = 0; t < t_3_iterations; t++) {
		on_longest_length = 0;
		for (int i = 0; i < num_vertices; i++) {
			found[i] = false;
		}
		for (int i = 0; i < num_vertices; i++) {
			for (int j = i; j < num_vertices; j++) {
				if ((orientation[i][j] == 1 && d_minus[i] + d_plus[j] == lambda) || (orientation[i][j] == -1 && d_minus[j] + d_plus[i] == lambda)) {
					if (!found[i]) {
						found[i] = true;
						on_longest_path[on_longest_length] = i;
						on_longest_length++;
					}
					if (!found[j]) {
						found[j] = true;
						on_longest_path[on_longest_length] = j;
						on_longest_length++;
					}
				}
			}
		}
		bool move_made = false;
		bool initial = true;
		int b_lambda; int b_v; int b_j;
		for (int i = 0; i < on_longest_length; i++) {
			for (int j = 0; j < 2; j++) {
				int v = on_longest_path[i];
				copy(&orientation[0][0], &orientation[0][0] + num_vertices * num_vertices, &orientation_temp[0][0]);
				make_move_3(orientation_temp, v, j);
				int t_l = populate_distances(orientation_temp, d_p_temp, d_m_temp);
				if (t_l < lambda) {
					tabu_list[v][1 - j] = t + random_L(seed) + tabu_tenure * on_longest_length;
					copy(&orientation_temp[0][0], &orientation_temp[0][0] + num_vertices * num_vertices, &orientation[0][0]);
					copy(begin(d_p_temp), end(d_p_temp), begin(d_plus));
					copy(begin(d_m_temp), end(d_m_temp), begin(d_minus));
					lambda = t_l;
					move_made = true;
					break;
				} else if (tabu_list[v][j] <= t && (initial || t_l < b_lambda)) {
					initial = false;
					b_lambda = t_l; b_v = v; b_j = j;
				}
			}
			if (move_made) {
				break;
			}
		}
		if (!move_made && !initial) {
			make_move_3(orientation, b_v, b_j);
			lambda = b_lambda;
			tabu_list[b_v][1 - b_j] = t + random_L(seed) + tabu_tenure * on_longest_length;
		}
	}
}

int num_colours(int * x) {
	bool found[num_vertices];
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

int t;
ofstream ofile;
bool find_colouring() {
	for (int i = 0; i < num_vertices; i++) {
		colouring[i] = random_colour(seed);
	}
	while(t < num_iterations){
	//while (chrono::duration_cast<chrono::minutes>(chrono::high_resolution_clock::now() - start) < duration) {
		tabucol_1();
		if (t % 10 == 0) {
			ofile << num_colours(colouring) << endl;
		}
		t++;
		if (f1(colouring) == 0) {
			copy(begin(colouring), end(colouring), begin(best_colouring));
			return true;
		}
		t13(colouring);
		tabucol_3();
		t32(colouring);
		tabucol_2();
		t21(colouring);
	}
	return false;
}

bool valid(int v, int c, int * col) {  // Returns whether or not vertex v can be coloured colour c in colouring col (legally)
	for (int i = 0; i < adj_list_length[v]; i++) {
		if (col[adj_list[v][i]] == c) {  // If v is adjacent to some vertex i coloured c then this is not a valid assignment
			return false;
		}
	}
	return true;  // If no such i can be found then the assignment is valid
}

int chromatic_bound() {
	int ret = 0;
	for (int i = 0; i < num_vertices; i++) {
		int c = 0;
		while (true) {
			if (valid(i, c, best_colouring)) {
				best_colouring[i] = c;
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

int main(){
	cout << "VSS\n";
	read_graph("flat300_26.col");
	ofile.open(results_directory + "flat300_26_vss.txt");
	// Populate order array
	for (int i = 0; i < num_vertices; i++) {
		order[i] = i;
	}
	k = chromatic_bound() - 1;
	// Generate initial solution
	bool found_colouring = true;
	start = chrono::high_resolution_clock::now();
	while (found_colouring) {
		random_colour = uniform_int_distribution<int>(0, k - 1);
		found_colouring = find_colouring();
		k = num_colours(best_colouring) - 1;
	}
	ofile.close();
	for (int i = 0; i < num_vertices; i++) {
		cout << best_colouring[i] << " ";
	}
	cout << endl << "Number of colours: " << num_colours(best_colouring) << endl;
	cout << endl << "Number of conflicts: " << f1(best_colouring) << endl;
	return 0;
}