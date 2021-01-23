#include <iostream>
#include <algorithm>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>

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

const int k = 2;

const int num_iterations = 3000;

mt19937 seed;
uniform_real_distribution<float> uni(0, 1);
uniform_int_distribution<int> random_colour(0, k - 1);

int colouring[num_vertices];

int orientation[num_vertices][num_vertices];
int d_plus[num_vertices];
int d_minus[num_vertices];

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
int colour_counts[k];
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
const float lambda = 0.6;
const int I_t = 10;
const int tabucol_iterations = 100;

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
				tabu_list[v][c] = t + random_L(seed) + lambda * num_critical;
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
		tabu_list[best_v][best_c] = t + random_L(seed) + lambda * num_critical;
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
				tabu_list[v][c] = lambda * num_uncoloured + random_L(seed);
				return true;
			}
			else if (tabu_list[v][c] <= t && (initial || f < best_f)) {
				initial = false;
				best_f = f; best_v = v; best_c = c;
			}
		}
	}
	if (!initial) {
		tabu_list[best_v][best_c] = lambda * num_uncoloured + random_L(seed);
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

void tabucol_3() {

}

int main(){
	// Populate order array
	for (int i = 0; i < num_vertices; i++) {
		order[i] = i;
	}
	// Generate initial solution
	for (int i = 0; i < num_vertices; i++) {
		colouring[i] = random_colour(seed);
	}
	for (int t = 0; t < num_iterations; t++) {
		tabucol_1();
		t13(colouring);
		//tabucol_3();
		t32(colouring);
		//tabucol_2();
		//t21(colouring);
	}
	for (int i = 0; i < num_vertices; i++) {
		cout << colouring[i] << " ";
	}
	cout << endl << "Number of conflicts: " << f1(colouring) << endl;
	return 0;
}