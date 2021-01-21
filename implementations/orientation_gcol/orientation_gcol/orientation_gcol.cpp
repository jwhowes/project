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

int adj_list[num_vertices][num_vertices] = {
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
};

int adj_list_length[num_vertices] = { 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 };

int tabu_list[num_vertices][2];

const int num_iterations = 100;

mt19937 seed;
uniform_int_distribution<int> random_L(0, 9);
const float tabu_tenure = 0.6;

void inline generate_initial_solution(int orientation[num_vertices][num_vertices]) {
	for (int i = 0; i < num_vertices; i++) {
		for(int j = i; j < num_vertices; j++){
			if (adj_matrix[i][j] == 1) {
				// j > i (hence orientation i <- j)
				orientation[i][j] = -1;
				orientation[j][i] = 1;
			} else {
				orientation[i][j] = orientation[j][i] = 0;
			}
		}
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
			} else {
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

void make_move(int orientation[num_vertices][num_vertices], int v, int j) {
	if (j == 0) {
		j = -1;
	}
	// Find all vertices on a longest path with an arc involving v and orientation j, then flip that orientation
	for (int i = 0; i < on_longest_length; i++) {
		int u = on_longest_path[i];
		if (orientation[u][v] == j) {
			orientation[u][v] = -j;
		}
	}
}

void copy_orientation(int src[num_vertices][num_vertices], int dst[num_vertices][num_vertices]) {
	for (int i = 0; i < num_vertices; i++) {
		for (int j = 0; j < num_vertices; j++) {
			dst[i][j] = src[i][j];
		}
	}
}

int main(){
	bool found[num_vertices];

	int orientation[num_vertices][num_vertices];
	int best_orientation[num_vertices][num_vertices];

	int d_plus[num_vertices];
	int d_minus[num_vertices];

	generate_initial_solution(orientation);
	//copy(&orientation[0][0], &orientation[0][0] + num_vertices * num_vertices, &best_orientation[0][0]);
	copy_orientation(orientation, best_orientation);

	int d_p_temp[num_vertices];
	int d_m_temp[num_vertices];
	int orientation_temp[num_vertices][num_vertices];

	for (int i = 0; i < num_vertices; i++) {
		tabu_list[i][0] = 0;
		tabu_list[i][1] = 0;
	}
	int lambda = populate_distances(orientation, d_plus, d_minus);
	int best_lambda = lambda;
	for (int t = 0; t < num_iterations; t++) {
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
		// Iterate over all possible moves
		bool move_made = false;
		bool initial = true;
		int b_lambda; int b_v; int b_j;
		for (int i = 0; i < on_longest_length; i++) {
			for (int j = 0; j < 2; j++) {
				int v = on_longest_path[i];
				copy(&orientation[0][0], &orientation[0][0] + num_vertices * num_vertices, &orientation_temp[0][0]);
				// Make the move on orientation_temp
				make_move(orientation_temp, v, j);
				int t_l = populate_distances(orientation_temp, d_p_temp, d_m_temp);
				if (t_l < best_lambda) {
					// Add (v, 1 - j) to TL
					tabu_list[v][1 - j] = t + random_L(seed) + tabu_tenure * on_longest_length;
					//copy(&orientation_temp[0][0], &orientation_temp[0][0] + num_vertices * num_vertices, &orientation[0][0]);
					//copy(&orientation[0][0], &orientation[0][0] + num_vertices * num_vertices, &best_orientation[0][0]);
					copy_orientation(orientation_temp, orientation);
					copy_orientation(orientation, best_orientation);
					best_lambda = lambda = t_l;
					copy(begin(d_p_temp), end(d_p_temp), begin(d_plus));
					copy(begin(d_m_temp), end(d_m_temp), begin(d_minus));
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
			// Make the move on orientation
			make_move(orientation, b_v, b_j);
			lambda = populate_distances(orientation, d_plus, d_minus);
			// Make the move tabu
			tabu_list[b_v][1 - b_j] = t + random_L(seed) + tabu_tenure * on_longest_length;
		}
	}
	for (int i = 0; i < num_vertices; i++) {
		for (int j = 0; j < num_vertices; j++) {
			cout << best_orientation[i][j] << " ";
		}
		cout << endl;
	}
	cout << best_lambda << endl;
	return 0;
}