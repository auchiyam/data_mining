#include <iostream>
#include "kdtree_search.h"
#include "kdtree_clustering.h"
#include <cmath>

int kdsearch(int dim, int ndata, double* data, int kk, int* cluster_start, int* cluster_size, double** cluster_bdry, double* query_pt, double* result_pt) {
	int k = find_closest_cluster(dim, kk, query_pt, cluster_bdry);
	double distance = -1;
	int closest_point = -1;
	int data_checked = 0;
	get_smallest_distance(dim, k, data, query_pt, cluster_start, cluster_size, &distance, &closest_point, &data_checked);
	int cluster_remaining = kk - 1;
	
	bool* clusters = new bool[kk];

	for (int i = 0; i < kk; i++) {
		clusters[i] = !(k == i);
	}

	while (cluster_remaining > 0) {
		eliminate_clusters(distance, dim, kk, cluster_bdry, query_pt, clusters, &cluster_remaining);
		for (int i = 0; i < kk; i++) {
			bool shorter_found = false;
			if (clusters[i]) {
				double curr_distance = -1;
				int curr_point = -1;

				get_smallest_distance(dim, k, data, query_pt, cluster_start, cluster_size, &curr_distance, &curr_point, &data_checked);

				if (curr_distance != -1 && curr_distance < distance) {
					distance = curr_distance;
					closest_point = curr_point;
					shorter_found = true;
				}

				clusters[i] = false;
				cluster_remaining--;
			}

			if (shorter_found) {
				break;
			}
		}
	}

	double* val = get_data(data, dim, closest_point);

	for (int i = 0; i < dim; i++) {
		result_pt[i] = val[i];
	}

	return data_checked;
}

int find_closest_cluster(int dim, int kk, double* query, double** cluster_bdry) {
	double d_min = -1;
	double distance = -1;
	int closest_cluster = -1;

	for (int i = 0; i < kk; i++) {
		distance = get_cluster_distance(dim, i, query, cluster_bdry);

		if (d_min == -1 || distance < d_min) {
			d_min = distance;
			closest_cluster = i;
		}
	}

	return closest_cluster;
}

int get_smallest_distance(int dim, int k, double* data, double* query, int* cluster_start, int* cluster_size, double* distance, int* closest_point, int* data_checked) {
	for (int i = 0; i < cluster_size[k]; i++) {
		double* data_point = get_data(data, dim, cluster_start[k] + i);

		double d = get_distance(dim, data_point, query);

		if (*distance == -1 || d < *distance) {	
			*distance = d;
			*closest_point = cluster_start[k] + i;
		}

		*data_checked += 1;
	}

	return 0;
}

int eliminate_clusters(int distance, int dim, int kk, double** cluster_bdry, double* query, bool* clusters, int* clusters_remaining) {
	for (int i = 0; i < kk; i++) {
		if (clusters[i]) {
			double d = get_cluster_distance(dim, i, query, cluster_bdry);

			if (d < distance) {
				clusters[i] = false;
				clusters_remaining--;
			}
		}
	}

	return 0;
}

double get_cluster_distance(int dim, int i, double* query, double** cluster_bdry) {
	double sum = 0;

	for (int j = 0; j < dim; j++) {
		double min = cluster_bdry[i][j*2];
		double max = cluster_bdry[i][j*2+1];

		if (query[j] < min) {
			sum += std::pow(min - query[j], 2);
		}
		if (query[j] > max) {
			sum += std::pow(max - query[j], 2);
		}
	}

	return std::sqrt(sum);
}

double get_distance(int dim, double* pt_A, double* pt_B) {
	double sum = 0;
	for (int i = 0; i < dim; i++) {
		sum += std::pow(pt_A[i] - pt_B[i], 2);
	}

	return std::sqrt(sum);
}
