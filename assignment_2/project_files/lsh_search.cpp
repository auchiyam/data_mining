#include<iostream>
#include<cmath>
#include<stdexcept>
#include"lsh_cluster.h"
#include"lsh_search.h"

int lsh_search(int dim, int ndata, double* data, int m, double w, double** h, int nclusters, int* cluster_start, int* cluster_size, cluster_tree* cluster_hashval, double* query, double* result) {
	int* hashval = find_hash_val(h, w, m, dim, query, 0);

	int cluster_num = find_cluster(cluster_hashval, hashval, m);

	double min_dist = -1;
	int res = -1;

	for (int i = cluster_start[cluster_num]; i < cluster_start[cluster_num] + cluster_size[cluster_num]; i++) {
		double dist = 0;
		for (int j = 0; j < dim; j++) {
			dist += std::pow(data[i*dim+j] - query[j], 2);
		}

		dist = std::sqrt(dist);

		if (min_dist == -1 || dist < min_dist) {
			min_dist = dist;
			res = i;
		}
	}

	for (int i = 0; i < dim; i++) {
		result[i] = data[res*dim+i];
	}

	return 0;
}

//assumtion: the queried value will always fit in one of the clusters
int find_cluster(cluster_tree* cl_tr, int* hash_val, int m) {
	cluster_tree* tmp = cl_tr;

	bool found = false;

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < cl_tr->length_nodes; j++) {
			if (hash_val[i] == tmp->vals[j]) {
				found = true;
				tmp = tmp->next[j];
				break;
			}
		}

		if (!found) {
			throw std::invalid_argument("Could not find any cluster the query belonged in.");
		}
	}

	return tmp->cluster_num;
}
