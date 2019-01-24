#include<iostream>
#include"lsh_cluster.h"

//roughly O(n*ncluster*m)
int lsh(int dim, int ndata, double* data, int m, double w, double** h, int* ncluster_ptr, int** cluster_start_ptr, int** cluster_size_ptr, cluster_tree** cluster_hashval_ptr) {
	//initialize
	*cluster_hashval_ptr = new cluster_tree();
	int max_clusters = 20;
	double** clusters = new double*[max_clusters];
	int* temp_cluster_size = new int[max_clusters];

	*ncluster_ptr = 0;

	for (int i = 0; i < max_clusters; i++) {
		clusters[i] = new double[ndata * dim];
		temp_cluster_size[i] = 0;
	}

	//O(n)
	for (int i = 0; i < ndata; i++) {
		int* hash_val = find_hash_val(h, w, m, dim, data, i); //O(dim*m)
		int cluster_num = add_hash_node(*cluster_hashval_ptr, hash_val, m, ncluster_ptr); //O(ncluster*m)

		add_cluster(&clusters, &max_clusters, dim, ndata, data, i, cluster_num, &temp_cluster_size); //O(dim * cluster_size[cluster_num]) if linked list, O(max_cluster) if array list
	}

	adjust_data(clusters, temp_cluster_size, dim, ndata, data, *ncluster_ptr, cluster_start_ptr, cluster_size_ptr); //O(n*dim)

	return 0;
}

int* find_hash_val(double** hash_func, double weight, int m, int dim, double* data, int loc) {
	int* hash_val = new int[m];
	for (int i = 0; i < m; i++) {
		double sum = 0;

		for (int j = 0; j < dim; j++) {
			sum += (data[loc*dim+j] * hash_func[i][j]);
		}

		hash_val[i] = (int)(sum / weight);
	}

	return hash_val;
}

int add_hash_node(cluster_tree* root, int* hash_val, int m, int* ncluster_ptr) {
	cluster_tree* iter = root;

	bool is_new = true;

	for (int i = 0; i < m; i++) {
		is_new = true;
		for (int j = 0; j < iter->length_nodes; j++) {
			if (hash_val[i] == iter->vals[j]) {
				is_new = false;
				iter = iter->next[j];
				break;
			}
		}

		if (is_new) {
			cluster_tree* new_node = new cluster_tree();

			iter->add_node(hash_val[i], new_node);

			iter = new_node;
		}
	}

	if (is_new) {
		iter->cluster_num = *ncluster_ptr;
		*ncluster_ptr += 1;
	}

	return iter->cluster_num;
}

void add_cluster(double*** clusters_ptr, int* max_clusters, int dim, int ndata, double* data, int loc, int cluster_num, int** cluster_size_ptr) {
	//the total clusters exceeded the maximum amount of clusters the program can store, so increase the size
	if (cluster_num >= *max_clusters) {
		*max_clusters *= 2;
		double** new_clusters = new double*[*max_clusters];
		int* new_cluster_size = new int[*max_clusters];

		for (int i = 0; i < cluster_num-1; i++) {
			double* new_data = new double[dim * ndata];

			for (int j = 0; j < dim; j++) {
				new_data[j] = *clusters_ptr[i][j];
			}

			new_clusters[i] = new_data;
			new_cluster_size[i] = *cluster_size_ptr[i];
		}

		for (int i = cluster_num; i < *max_clusters; i++) {
			new_clusters[i] = new double[dim * ndata];
			new_cluster_size[i] = 0;
		}

		delete clusters_ptr;
		delete cluster_size_ptr;

		*clusters_ptr = new_clusters;
		*cluster_size_ptr = new_cluster_size;
	}

	for (int i = 0; i < dim; i++) {
		(*clusters_ptr)[cluster_num][(*cluster_size_ptr)[cluster_num]*dim+i] = data[loc*dim+i];
	}

	(*cluster_size_ptr)[cluster_num] += 1;
}

void adjust_data(double** clusters, int* cluster_size_holder, int dim, int ndata, double* data, int nclusters, int** cluster_start, int** cluster_size) {
	int starting_point = 0;

	*cluster_start = new int[nclusters];
	*cluster_size = new int[nclusters];

	for (int i = 0; i < nclusters; i++) {
		(*cluster_start)[i] = starting_point;
		(*cluster_size)[i] = cluster_size_holder[i];

		for (int j = 0; j < cluster_size_holder[i]; j++) {
			for (int k = 0; k < dim; k++) {
				data[starting_point*dim+j*dim+k] = clusters[i][j*dim+k];
			}
		}

		starting_point += cluster_size_holder[i];
	}
}

cluster_tree::cluster_tree() {
	length_nodes = 0;
	max_nodes = 10;

	next = new cluster_tree*[max_nodes];
	vals = new int[max_nodes];
}

void cluster_tree::add_node(int hash_val, cluster_tree* next_node) {
	if (length_nodes == max_nodes) {
		max_nodes *= 2;

		cluster_tree** new_next = new cluster_tree*[max_nodes];
		int* new_vals = new int[max_nodes];

		for (int i = 0; i < length_nodes; i++) {
			new_next[i] = next[i];
			new_vals[i] = vals[i];
		}

		delete next;
		delete vals;

		next = new_next;
		vals = new_vals;
	}

	next[length_nodes] = next_node;
	vals[length_nodes] = hash_val;

	length_nodes++;
}
