#ifndef LSH_SEARCH
#define LSH_SEARCH

int lsh_search(int dim, int ndata, double* data, int m, double w, double** h, int nclusters, int* cluster_start, int* cluster_size, cluster_tree* cluster_hashval, double* query, double* result);
int find_cluster(cluster_tree* cl_tr, int* hash_val, int m);

#endif
