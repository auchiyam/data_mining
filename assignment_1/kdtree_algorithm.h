#ifndef KDTREE_CLUSTERING
#define KDTREE_CLUSTERING

//kdtree_clustering
int kdtree(int, int, double*, int, int*, int*, double**, double**);
double* find_centroid(int, double*, int, int);
double* find_variance(int, double*, int, int, double*);
double* calculate_criteria(int, double*, double*);
int bipartition(int, double*, int, int, double*, int[2], int[2]);
double* calculate_boundary(int, double*, int, int);
int assign_values(int*, int[2], int*, int[2], double**, double*[2], double**, double*[2], int*, int*);
double* get_data(double*, int, int);

//kdtree_search
int kdsearch(int, int, double*, int, int*, int*, double**, double*, double*);
int find_closest_cluster(int, int, double*, double**);
int get_smallest_distance(int, int, double*, double*, int*, int*, double*, int*, int*);
int eliminate_clusters(int, int, int, double**, double*, bool*, int*);
double get_cluster_distance(int, int, double*, double**);
double get_distance(int, double*, double*);

#endif
