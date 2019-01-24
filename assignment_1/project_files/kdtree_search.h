#ifndef KDTREE_SEARCH
#define KDTREE_SEARCH

//kdtree_search
int kdsearch(int, int, double*, int, int*, int*, double**, double*, double*);
int find_closest_cluster(int, int, double*, double**);
int get_smallest_distance(int, int, double*, double*, int*, int*, double*, int*, int*);
int eliminate_clusters(int, int, int, double**, double*, bool*, int*);
double get_cluster_distance(int, int, double*, double**);
double get_distance(int, double*, double*);

#endif
