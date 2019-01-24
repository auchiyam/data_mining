#ifndef LSH_CLUSTER
#define LSH_CLUSTER

class cluster_tree {
private:
	int max_nodes;
public:
	int length_nodes;
	cluster_tree** next;
	int* vals;
	int cluster_num;

	cluster_tree();
	void add_node(int hash_val, cluster_tree* next_node);
};

int lsh(int, int, double*, int, double, double**, int*, int**, int**, cluster_tree**);
int* find_hash_val(double**, double, int, int, double*, int);
int add_hash_node(cluster_tree*, int*, int, int*);
void add_cluster(double***, int*, int, int, double*, int, int, int**);
void adjust_data(double**, int*, int, int, double*, int, int**, int**);

#endif
