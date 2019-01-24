kdtree: main.o kdcluster.o kdsearch.o
	g++ -o kdtree main.o kdcluster.o kdsearch.o

main.o: main.cpp
	g++ -c main.cpp

kdcluster.o: kdtree_clustering.cpp kdtree_clustering.h
	g++ -c kdtree_clustering.cpp -o kdcluster.o

kdsearch.o: kdtree_search.cpp kdtree_search.h
	g++ -c kdtree_search.cpp -o kdsearch.o

clean:
	rm -f *.o
