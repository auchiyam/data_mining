lsh: main.o lsh_cluster.o lsh_search.o
	g++ -o lsh main.o lsh_cluster.o lsh_search.o

main.o: main.cpp
	g++ -c main.cpp

lsh_cluster.o: lsh_cluster.cpp lsh_cluster.h
	g++ -c lsh_cluster.cpp -o lsh_cluster.o

lsh_search.o: lsh_search.cpp lsh_search.h
	g++ -c lsh_search.cpp -o lsh_search.o

clean:
	rm -f *.o
