CPPFLAGS=-std=c++17 -O3 -Wall -g -Wno-sign-compare

all: set_trie_test graph_test main

main: main.o graph.o set_trie.o
	g++ -o $@ $^

set_trie_test: set_trie_test.o set_trie.o
	g++ -o $@ $^

graph_test: graph_test.o graph.o
	g++ -o $@ $^

clean:
	rm set_trie_test *.o graph_test main *.d

%.o : %.cpp
	g++ $(CPPFLAGS) -c $< -o $@
	g++ -MM $(CPPFLAGS) $*.cpp > $*.d

-include *.d
