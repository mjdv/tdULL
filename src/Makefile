
CPPFLAGS=-std=c++17 -O3 -Wall -g -Wno-sign-compare

N=001
INPUT_DIR=../input/exact
OUTPUT_DIR=../output

all: set_trie_test graph_test treedepth_test main

instance: main
	./main $(INPUT_DIR)/exact_$(N).gr $(OUTPUT_DIR)/exact_$(N).tree

until: main
	bash until.sh $N

main: main.o graph.o set_trie.o
	g++ -o $@ $^

set_trie_test: set_trie_test.o set_trie.o
	g++ -o $@ $^

graph_test: graph_test.o graph.o
	g++ -o $@ $^

treedepth_test: treedepth_test.o graph.o set_trie.o
	g++ -o $@ $^

clean:
	rm set_trie_test *.o graph_test main *.d treedepth_test

%.o : %.cpp
	g++ $(CPPFLAGS) -c $< -o $@
	g++ -MM $(CPPFLAGS) $*.cpp > $*.d

-include *.d
