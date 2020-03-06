
CPPFLAGS=-std=c++17 -O3 -Wall -g -Wno-sign-compare

N=001
F=001
T=199

INPUT_DIR=../input/exact
OUTPUT_DIR=../output

all: set_trie_test graph_test treedepth_test main verify graph_generate

instance: main
	./main $(INPUT_DIR)/exact_$(N).gr $(OUTPUT_DIR)/exact_$(N).tree

until: main
	bash from_to.sh $(N)

fromto: main
	bash from_to.sh $(F) $(T)

main: main.o graph.o set_trie.o exact_cache.o
	g++ -o $@ $^

set_trie_test: set_trie_test.o set_trie.o
	g++ -o $@ $^

graph_test: graph_test.o graph.o
	g++ -o $@ $^

treedepth_test: treedepth_test.o graph.o set_trie.o exact_cache.o
	g++ -o $@ $^

graph_generate: graph.o graph_generate.o set_trie.o exact_cache.o
	g++ -o $@ $^

verify: verify.o
	g++ -o $@ $^

clean:
	rm set_trie_test *.o graph_test main *.d treedepth_test graph_generate

%.o : %.cpp
	g++ $(CPPFLAGS) -c $< -o $@
	g++ -MM $(CPPFLAGS) $*.cpp > $*.d

-include *.d
