
CPPFLAGS=-std=c++17 -O3 -Wall -g -Wno-sign-compare

N=001
F=001
T=199

INPUT_DIR=../input/exact
OUTPUT_DIR=../output

all: set_trie_test graph_test treedepth_test main verify generate_exact_cache

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

generate_exact_cache: graph.o generate_exact_cache.o set_trie.o exact_cache.o
	g++ -o $@ $^

verify: verify.o
	g++ -o $@ $^

clean:
	rm set_trie_test *.o graph_test main *.d treedepth_test generate_exact_cache

%.o : %.cpp
	g++ $(CPPFLAGS) -c $< -o $@
	g++ -MM $(CPPFLAGS) $*.cpp > $*.d

-include *.d
