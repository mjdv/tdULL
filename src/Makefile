
CPPFLAGS=-std=c++17 -O3 -Wall -g -Wno-sign-compare -march=native

N=001
F=001
T=199

INPUT_DIR=../input/exact
OUTPUT_DIR=../output/exact

all: set_trie_test graph_test treedepth_test main verify generate_exact_cache centrality_test heuristic

instance: main verify
	./main $(INPUT_DIR)/exact_$(N).gr $(OUTPUT_DIR)/exact_$(N).tree
	./verify $(INPUT_DIR)/exact_$(N).gr $(OUTPUT_DIR)/exact_$(N).tree

until: main verify
	bash from_to.sh $(T)

fromto: main verify
	bash from_to.sh $(F) $(T)

timing: main verify
	mkdir -p timings
	echo "fn,treedepth,time,error" > timings/computed_treedepths_$(shell git rev-parse --abbrev-ref HEAD)_$(shell git rev-parse --short HEAD).csv
	bash from_to.sh 001 199 2>> timings/computed_treedepths_$(shell git rev-parse --abbrev-ref HEAD)_$(shell git rev-parse --short HEAD).csv > timings/computed_treedepths_$(shell git rev-parse --abbrev-ref HEAD)_$(shell git rev-parse --short HEAD).log


main: main.o graph.o separator.o set_trie.o exact_cache.o centrality.o
	g++ -o $@ $^

heuristic: main_heuristic.o graph.o separator.o set_trie.o exact_cache.o centrality.o
	g++ -o $@ $^

set_trie_test: set_trie_test.o set_trie.o
	g++ -o $@ $^

graph_test: graph_test.o graph.o separator.o
	g++ -o $@ $^

treedepth_test: treedepth_test.o graph.o separator.o set_trie.o exact_cache.o centrality.o
	g++ -o $@ $^

centrality_test: centrality_test.o graph.o centrality.o
	g++ -o $@ $^

generate_exact_cache: graph.o separator.o generate_exact_cache.o set_trie.o exact_cache.o centrality.o
	g++ -o $@ $^

verify: verify.o
	g++ -o $@ $^

clean:
	rm set_trie_test *.o graph_test main *.d treedepth_test generate_exact_cache verify || true

exact_cache.o : exact_caches/exact_cache_8.bin
exact_caches/exact_cache_8.bin: exact_caches/exact_cache_8.bin.gz
	gunzip < $^ > $@



%.o : %.cpp
	g++ $(CPPFLAGS) -c $< -o $@
	g++ -MM $(CPPFLAGS) $*.cpp > $*.d

-include *.d
