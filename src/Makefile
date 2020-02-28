CPPFLAGS=-std=c++17 -O3 -Wall

set_trie_test: set_trie_test.o set_trie.o
	g++ -o $@ $^

clean:
	rm set_trie_test *.o

%.o : %.cpp
	g++ $(CPPFLAGS) -c $< -o $@
	g++ -MM $(CPPFLAGS) $*.cpp > $*.d

-include *.d
