#include <exception>

#include "treedepth.hpp"

std::string ExtractFileName(const std::string& str) {
  return str.substr(str.find_last_of("/") + 1);
}

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "Expecting 2 arguments." << std::endl;
    return 1;
  }

  std::ifstream input;
  input.open(argv[1], std::ios::in);
  LoadGraph(input);
  input.close();

  time_t start, end;
  time(&start);

  std::cout << "Calculating treedepth for " << ExtractFileName(argv[1])
            << std::endl;
  try {
    //    auto seperator_gen = SeparatorGenerator(full_graph_as_sub);
    //    size_t total_count = 0;
    //    while (seperator_gen.HasNext()) {
    //      total_count += seperator_gen.Next(100000).size();
    //      time(&end);
    //      std::cout << "Total number of separators is " << total_count
    //                << ". Speed is " << double(total_count) / difftime(end,
    //                start)
    //                << " seps / s.\n";
    //    }
    auto ap = full_graph.ArticulationPoints();
    assert(std::set<int>(ap.begin(), ap.end()).size() == ap.size());
    for (int v : ap) assert(v >= 0 && v < full_graph.N);
    // auto cc = full_graph.WithoutVertices(ap);
    // std::sort(cc.begin(), cc.end(),
    //          [](auto c1, auto c2) { return c1.N > c2.N; });
    std::cerr << ExtractFileName(argv[1]) << " has " << full_graph.N
              << " vertices, and " << ap.size() << " articulation points"
              << std::endl;
    auto core = full_graph.kCore(2);
    assert(core.size() == 1);
    auto core_ap = core[0].ArticulationPoints();
    std::cerr << ExtractFileName(argv[1]) << " 2core has " << core[0].N
              << " vertices, and " << core_ap.size() << " articulation points"
              << std::endl;
    if (core_ap.size()) {
      std::cerr
          << "Upon removing these articulation points from the 2core, the "
             "component sizes are:"
          << std::endl;
      auto cc = core[0].WithoutVertices(core_ap);
      std::sort(cc.begin(), cc.end(),
                [](auto c1, auto c2) { return c1.N > c2.N; });
      for (auto c : cc) std::cerr << c.N << ", ";
      std::cerr << std::endl;

      SeparatorGenerator sep_gen(cc[0]);
      auto seps = sep_gen.Next(100000).size();
      if (!sep_gen.HasNext())
        std::cerr << "Largest component has " << seps << " separators"
                  << std::endl;
      else
        std::cerr << "Largest component has more than " << seps
                  << " separators " << std::endl;
    }
    std::cerr << std::endl;
    return 0;
    auto [td, tree] = treedepth(full_graph);
    time(&end);
    std::cout << "Treedepth is: " << td << std::endl;
    std::cout << "Elapsed time is " << difftime(end, start) << " seconds.\n";

    std::ofstream output;
    output.open(argv[2], std::ios::out);
    output << td << std::endl;
    for (int parent : tree) output << parent << std::endl;
    output.close();
    std::cout << "Saved the tree to '" << argv[2] << "'" << std::endl;
    std::cerr << ExtractFileName(argv[1]) << "," << td << ","
              << difftime(end, start) << ", " << std::endl;
    return 0;
  } catch (std::exception& e) {
    time(&end);
    std::cout << "Failed! Encountered exception:\"" << e.what() << "\"."
              << std::endl;
    std::cerr << ExtractFileName(argv[1]) << "," << -1 << ","
              << difftime(end, start) << "," << e.what() << std::endl;
    return 1;
  }
}
