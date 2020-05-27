#include <chrono>
#include <cmath>
#include <exception>

#include "treedepth.hpp"

std::string ExtractFileName(const std::string& str) {
  return str.substr(str.find_last_of("/") + 1);
}

int main(int argc, char** argv) {
  LoadGraph(std::cin);

  auto start = std::chrono::steady_clock::now();
  try {
    //    auto seperator_gen = SeparatorGenerator(full_graph);
    //    std::cerr << typeid(seperator_gen.done).name() << std::endl;
    //    size_t total_count = 0;
    //    auto last_print = std::chrono::steady_clock::now();
    //    while (seperator_gen.HasNext()) {
    //      total_count += seperator_gen.Next(10).size();
    //      if (std::chrono::duration<double>(std::chrono::steady_clock::now() -
    //                                        last_print)
    //              .count() > 15) {
    //        last_print = std::chrono::steady_clock::now();
    //        double time_elapsed = std::chrono::duration<double>(
    //                                  std::chrono::steady_clock::now() -
    //                                  start) .count();
    //        std::cerr << "Total time spent: " << 0.1 * std::round(10 *
    //        time_elapsed)
    //                  << "s. Total number of separators is " << total_count
    //                  << ". Speed is " << double(total_count) / time_elapsed
    //                  << " seps / s.\n";
    //      }
    //    }
    //    double time_elapsed =
    //        std::chrono::duration<double>(std::chrono::steady_clock::now() -
    //        start)
    //            .count();
    //    std::cerr << "Total number of separators is " << total_count
    //              << ". Total time is " << time_elapsed << "s. Speed is "
    //              << double(total_count) / time_elapsed << " seps / s.\n";
    //    return 0;

    auto [td, tree] = treedepth(full_graph);
    double time_elapsed =
        0.1 * std::round(10 * std::chrono::duration<double>(
                                  std::chrono::steady_clock::now() - start)
                                  .count());
    std::cerr << "Treedepth is: " << td << std::endl;
    std::cerr << "Elapsed time is " << time_elapsed << " seconds.\n";

    std::cout << td << std::endl;
    for (int parent : tree) std::cout << parent << std::endl;

    std::cerr << td << "," << time_elapsed << ", " << std::endl;
    return 0;
  } catch (std::exception& e) {
    double time_elapsed =
        0.1 * std::round(10 * std::chrono::duration<double>(
                                  std::chrono::steady_clock::now() - start)
                                  .count());
    std::cerr << "Failed! Encountered exception:\"" << e.what() << "\"."
              << std::endl;
    std::cerr << -1 << "," << time_elapsed << "," << e.what() << std::endl;
    return 1;
  }
}
