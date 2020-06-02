# tdULL

tdULL is a submission for the [PACE 2020 Treedepth
Challenge](https://pacechallenge.org/2020/) by Ruben Brokkelkamp, Raymond van
Venetië, Mees de Vries and Jan Westerdiep. It is a tool for computing the exact
[treedepth](https://en.wikipedia.org/wiki/Tree-depth) of a general graph.

## Usage

This tool uses the C++ Boost library. For detailed installation instructions,
see [the Boost website](https://www.boost.org/).

This tool uses Gregory Popovitch's [parallel hash
map](https://github.com/greg7mdp/parallel-hashmap). After cloning the
repository, before compiling, run `git init submodule` and then `git submodule
update` to include this code.

To build the tool, navigate to the `src` directory and run `make`. This
produces an executable called `main`, which can be used to compute treedepth:
it expects a graph in the [PACE Challenge input
format](https://pacechallenge.org/2020/td/) from `stdin`, and will output an
optimal treedepth decomposition of the graph in PACE Challenge output format to
`stdout`. Additionally, it will output some general information to `stderr`.

The executable `treedepth_test` (also produced by `make`) applies tdULL to a
selection of the public inputs of the PACE 2020 challenge. This test should
take no more than a few minutes.

## Status

By the nature of the challenge, tdULL is the result of weeks of rapid
experimentation. The resulting code is not especially legible or clean. It
would not normally be published yet, if that were not a requirement for the
PACE Challenge.
