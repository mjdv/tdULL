# tdULL

tdULL is a submission for the [PACE 2020 Treedepth
Challenge](https://pacechallenge.org/2020/) by Ruben Brokkelkamp, Raymond van
Venetië, Mees de Vries and Jan Westerdiep. It is a tool for computing the exact
[treedepth](https://en.wikipedia.org/wiki/Tree-depth) of a general graph.

## Usage

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
