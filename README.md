*Project forked from [here](https://github.com/gerstung-lab/clonex). Contributions include adding the possibility of catastrophic events leading to decreases in population size.*

# clonex
This archive contains a C programme to calculate a Wright-Fisher model of cancer development.

## Background
This tools models a population growing exponentially from size `n` to `N` in `g` generations.
It allows to specify `d` driver sites with selective advantage `s`; `p` passenger sites, which are neutral and
`o` other sites with selective advantage `t`. Sites mutate at rates `u` (drivers, other sites) and `v`, respectively.

## Installation
As easy as
```{bash}
$ make
```

## Running
As easy as
```{bash}
$ ./clonex -h
usage: clonex [-N:n:u:v:s:t:g:R:f:r:p:d:o:wh]
  N - Maximal population size (default = 1000000000)
  n - Initial population size (default = 1)
  a - Initial growth rate (default = 2)
  u - Mutation rate (default = 1e-07)
  v - Mutation rate passengers (default = u)
  s - Selective advantage (default = 0.01)
  t - Selective advantage of other drivers (default = 0.015)
  g - Number of generations (default = 1800)
  X - generation at which population starts decreasing (default = 0)
  Y - generation at which population stops decreasing (default = 0)
  Z - rate at which population decreases (default = 1e-05)
  L - Population increases after decreasing (default = 0 (no))
  d - Number of drivers (default = 1000)
  p - Number of passengers (default = 0)
  o - Number of other drivers (default = 0)
  R - Replicates (default = 1)
  r - Random seed (default = time)
  f - File directory (Required! Make sure that the directory exists!)
  G - Output every G generations (default = g)
  h - This help

$ ./clonex -N 1000000 -n 1000000 -s 0.01 -d 1000 -u 1e-9 -f foo -p 1000000 -v 1e-8 -g 2000
```
