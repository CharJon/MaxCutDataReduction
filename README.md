# Max Cut Preprocessing

Code for our paper "Separator based Data Reduction for the Maximum Cut Problem", which will be presented at Symposium on Experimental Algorithms (SEA) 2024.

## External dependencies

All dependencies of this project are manged by git (via submodules).

## Building

1) clone project
2) run ```git submodule update --init --recursive```
3) now cmake building should work:

```bash
mkdir build && cd build
cmake .. -DGUROBI=ON -DGUROBI_DIR="/abs/path/to/gurobi/gurobi1001/linux64/"
make -j 10 # -j 10 makes make use 10 threads
```
## Run

For larger graphs, set callstack size on linux to a large value via set ulimit:
```bash
ulimit -s 32000
```
