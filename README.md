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


## Code checklist
Flag merge requests ready, as soon as everything is done, but before

- Double check changed files (see "Changes" in the merge request window)
  - Did files get committed that should not be pushed?
  - Is the CMake clean (no hpps under sources etc.)
- Check clang-tidy, CLion remarks on spelling etc. and styleguide compliance (see above.)
- Code formatting / clang-format: On CLion hit ctrl+alt+l for a nice autoformat on all changed files
- Check compiler warnings, there should be none
- Check for unused imports

## Run benchmarks

```bash
simex e launch --instset test --experiment gurobi-solve
python3 eval.py eval_out
python3 check_correctness.py eval_out
```

```bash
rm -rf builds/* aux/* output/* # hard reset simexpal
```




