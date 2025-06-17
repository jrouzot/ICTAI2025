# overlapping Memory Dumping Problem with interruptions

## Introduction

In the context of space missions, an efficient data transfer strategy is crucial to avoid onboard buffer overflow. In many missions, such as ESA's Rosetta mission, the data produced by each instrument is temporarily stored in a dedicated buffer and then dumped to Earth during downlink windows under a limited bandwidth.
Previous work has considered the overlapping memory dumping problem (oMDP), which consists in assigning transfer priorities to the memory buffer to avoid data loss and minimize the peak in memory usage. The oMDP is NP-hard and has been tackled using heuristic methods that demonstrated their efficiency on real instances.
Here, we consider additional decisions in the memory dumping plans that are implementable in practice: data transfer from each buffer can be interrupted after a given time, once per downlink window, preventing it from dumping data until the next window. The new problem is called oMDPi (oMDP with interruptions).
We propose a hybrid heuristic to solve the general oMDPi, embedding a flow relaxation and the single-window heuristic. The results on both real and realistic generated instances show that our heuristic achieves a significant reduction of memory peaks in a reasonable time compared to previous work, making the new policy attractive for future space missions.

## Setup

- Download Cplex from the official website: [Install cplex](https://www.ibm.com/fr-fr/products/ilog-cplex-optimization-studio/cplex-optimizer)
- Update the path to cplex and concert in CMakeList.txt:

```cmake
# Add paths for CPLEX
# set(MAIN_DIR "/home/jrouzot/")
# ==>
set(MAIN_DIR "/home/your/main/dir")
```

- Make a build folder:

```bash
mkdir build
```

- Open a terminal to build the project:

```bash
cmake -S . -B build 
cmake --build build
```

- Run the example:

```bash
./build/omdpi instances/MTP/MTP011
```

- To reproduce the experiments:

```bash
python runexpe.py instances/MTP results/MTP
# AND/OR
python runexpe.py instances/generated results/generated
```