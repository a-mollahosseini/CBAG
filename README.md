# CBAG: state-of-the-art algorithm for Graph Burning
Official implementation of Centrality BAsed Genetic-algorithm ([CBAG](https://arxiv.org/abs/2208.01008)): state-of-the-art heuristic for solving the graph burning problem.

## How to Use
1. Install the [Boost Library](https://www.boost.org/).
2. Place your your dataset in `./Data` directory. Refer to below section for more information.
3. Change the value of `Path_To_Your_Boost_Here` in `Makefile` file to the directory where your Boost library is installed in.
4. Place an approximation algorithm for betweenness centrality in the `./Code/BetweennessApprox` directory. You may use [1](https://github.com/ecrc/BeBeCA/tree/master/Source_Code/src) or [2](https://networkit.github.io/). Refer to `Solver.hpp` and `Solver.cpp` to see the usage format.
5. Assuming your file is named `graph.csv` and you wish to find a solution with burning length `10`, run the following commands in your terminal: 
```bash
$ make
$ ./CBAG graph.csv 10
```
## Dataset
In order to use the algorithm for any desired graph:
1. Prepare the graph in a file with the following format.
    - First line contains the name of the graph with no spaces. 
    - Second line contains four numbers: the number of nodes, the number of edges, whether the nodes are `1-indexed` or not (`0` if `0-indexed` and `1` if `1-indexed`), and whether the file contains edge weights or not (`0` if not, `1` if yes) respectively. 
    - The next edge number of lines each contains two or three numbers: each of the two nodes that are connected by that edge, and the weight of the edge optionally.
    
    For example, all of the three files below are denoting the K4 graph.
    ```
    Square-1
    4 6 0 0
    0 1
    0 2
    0 3
    1 2
    1 3
    2 3
    ```
    or
    ```
    Square-2
    4 6 1 0
    1 2
    1 3
    1 4
    2 3
    2 4
    3 4
    ```
    or even
    ```
    Square-3
    4 6 1 1
    1 2 23.4
    1 3 12.1
    1 4 24.4
    2 3 -34.3
    2 4 -12.3
    3 4 56.7
    ```
2. Place your file in the `./Data` folder and run the algorithm. 
## Updates
4th June 2024: dataset and external libraries are replaced with links to them.
## License
This project is licensed under the terms of the MIT license.
