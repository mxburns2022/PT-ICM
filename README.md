# PT-ICM

[![CMake](https://github.com/mxburns2022/PT-ICM/actions/workflows/build.yml/badge.svg)](https://github.com/mxburns2022/PT-ICM/actions/workflows/build.yml)

An implementation of Parallel Tempering with Isoenergetic Cluster Moves (ICM). The parallel tempering code implements the (NN)_a / M-H scheme described in Malakis et al. 2013. If the `--preprocess` flag is passed, then the code will use a multihistogram method to estimate the density of states and space the temperatures to ensure a constant acceptance rate (described in Bittner et. al 2008), otherwise it will linearly space between the temperatures given. ICMs were proposed and described in Zhu et al. 2015 "Efficient Cluster Algorithm for Spin Glasses in Any Space Dimension". The implementation here still requires testing for the expected efficiency gains.

## Build

### Dependencies

#### Method 1: Docker

1. [Install Docker](https://docs.docker.com/install/) and VS Code [Dev containers](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers) extension.

2. Open the workspace, switch to the branch you want to work on, and VS Code will prompt you to build the container. This will take a while the first time, but subsequent builds will be faster. If prompt is not shown, you can build the container manually by running the command `Remote-Containers: Rebuild Container` from the command palette.

Now you have the environment with `mpc++`, `cmake`. It won't mess up your local environment. Also, you can install your favorite VS Code extensions in the container.

#### Method 2: Local / Bluehive

Load neccessary modules.

### Build the project

Generate necessary files and build the project. At the root of the project, run:

```bash
cmake -B build/
cmake --build build/ -j8
```

It will also build a static library `libreplica.a` for general Metropolis-Hastings replica utility use.

## Run

### Parameters

- `-g` Path to the input graph file
- `--preprocess` the temperatures will be linearly spaced
- `--inverse` enable linearly spacing the inverse temperatures instead of the temperatures themselves
- `-T0` and `-T1` the lowest and highest temperatures to use. If `--inverse` is passed, they should be converted to the inverse temperatures. (2.5 -> 0.4)
- `-pts` parallel tempering steps, 20,000 - 40,000 is enough for G022 (200 nodes)
- `--icm` disable or enable ICM

### Example 1 - Run `pt_icm` locally

```bash
./build/pt_icm -p example_graphs/G022 -T0 0.5 -T1 15 --pts 20000 --icm
```

### Example 2 - Sweep parameters on Bluehive

Modify the parameters in `scripts/run.py`. Then just run the script.