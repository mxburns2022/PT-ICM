# PT-ICM

An implementation of Parallel Tempering with Isoenergetic Cluster Moves (ICM). The parallel tempering code implements the (NN)_a / M-H scheme described in Malakis et al. 2013. If the `--preprocess` flag is passed, then the code will use a multihistogram method to estimate the density of states and space the temperatures to ensure a constant acceptance rate (described in Bittner et. al 2008), otherwise it will linearly space between the temperatures given. ICMs were proposed and described in Zhu et al. 2015 "Efficient Cluster Algorithm for Spin Glasses in Any Space Dimension". The implementation here still requires testing for the expected efficiency gains.

## Build & Run

### Dependencies

1. [Install Docker](https://docs.docker.com/install/) and VS Code [Dev containers](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers) extension.

2. Open the workspace, switch to the branch you want to work on, and VS Code will prompt you to build the container. This will take a while the first time, but subsequent builds will be faster. If prompt is not shown, you can build the container manually by running the command `Remote-Containers: Rebuild Container` from the command palette.

Now you have the environment with `mpc++`, `cmake`. It won't mess up your local environment. Also you can install your favorite VS Code extensions in the container.

### Building

1. From PT_ICM directory: `mkdir build;cd build`. Creates a `build` directory and enters into it.
2. `cmake ..`: Builds in "Release" mode by default, will look for Boost libraries for arg parsing
3. `make`: Build the project. It will also build a static library `libreplica.a` for general Metropolis-Hastings replica utility use.
