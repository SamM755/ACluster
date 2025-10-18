# Cluster-based Lossless Encoding Method (C++ Implementation)

This repository contains the C++ implementation of the core algorithms presented in our paper. This version is designed for high performance and serves as the basis for the experimental results reported.

The methods have been designed for integration with data systems like Apache TsFile and Apache IoTDB, ensuring their robustness and scalability.

This paper introduces a novel cluster-based lossless compression framework for multi-dimensional numerical data. Our approach pivots from traditional local, single-dimensional references to using **global, multi-dimensional reference points** identified through our proposed algorithms. By capturing the underlying spatial structure of the entire dataset, our methods, **KCluster** and **ACluster**, achieve superior compression ratios while maintaining high time efficiency.

### Code Structure

For readers interested in the specific implementation details, the core logic of our proposed methods is organized as follows:

-   **`encoding_benchmark/algorithms/`**: This directory contains the implementation of all encoding methods discussed in the paper.
    -   **`kcluster.cpp`** / **`kcluster.h`**: The main implementation of our **KCluster** algorithm.
    -   **`acluster.cpp`** / **`acluster.h`**: The main implementation of our **ACluster** algorithm.
    -   **`cluster_encoder_logic.cpp`** / **`cluster_encoder_logic.h`**: Contains the shared logic for both KCluster and ACluster, including residual calculation and data serialization/deserialization.
    -   **`cluster_common.h`**: Defines the common data structures and constants used across our cluster-based methods.

-   **`encoding_benchmark/main.cpp`**: The main entry point of the benchmark program, which orchestrates the experiments and performance measurements.

## Artifact Evaluation: Building and Running Experiments

This repository provides a standalone C++ environment for building the project and reproducing the key performance results presented in our paper.

### Environment & Dependencies

To build and run the code in this repository, you will need the following environment:

- **C++ Compiler**: A C++17 compliant compiler is required. The experiments in our paper were conducted using:
  - **MinGW-w64 GCC version 15.2.0** on Windows 11.
  - Other modern compilers like `g++` on Linux or `Clang` on macOS should also be compatible.

- **Build System (Optional but Recommended)**:
  - **CMake (version 3.10 or higher)**: The project includes a `CMakeLists.txt` file for easy cross-platform building.
  - **Make** (or an equivalent build tool like Ninja).

### How to Build

This project is configured for building with CMake. The following instructions are tailored for a **Windows environment using MinGW-w64 and CMake**. Users on other platforms (e.g., Linux or macOS) may need to adjust the generator (`-G`) flag or use the default build system.

1.  **Clone the repository:**
    ```shell
    git clone https://github.com/SamM755/ACluster.git
    cd ACluster/encoding_benchmark
    ```

2.  **Create a build directory:**
    ```shell
    mkdir build
    cd build
    ```

3.  **Run CMake and build the project:**
    ```shell
    # Generate build files
    cmake -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release ..
    
    # Compile the source code
    cmake --build .
    ```
    On Windows with MinGW, you might use `mingw32-make` instead of `make`.

    After a successful build, an executable (e.g., `benchmark`) will be created in the `build` directory.

### How to Run

The main entry point for testing is the `benchmark` executable. You can run it from the `build` directory:

```shell
# From within the 'build' directory
./benchmark
```

This will execute the pre-configured performance tests on the datasets included in the repository. The output will display the compression ratio, encoding time, and decoding time for each method, mirroring the results presented in our paper.

### Reproducing Experimental Results

The main source file (e.g., `main.cpp` or `benchmark.cpp`) is pre-configured with the experimental setups described in our paper. By examining this file, you can see how different methods (KCluster, ACluster) and datasets are tested. To reproduce specific results, you can modify the main function to run only the desired tests, following the parameters detailed in the paper's **"Experimental Setup"** section.


### System Deployment in IoTDB and TsFile

The core algorithms presented in our paper are not just theoretical; they have been fully implemented and integrated into Apache TsFile and Apache IoTDB projects. 
This ensures their robustness, scalability, and availability to the wider data management community.

-   **Apache TsFile Implementation:** 
    -   **Link:** [https://github.com/apache/tsfile/tree/research/cluster-compress](https://github.com/apache/tsfile/tree/research/cluster-compress) 

-   **Apache IoTDB Implementation:** 
    -   **Link:** [https://github.com/apache/iotdb/tree/research/cluster-compress](https://github.com/apache/iotdb/tree/research/cluster-compress) 