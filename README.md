# Cluster-based Lossless Compression

The core algorithms presented in our paper are not just theoretical; they have been fully implemented and integrated into Apache TsFile and Apache IoTDB projects. 
This ensures their robustness, scalability, and availability to the wider data management community.

-   **Apache TsFile Implementation:** 
    -   **Link:** [https://github.com/apache/tsfile/tree/research/cluster-compress](https://github.com/apache/tsfile/tree/research/cluster-compress) 

-   **Apache IoTDB Implementation:** 
    -   **Link:** [https://github.com/apache/iotdb/tree/research/cluster-compress](https://github.com/apache/iotdb/tree/research/cluster-compress) 


This paper introduces a novel cluster-based lossless compression framework designed for multi-dimensional numerical data. 
Our approach pivots from traditional local, single-dimensional references to using **global, multi-dimensional reference points** identified through advanced cluster-based algorithms. 
By capturing the underlying spatial structure of the entire dataset, our methods, KCluster and ACluster, achieve superior compression ratios while maintaining high time efficiency.


## Artifact Evaluation: Quick Tests and Exp Reproducibility

This repository provides a lightweight, standalone environment for quickly evaluating the compression performance and Approximate Nearest Neighbor (ANN) search accuracy of our proposed methods.

### Environment Setup

To run the code in this repository, you will need the following environment:

-   **Java Development Kit (JDK):** Version **21 or higher** is required, for **Java Vector API** for SIMD (Single Instruction, Multiple Data) acceleration to significantly speed up distance calculations, which is a stable feature in Java 21.

-   **VM Options:** You must configure your IDE or command line to enable the Vector API module.
    -   Add the following to your VM options:
        ```
        --add-modules jdk.incubator.vector
        ```

    -   **In IntelliJ IDEA:**
        1.  Go to `Run` -> `Edit Configurations...`.
        2.  Find the run configuration for `QuickTest`.
        3.  In the `VM options` field, paste `--add-modules jdk.incubator.vector`.
        4.  Click `Apply` and `OK`.

    -   **From the command line:**
        When you run the Java application, include the flag:
        ```shell
        java --add-modules jdk.incubator.vector -cp <your_classpath> your.main.Class
        ```

### How to Run

The main entry point for testing is the `QuickTest.java` file, located in the `src/test/` directory. 
This class contains example method calls for both compression and ANN search functionalities, allowing you to quickly see the algorithms in action.

### Reproducing Experimental Results

The `QuickTest` class is pre-configured with examples that mirror the experimental setups described in our paper. By following the parameters and dataset configurations detailed in the paper's **"Experimental Settings"** section, you can use the code in this repository to reproduce the key results presented.