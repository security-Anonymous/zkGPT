# zkGPT

## Introduction

This is the implementation of zkGPT, which is a SNARK for LLM inference. 
This project partially references code from the work of zkCNN and Lasso.
**Current implementation only supports GPT-2.** 


## Requirement
### Software Requirement
- C++14
- cmake >= 3.10
- GMP library

### Recommended Server Configuration
To ensure smooth execution, we recommend using a server with the following specifications:
- Operating System (OS): Linux (e.g., Ubuntu 18.04 or later)
- Processor (CPU): Multi-core CPU, preferably with 16 cores or more
- Memory (RAM): at least 200GB


## Experiment Script
### Clone the repo
To run the code, make sure you clone with
``` bash
git clone --recurse-submodules git@github.com:security-Anonymous/zkGPT.git
```
since the mcl library is included as a submodule.

### Install GMP Library
```
sudo apt install libgmp-dev
```

### Run a demo of proving LLM inference
The script to run LLM inference proving.
``` bash
./llm.sh
```
### Performance Notes
The circuit initiation phase may take minutes because brute-forcely computing plaintext matrix multiplication using CPU single thread. We will optimize this part in future. Note this part will not become bottleneck in real world, due to it can be computed extremely efficiently on GPU.

## Reference
- [zkCNN: Zero knowledge proofs for convolutional neural network predictions and accuracy](https://doi.org/10.1145/3460120.3485379).

- [Unlocking the lookup singularity with Lasso](https://eprint.iacr.org/2023/1216)

- [mcl](https://github.com/herumi/mcl)

