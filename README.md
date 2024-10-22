# 🚀 High-Performance Computing (HPC) Course Projects

In this course, I explored parallel programming techniques to optimize performance for computationally heavy tasks. Below are two key assignments where I used different approaches to parallelization in **C**.

## 🖥️ OpenMP Parallelization with SIMD in C
In this assignment, I used **OpenMP** to parallelize a computational task in **C**, combining traditional multithreading with **SIMD (Single Instruction, Multiple Data)** for vectorized operations. This further enhanced the efficiency of the program by utilizing data-level parallelism.

### Key Features:
- **Shared Memory Parallelism**: Implemented parallelism in shared memory using OpenMP directives.
- **SIMD Optimization**: Used SIMD to vectorize operations, allowing the program to process multiple data points simultaneously within a single CPU instruction.
- **Task Decomposition**: Broke down complex operations into smaller, independent tasks for better performance.
- **Performance Gains**: Achieved significant reduction in computation time by parallelizing loops and applying SIMD to handle operations on large data sets more efficiently.

---

## 🔗 C11 Threads Parallelization in C
For this assignment, I leveraged **C11 threads** to manually manage parallel execution, creating and synchronizing threads for an efficient distribution of tasks.

### Key Features:
- **Low-level Thread Management**: Used C11’s thread library to create, join, and manage threads.
- **Concurrency Control**: Implemented synchronization mechanisms (e.g., mutexes) to handle shared data access and prevent race conditions.
- **Fine-grained Control**: Provided greater control over thread execution compared to OpenMP.

---

### Technologies Used:
- **C**
- **OpenMP** with SIMD
- **C11 Threads**

