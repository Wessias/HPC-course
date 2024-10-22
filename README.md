# 🚀 High-Performance Computing (HPC) Course Projects

In this course, I explored parallel programming techniques to optimize performance for computationally heavy tasks. Below are two key assignments where I used different approaches to parallelization in **C**.

## 🔗 C11 Threads Parallelization in C

For this assignment, I used **C11 threads** to parallelize the implementation of **Newton's method** for finding the roots of the equation \( z^d + 1 = 0 \). The solution involved calculating which root each pixel converges to, and then generating images based on the results. The images are colored in two ways:
1. **By the number of iterations** required for convergence.
2. **By the root** each pixel converges to.

The example images below are generated for \( z^7 + 1 \).

### Key Features:
- **Newton's Method**: Implemented to find the complex roots of \( z^d + 1 \), with a focus on visualizing the iterative process.
- **Multithreading**: Used C11 threads to parallelize the computation of convergence for each pixel across multiple cores.
- **Image Generation**: Created two types of images:
  1. An image where the colors represent the **number of iterations** needed for each pixel to converge to a root.
  2. An image where the colors represent **which root** each pixel converges to.
- **Concurrency Control**: Managed thread synchronization and shared resources to ensure smooth parallel computation.

### Image Examples for \( z^7 + 1 \)
<div align="center">
  <img src="https://i.imgur.com/DU6y0gn.png" alt="Iterations to Convergence" width="400" />
  <img src="https://i.imgur.com/FH34I9O.png" alt="Root Convergence" width="400" />
</div>

---

## 🖥️ OpenMP Parallelization with SIMD in C
In this assignment, I used **OpenMP** to parallelize a computational task in **C**, combining traditional multithreading with **SIMD (Single Instruction, Multiple Data)** for vectorized operations. This further enhanced the efficiency of the program by utilizing data-level parallelism.

### Key Features:
- **Shared Memory Parallelism**: Implemented parallelism in shared memory using OpenMP directives.
- **SIMD Optimization**: Used SIMD to vectorize operations, allowing the program to process multiple data points simultaneously within a single CPU instruction.
- **Task Decomposition**: Broke down complex operations into smaller, independent tasks for better performance.
- **Performance Gains**: Achieved significant reduction in computation time by parallelizing loops and applying SIMD to handle operations on large data sets more efficiently.

---

### Technologies Used:
- **C**
- **C11 Threads**
- **OpenMP** with SIMD
