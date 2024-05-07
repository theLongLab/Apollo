
![2](https://github.com/theLongLab/Apollo/assets/55466094/51519eae-aacc-4afb-87aa-8a6ce20656ef)
# Apollo

Apollo is a comprehensive simulation software powered by **CATE** our large scale parallel processing architecture [MEE, CATE: A fast and scalable CUDA implementation to conduct highly parallelized evolutionary tests on large scale genomic data](https://doi.org/10.1111/2041-210X.14168)

---
#### Description

Apollo is a stochastic evolutionary simulator for haploid viruses, designed to study the progression of a disease across a population. Apollo’s simulations move forward in time, and it aims to facilitate the design of highly dynamic and robust simulation models for studying within-host viral evolution and disease spread. 

---
#### Prerequisites

1. CUDA capable hardware
2. LINUX or UNIX based kernel
3. NVIDIA's CUDA toolkit (nvcc compiler)
4. C++ compiler (gcc compiler)

---

#### How to INSTALL

To install Apollo you may have to compile the code using an nvcc compiler. If so execute the following on the terminal:

Download the repository:
````
git clone "https://github.com/theLongLab/Apollo/"
````
````
cd CATE/
````
*cuda 11.3.0 or higher*
````
module load cuda/11.3.0
````

Finally, compile the project:
````
nvcc -std=c++17 *.cu *.cpp -o "Apollo"
````
---
