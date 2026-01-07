# PARCO Computing 2026 – Deliverable 2  
## SpMV Strong and Weak Scaling with MPI and Hybrid (MPI+OpenMP)

**Author:** Alessandro Turri  
**Course:** Parallel Computing 2025/2026 – University of Trento  

---

## 1. Overview

This project evaluates the scalability of **Sparse Matrix–Vector Multiplication (SpMV)** using the **CSR (Compressed Sparse Row)** format on a distributed-memory HPC system.

The objective of **Deliverable 2** is to analyse:

- MPI-only parallel execution
- Hybrid MPI + OpenMP execution
- Strong scaling behaviour
- Weak scaling behaviour
- Communication vs computation trade-offs

All experiments were executed on the **UniTN HPC cluster** using PBS job scripts.

---

## 2. Repository Structure

'''
repo/
├── README.md
├── .gitignore
├── src/
│   ├── smpv_d2.cpp
│   └── smpv_d2
├── scripts/
│   ├── strong_mpi.pbs
│   ├── weak_mpi.pbs
│   ├── strong_hybrid_128.pbs
│   └── weak_hybrid_128.pbs
├── matrix/
│   └── default.txt
├── results/
│   └── default.txt
├── plots/
│   └── default.txt
└── .git/
'''

Directories `matrix/`, `results/`, and `plots/` contain placeholder files and are populated during experiments.

---

## 3. Matrices

### Strong Scaling

Strong scaling experiments use real sparse matrices from the **SuiteSparse Matrix Collection**, provided externally and placed in the following layout:

'''
matrix/<name>/<name>.mtx
'''

Matrices used in this project:

- xenon2 (Ronis)
- SiO2 (PARSEC)
- cage12 (vanHeukelum)
- stocF-1469 (GHS_indef)

### Weak Scaling

Weak scaling experiments do **not** use real matrices.  
Matrices are generated at runtime using a synthetic generator, keeping a constant workload per MPI rank.

---

## 4. Parallel Configurations

### MPI-only

- One MPI process per core
- OpenMP disabled at runtime
- `OMP_NUM_THREADS=1`
- Explicit `--no-omp` flag passed to the executable

### Hybrid MPI + OpenMP

- Fixed total number of cores: **128**
- Exact `(P, T)` pairs such that:

'''
P × T = 128
'''

Pairs used:

'''
(128,1), (64,2), (32,4), (16,8), (8,16), (4,32), (2,64)
'''

This keeps the total computational resources constant while varying the MPI/OpenMP balance.

---

## 5. Compilation

Compilation is performed **inside each PBS job**, ensuring reproducibility with the cluster toolchain.

The following command is executed within the job scripts:

'''
mpicxx -O3 -std=c++17 -fopenmp -march=native -o src/smpv_d2 src/smpv_d2.cpp
'''

The same executable is used for both MPI-only and hybrid runs.

---

## 6. Execution on the UniTN HPC Cluster

All jobs are submitted to the `short_cpuQ` queue and allocate up to four nodes:

'''
select=4:ncpus=96:mpiprocs=32
'''

### Strong Scaling

MPI-only:
'''
qsub -v MATRIX=matrix/<name>/<name>.mtx scripts/strong_mpi.pbs
'''

Hybrid:
'''
qsub -v MATRIX=matrix/<name>/<name>.mtx scripts/strong_hybrid_128.pbs
'''

### Weak Scaling

MPI-only:
'''
qsub scripts/weak_mpi.pbs
'''

Hybrid:
'''
qsub scripts/weak_hybrid_128.pbs
'''

---

## 7. Parameters

PBS scripts define default values for all parameters, which can be overridden via `qsub -v`.

Main parameters:

- `ITERS`   (default: 50)
- `WARMUP`  (default: 3)
- `RUNS`    (default: 10)
- `SEED`    (default: 1)

Example:

'''
qsub -v MATRIX=matrix/cage12/cage12.mtx,ITERS=80,RUNS=5 scripts/strong_mpi.pbs
'''

---

## 8. Output Files

All output files are written to the `results/` directory.

### Strong Scaling

- `strong_<matrix>_mpi.csv`
- `strong_<matrix>_hybrid_128.csv`

### Weak Scaling

- `weak_mpi.csv`
- `weak_hybrid_128.csv`

PBS standard output logs:

- `strong_mpi.log`
- `weak_mpi.log`
- `strong_hybrid_128.log`
- `weak_hybrid_128.log`

---

## 9. Notes

- SpMV is memory-bound
- MPI-only scalability is limited by communication overhead
- Hybrid MPI + OpenMP reduces communication volume
- Weak scaling shows better behaviour than strong scaling

