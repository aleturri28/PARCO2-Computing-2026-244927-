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

repo/
├── README.md
├── .gitignore
├── src/
│ ├── smpv_d2.cpp
│ └── smpv_d2
├── scripts/
│ ├── strong_mpi.pbs
│ ├── weak_mpi.pbs
│ ├── strong_hybrid_128.pbs
│ └── weak_hybrid_128.pbs
├── matrix/
│ └── default.txt
├── results/
│ └── default.txt
├── plots/
│ └── default.txt
└── .git/

yaml
Copia codice

Directories `matrix/`, `results/`, and `plots/` contain placeholder files and are populated during experiments.

---

## 3. Matrices

### Strong Scaling

Strong scaling experiments use real sparse matrices from the **SuiteSparse Matrix Collection**, downloaded externally and placed in the following layout:

matrix/<name>/<name>.mtx

yaml
Copia codice

Matrices used:

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
- `--no-omp` flag passed to the executable  

### Hybrid MPI + OpenMP

- Fixed total number of cores: **128**
- Exact `(P, T)` pairs such that:

P × T = 128

yaml
Copia codice

Pairs used:

(128,1), (64,2), (32,4), (16,8), (8,16), (4,32), (2,64)

yaml
Copia codice

This keeps the total computational resources constant while varying the MPI/OpenMP balance.

---

## 5. Compilation

Compilation is performed **inside each PBS job**, ensuring reproducibility with the cluster toolchain.

The following command is executed within the job scripts:

mpicxx -O3 -std=c++17 -fopenmp -march=native -o src/smpv_d2 src/smpv_d2.cpp

yaml
Copia codice

The same executable is used for both MPI-only and hybrid runs.

---

## 6. Execution on the UniTN HPC Cluster

All jobs are submitted to the `short_cpuQ` queue and allocate up to four nodes:

select=4:ncpus=96:mpiprocs=32

makefile
Copia codice

### Strong Scaling

MPI-only:
qsub -v MATRIX=matrix/<name>/<name>.mtx scripts/strong_mpi.pbs

makefile
Copia codice

Hybrid:
qsub -v MATRIX=matrix/<name>/<name>.mtx scripts/strong_hybrid_128.pbs

makefile
Copia codice

### Weak Scaling

MPI-only:
qsub scripts/weak_mpi.pbs

makefile
Copia codice

Hybrid:
qsub scripts/weak_hybrid_128.pbs

yaml
Copia codice

---

## 7. Parameters

PBS scripts define default values for all parameters, which can be overridden via `qsub -v`.

Main parameters:

- `ITERS`   (default: 50)
- `WARMUP`  (default: 3)
- `RUNS`    (default: 10)
- `SEED`    (default: 1)

Example:

qsub -v MATRIX=matrix/cage12/cage12.mtx,ITERS=80,RUNS=5 scripts/strong_mpi.pbs

yaml
Copia codice

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

## 9. One-shot Reproducibility (copy & paste)

The following commands clone the repository, download all required matrices, and submit all PBS jobs.  
They are intended to be executed on the UniTN HPC cluster.

### 9.1 Clone the repository

git clone https://github.com/aleturri28/PARCO2-Computing-2026-244927-.git
cd PARCO2-Computing-2026-244927-

graphql
Copia codice

### 9.2 Download and prepare matrices

mkdir -p matrix
cd matrix

markdown
Copia codice

**xenon2 (Ronis)**

wget https://sparse.tamu.edu/MM/Ronis/xenon2.tar.gz
tar -xzf xenon2.tar.gz
mkdir -p xenon2
find . -name "xenon2.mtx" -exec mv {} xenon2/ ;
rm -rf xenon2.tar.gz

markdown
Copia codice

**SiO2 (PARSEC)**

wget https://sparse.tamu.edu/MM/PARSEC/SiO2.tar.gz
tar -xzf SiO2.tar.gz
mkdir -p SiO2
find . -name "SiO2.mtx" -exec mv {} SiO2/ ;
rm -rf SiO2.tar.gz

markdown
Copia codice

**cage12 (vanHeukelum)**

wget https://sparse.tamu.edu/MM/vanHeukelum/cage12.tar.gz
tar -xzf cage12.tar.gz
mkdir -p cage12
find . -name "cage12.mtx" -exec mv {} cage12/ ;
rm -rf cage12.tar.gz

scss
Copia codice

**stocF-1469 (GHS_indef)**

wget https://sparse.tamu.edu/MM/GHS_indef/stocF-1469.tar.gz
tar -xzf stocF-1469.tar.gz
mkdir -p stocF-1469
find . -name "stocF-1469.mtx" -exec mv {} stocF-1469/ ;
rm -rf stocF-1469.tar.gz

vbnet
Copia codice

Return to the repository root:

cd ..

shell
Copia codice

### 9.3 Submit PBS jobs

Weak scaling (submitted once):

qsub scripts/weak_mpi.pbs
qsub scripts/weak_hybrid_128.pbs

css
Copia codice

Strong scaling (for each matrix):

qsub -v MATRIX=matrix/xenon2/xenon2.mtx scripts/strong_mpi.pbs
qsub -v MATRIX=matrix/xenon2/xenon2.mtx scripts/strong_hybrid_128.pbs

qsub -v MATRIX=matrix/SiO2/SiO2.mtx scripts/strong_mpi.pbs
qsub -v MATRIX=matrix/SiO2/SiO2.mtx scripts/strong_hybrid_128.pbs

qsub -v MATRIX=matrix/cage12/cage12.mtx scripts/strong_mpi.pbs
qsub -v MATRIX=matrix/cage12/cage12.mtx scripts/strong_hybrid_128.pbs

qsub -v MATRIX=matrix/stocF-1469/stocF-1469.mtx scripts/strong_mpi.pbs
qsub -v MATRIX=matrix/stocF-1469/stocF-1469.mtx scripts/strong_hybrid_128.pbs

yaml
Copia codice

Job monitoring:

qstat -u $USER

yaml
Copia codice

---

## 10. Notes

- SpMV is memory-bound  
- MPI-only scalability is limited by communication overhead  
- Hybrid MPI + OpenMP reduces communication volume  
- Weak scaling shows better behaviour than strong scaling  
