# PARCO Computing 2026 – Deliverable 2  
## SpMV Strong and Weak Scaling with MPI and Hybrid (MPI+OpenMP)

**Author:** Alessandro Turri  
**Course:** Parallel Computing 2025/2026 – University of Trento  

---

## 1. Overview

This project analyses the scalability of **Sparse Matrix–Vector Multiplication (SpMV)** using the **CSR (Compressed Sparse Row)** format on a distributed-memory HPC system.

The goal of **Deliverable 2** is to study:

- MPI-only parallel execution
- Hybrid MPI + OpenMP execution
- Strong scaling behaviour
- Weak scaling behaviour
- Communication vs computation impact

All experiments were executed on the **UniTN HPC cluster** using PBS.

---

## 2. Repository Structure

repo/  
├── README.md  
├── src/  
│   ├── smpv_d2.cpp  
│   └── smpv_d2  
├── scripts/  
│   ├── strong_mpi.pbs  
│   ├── weak_mpi.pbs  
│   ├── strong_hybrid_128.pbs  
│   └── weak_hybrid_128.pbs  
├── matrix/  
│   ├── xenon2/xenon2.mtx  
│   ├── SiO2/SiO2.mtx  
│   ├── StocF-1465/StocF-1465.mtx  
│   └── cage14/cage14.mtx  
├── results/  
│   ├── strong_<matrix>_mpi.csv  
│   ├── strong_<matrix>_hybrid_128.csv  
│   ├── weak_mpi.csv  
│   └── weak_hybrid_128.csv  
└── plots/  

---

## 3. Matrices

Strong scaling experiments use real matrices from the SuiteSparse collection:

- xenon2  
- SiO2  
- StocF-1465  
- cage14  

Matrices are stored as:

matrix/<name>/<name>.mtx  

Weak scaling experiments use a **synthetic matrix generator** to keep a constant workload per MPI rank.

---

## 4. Parallel Configurations

### MPI-only
- OpenMP disabled at runtime
- OMP_NUM_THREADS=1

### Hybrid MPI + OpenMP
- Fixed total number of cores: 128
- Exact (P,T) pairs such that P × T = 128:

(128,1), (64,2), (32,4), (16,8), (8,16), (4,32), (2,64)

---

## 5. Compilation

Compilation is performed **inside each PBS job** to ensure reproducibility with the cluster toolchain:

mpicxx -O3 -std=c++17 -fopenmp -march=native -o src/smpv_d2 src/smpv_d2.cpp

The same executable is used for both MPI-only and hybrid runs.

---

## 6. Running on the UniTN HPC Cluster

All jobs use the queue `short_cpuQ` and allocate up to 4 nodes:

select=4:ncpus=96:mpiprocs=32  

### Strong scaling (example: xenon2)

MPI-only:
qsub -v MATRIX=matrix/xenon2/xenon2.mtx scripts/strong_mpi.pbs  

Hybrid:
qsub -v MATRIX=matrix/xenon2/xenon2.mtx scripts/strong_hybrid_128.pbs  

### Weak scaling

MPI-only:
qsub scripts/weak_mpi.pbs  

Hybrid:
qsub scripts/weak_hybrid_128.pbs  

---

## 7. Parameters

PBS scripts support parameter overrides:

- ITERS (default 50)  
- WARMUP (default 3)  
- RUNS (default 10)  
- SEED (default 1)  

Example:
qsub -v MATRIX=matrix/cage14/cage14.mtx,ITERS=80,RUNS=5 scripts/strong_mpi.pbs  

---

## 8. Output Files

Results are written in the `results/` directory:

- strong_<matrix>_mpi.csv  
- strong_<matrix>_hybrid_128.csv  
- weak_mpi.csv  
- weak_hybrid_128.csv  

PBS logs:

- strong_mpi.log  
- weak_mpi.log  
- strong_hybrid_128.log  
- weak_hybrid_128.log  

---

## 9. One-shot Reproducibility Script

The following script clones the repository, downloads all matrices, and submits all PBS jobs.

Create the script:

cat > run_all_d2.sh << EOF
#!/bin/bash
set -euo pipefail

REPO_URL="https://github.com/aleturri28/PARCO-Computing-2026--244927-.git"
MAXJOBS=28

git clone \$REPO_URL repo
cd repo

mkdir -p matrix
cd matrix

wget https://sparse.tamu.edu/MM/Boeing/xenon2.tar.gz
tar -xzf xenon2.tar.gz
mkdir -p xenon2
find . -name "xenon2.mtx" -exec mv {} xenon2/ \;
rm -rf xenon2.tar.gz

wget https://sparse.tamu.edu/MM/Simon/SiO2.tar.gz
tar -xzf SiO2.tar.gz
mkdir -p SiO2
find . -name "SiO2.mtx" -exec mv {} SiO2/ \;
rm -rf SiO2.tar.gz

wget https://sparse.tamu.edu/MM/Chen/StocF-1465.tar.gz
tar -xzf StocF-1465.tar.gz
mkdir -p StocF-1465
find . -name "StocF-1465.mtx" -exec mv {} StocF-1465/ \;
rm -rf StocF-1465.tar.gz

wget https://sparse.tamu.edu/MM/vanHeukelum/cage14.tar.gz
tar -xzf cage14.tar.gz
mkdir -p cage14
find . -name "cage14.mtx" -exec mv {} cage14/ \;
rm -rf cage14.tar.gz

cd ../scripts

qsub weak_mpi.pbs
qsub weak_hybrid_128.pbs

for M in xenon2 SiO2 StocF-1465 cage14; do
  while [ "\$(qstat -u \$USER | grep -c \$USER)" -ge "\$MAXJOBS" ]; do
    sleep 5
  done
  qsub -v MATRIX=matrix/\$M/\$M.mtx strong_mpi.pbs
  qsub -v MATRIX=matrix/\$M/\$M.mtx strong_hybrid_128.pbs
done
EOF

chmod +x run_all_d2.sh  

Run everything:
./run_all_d2.sh  

---

## 10. Notes

- SpMV is memory-bound
- MPI-only scaling is limited by communication
- Hybrid MPI + OpenMP reduces communication overhead
- Weak scaling behaves better than strong scaling
