# Paralelisasi Perhitungan Interaksi Antaratom pada Struktur Protein GroEL 1XCK Menggunakan CUDA

## Informasi Project
Project ini bertujuan untuk mengimplementasikan dan membandingkan perhitungan interaksi antaratom secara sekuensial (CPU) dan terparalelisasi (CUDA dan MPI) pada struktur protein GroEL (PDB ID: 1XCK).

## Dataset
- **PDB ID**: 1XCK
- **Protein**: GroEL
- **Jumlah Record ATOM**: 53.970
- **Jumlah Pasangan Atom Unik**: 1.456.353.465 (dihitung dengan N(N-1)/2)
- **Lokasi File**: `data/1XCK.pdb`

## Tujuan Project
- Mengukur dan membandingkan performa perhitungan interaksi antaratom antara CPU sekuensial, CUDA, dan MPI.
- Mengidentifikasi potensi speedup dan efisiensi paralelisasi pada dataset berukuran besar.

## Struktur Folder
```
protein-cuda-1xck/
├── data/
│   └── 1XCK.pdb
├── src/
│   ├── protein_cpu.c
│   ├── protein_cuda.cu
│   └── protein_mpi.c
├── notebook/
│   └── protein_1xck.ipynb
├── results/
├── docs/
├── README.md
└── .gitignore
```
