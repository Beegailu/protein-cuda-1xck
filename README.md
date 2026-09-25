# Paralelisasi Analisis Contact Map Struktur Protein GroEL 1XCK Menggunakan CUDA

Project ini melakukan analisis contact map pada struktur protein GroEL
(PDB ID 1XCK) menggunakan CUDA untuk memparalelkan perhitungan jarak
antarbagian struktur protein.

## Dataset

- Protein: GroEL
- PDB ID: 1XCK
- File: `data/1XCK.pdb`
- Jumlah atom: 53.970
- Jumlah pasangan atom unik: 1.456.353.465

## Struktur Project

- `data/` : dataset protein
- `src/` : source code CPU dan CUDA
- `notebook/` : notebook eksperimen
- `results/contact_map/` : hasil contact map
- `results/benchmark/` : hasil pengujian performa
- `docs/` : dokumentasi project

## Status

Project masih dalam tahap pengembangan.