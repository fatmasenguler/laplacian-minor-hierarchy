# Laplacian Minor Hierarchy for Many-Body Protein Communication

Code accompanying:

> Senguler Ciftci F. and Erman B. *The Geometry of Allostery: A Laplacian Minor Hierarchy for Many-Body Protein Communication.* (2026)

This repository provides code for computing higher-order effective-distance invariants of a protein contact network directly from PDB structures.

The main runnable file currently provided here is a Google Colab notebook for the third-order cooperation index, ( \chi_{ijk} ).

---

## What is in this repository?

| File                                             | What it does                                                                         | How to run it       |
| ------------------------------------------------ | ------------------------------------------------------------------------------------ | ------------------- |
| [`compute_chi_ijk.ipynb`](compute_chi_ijk.ipynb) | Computes the third-order cooperation index ( \chi_{ijk} ) for two uploaded PDB files | Run in Google Colab |

The notebook compares two protein structures. When you run it, it asks you to upload exactly two `.pdb` files.

---

## Before you start: a note about GitHub buttons

This page is hosted on GitHub. You do **not** need to know GitHub or Git to use the notebook.

However, GitHub has several buttons that can be confusing:

* **Copy** copies text to your clipboard. It does **not** download a runnable file.
* **Raw** shows the plain file contents in the browser. For a notebook, this may look like JSON text.
* **Download raw file** downloads the actual file to your computer.
* **Download ZIP** downloads the whole repository.

If you only want to run the notebook, use **Download raw file** or **Download ZIP**. Do not use the small **Copy** buttons unless you only want to copy a piece of text.

---

## Option 1 — Run the notebook in Google Colab

This is the easiest option for readers who do not use GitHub.

1. Open [`compute_chi_ijk.ipynb`](compute_chi_ijk.ipynb) on GitHub.

2. Click the download icon, usually labeled **Download raw file**, to save the notebook.

3. Make sure the downloaded file is named:

   ```text
   compute_chi_ijk.ipynb
   ```

4. Open [Google Colab](https://colab.research.google.com/).

5. In Colab, choose **File → Upload notebook**.

6. Upload `compute_chi_ijk.ipynb`.

7. Run the cells from top to bottom.

8. When prompted, upload exactly two `.pdb` files.

The notebook will compute ( \chi_{ijk} ), create output CSV files, and prepare downloadable results.

---

## Option 2 — Download the whole repository

Use this option if you want a local copy of everything in the repository.

1. Go to the main repository page.

2. Click the green **`<> Code`** button.

3. Choose **Download ZIP**.

4. Unzip the downloaded file on your computer.

5. Open the file:

   ```text
   compute_chi_ijk.ipynb
   ```

6. Upload that notebook to Google Colab using the steps in Option 1.

---

## Option 3 — Use Git, if you already know Git

If you already use Git, you can clone the repository:

```bash
git clone https://github.com/fatmasenguler/laplacian-minor-hierarchy.git
cd laplacian-minor-hierarchy
```

Then open `compute_chi_ijk.ipynb`.

This option is only necessary if you already use Git. If you do not know Git, use Option 1 or Option 2 instead.

---

## Running the notebook

The notebook is designed for Google Colab.

In Colab:

1. Run the notebook cells from top to bottom.
2. When the upload box appears, upload exactly two PDB files.
3. Wait for the calculations to finish.
4. Download the generated output files.

The notebook produces tables of the third-order cooperation index ( \chi_{ijk} ). Large outputs may be split into multiple CSV files so that they can be opened safely in spreadsheet software.

---

## Required input files

You need two protein structures in PDB format.

Each file should have the extension:

```text
.pdb
```

For the structures used in the manuscript, see the **Structures** section below.

---

## Python packages

If you run the notebook in Google Colab, the required packages are handled inside the notebook environment.

If you adapt the notebook for local Python use, the main Python packages are:

```bash
pip install numpy scipy pandas biopython
```

The notebook also uses Google Colab upload/download tools, so local use may require small changes to the file-upload parts of the code.

---

## The method in brief

The protein is represented as a ( C_\alpha ) contact graph with weighted Laplacian ( L ).

From the Moore-Penrose pseudoinverse

[
K = L^{+},
]

the pairwise effective distance is

[
R_{ab} = K_{aa} + K_{bb} - 2K_{ab}.
]

The overlap of two paths anchored at a reference residue ( i ) is

[
K^{(i)}_{ab}
============

\frac{1}{2}
\left(
R_{ia} + R_{ib} - R_{ab}
\right).
]

The third-order cooperation index is

[
\chi_{ijk}
==========

## 1

\frac{
\left(K^{(i)}*{jk}\right)^2
}{
R*{ij} R_{ik}
}.
]

The manuscript also discusses the fourth-order normalized invariant,

[
\chi_{ijkl}
===========

\frac{
\det G^{(i)}*{jkl}
}{
R*{ij} R_{ik} R_{il}
},
]

where

[
G^{(i)}_{jkl}
=============

\begin{pmatrix}
R_{ij} & K^{(i)}*{jk} & K^{(i)}*{jl} \
K^{(i)}*{jk} & R*{ik} & K^{(i)}*{kl} \
K^{(i)}*{jl} & K^{(i)}*{kl} & R*{il}
\end{pmatrix}.
]

Both indices lie in the interval ([0,1]). See the manuscript Methods and Appendix C for the full derivation.

---

## Parameters used in the paper

To reproduce the manuscript figures, use the same parameter choices described in the Methods section of the paper.

The relevant parameters include:

* Contact cutoff distance
* Edge-weighting rule
* Sequence-neighbor exclusion rule
* Chain selection
* Residue indexing convention

For consistency, use the same settings for all structures being compared.

---

## Structures

PSD-95 PDZ3 variants, all T-2F bound, from the [RCSB PDB](https://www.rcsb.org/):

| PDB ID | Variant       |
| ------ | ------------- |
| 5HED   | Wild type     |
| 5HF1   | G330T         |
| 5HFC   | H372A         |
| 5HFF   | G330T + H372A |

Download the PDB files from the RCSB PDB and upload them to the notebook when prompted.

---

## Troubleshooting

### I clicked Copy, but I did not get a Python or notebook file

The **Copy** button only copies text. It does not download the file.

Use **Download raw file** or **Download ZIP** instead.

### The notebook opens as strange text

A Jupyter notebook file, `.ipynb`, is stored internally as JSON. If you click **Raw**, GitHub may show that JSON text.

This is normal. To use the notebook, download the file and open it in Google Colab.

### Google Colab asks me to upload files

That is expected. Upload exactly two `.pdb` files when prompted.

### The output is split into several files

That is also expected for large proteins. The notebook may split large tables into smaller CSV files so that they can be opened in Excel or other spreadsheet programs.

---

## Citation and license

If you use this code, please cite the paper above.

Released under the MIT License. See [`LICENSE`](LICENSE).
