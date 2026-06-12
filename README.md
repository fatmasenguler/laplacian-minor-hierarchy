# Laplacian Minor Hierarchy for Many-Body Protein Communication

Code accompanying:

> Senguler Ciftci F. and Erman B. *The Geometry of Allostery: A Laplacian Minor Hierarchy for Many-Body Protein Communication.* (2026)

This repository provides Google Colab notebooks for computing higher-order effective-distance invariants of a protein contact network directly from PDB structures.

The notebooks compute:

* the third-order cooperation index ( \chi_{ijk} )
* the fourth-order normalized invariant ( R_{ijkl} )

---

## What is in this repository?

| File                                                   | What it does                                                                                                    | How to run it       |
| ------------------------------------------------------ | --------------------------------------------------------------------------------------------------------------- | ------------------- |
| [`compute_chi_ijk.ipynb`](compute_chi_ijk.ipynb)       | Computes the third-order cooperation index ( \chi_{ijk} ) for two uploaded PDB files                            | Run in Google Colab |
| [`compute_Rijklipynb.ipynb`](compute_Rijklipynb.ipynb) | Computes the fourth-order invariant ( R_{ijkl} ) for a fixed triple ( i,j,k ) and every remaining residue ( l ) | Run in Google Colab |

The first notebook, `compute_chi_ijk.ipynb`, compares two protein structures.

The second notebook, `compute_Rijklipynb.ipynb`, analyzes one protein structure at a time. It asks for a PDB file, a chain ID, three residue numbers ( i,j,k ), and then computes ( R_{ijkl} ) for every valid residue ( l ).

> **Important:** Make sure both notebook files are uploaded to the repository root. If `compute_Rijklipynb.ipynb` is missing, the second link above will not work.

---

## Before you start: a note about GitHub buttons

This page is hosted on GitHub. You do **not** need to know GitHub or Git to use the notebooks.

However, GitHub has several buttons that can be confusing:

* **Copy** copies text to your clipboard. It does **not** download a runnable file.
* **Raw** shows the plain file contents in the browser. For a notebook, this may look like JSON text.
* **Download raw file** downloads the actual file to your computer.
* **Download ZIP** downloads the whole repository.

If you only want to run the notebooks, use **Download raw file**, **Download ZIP**, or open the notebook directly in Google Colab. Do not use the small **Copy** buttons unless you only want to copy a piece of text.

---

## Option 1 — Run a notebook in Google Colab

This is the easiest option for readers who do not use GitHub.

1. Open the notebook you want to use:

   * [`compute_chi_ijk.ipynb`](compute_chi_ijk.ipynb)
   * [`compute_Rijklipynb.ipynb`](compute_Rijklipynb.ipynb)

2. If you see an **Open in Colab** button at the top of the notebook, click it.

3. If you do not see that button:

   * click the download icon, usually labeled **Download raw file**
   * save the notebook to your computer
   * open [Google Colab](https://colab.research.google.com/)
   * choose **File → Upload notebook**
   * upload the `.ipynb` file

4. Run the notebook cells from top to bottom.

---

## Option 2 — Download the whole repository

Use this option if you want a local copy of everything in the repository.

1. Go to the main repository page.

2. Click the green **`<> Code`** button.

3. Choose **Download ZIP**.

4. Unzip the downloaded file on your computer.

5. Open the notebook you want:

   ```text
   compute_chi_ijk.ipynb
   compute_Rijklipynb.ipynb
   ```

6. Upload the notebook to Google Colab using the steps in Option 1.

---

## Option 3 — Use Git, if you already know Git

If you already use Git, you can clone the repository:

```bash
git clone https://github.com/fatmasenguler/laplacian-minor-hierarchy.git
cd laplacian-minor-hierarchy
```

Then open one of the notebooks:

```text
compute_chi_ijk.ipynb
compute_Rijklipynb.ipynb
```

This option is only necessary if you already use Git. If you do not know Git, use Option 1 or Option 2 instead.

---

## Notebook 1: `compute_chi_ijk.ipynb`

Use this notebook to compute the third-order cooperation index ( \chi_{ijk} ).

### Input

The notebook asks you to upload exactly two PDB files.

Each file should have the extension:

```text
.pdb
```

### What it computes

For each uploaded structure, the notebook builds a protein contact network and computes the third-order cooperation index

[
\chi_{ijk}
]

for residue triples.

### Output

The notebook creates output tables, usually in CSV format. For large proteins, the output may be split into multiple files so that the tables can be opened safely in spreadsheet software such as Excel.

---

## Notebook 2: `compute_Rijklipynb.ipynb`

Use this notebook to compute the fourth-order invariant ( R_{ijkl} ).

This notebook analyzes one PDB structure at a time.

### Input

When you run the notebook, it asks for:

1. the PDB file path, for example:

   ```text
   5HED.pdb
   ```

2. the chain ID, for example:

   ```text
   A
   ```

   If you leave the chain ID blank, the notebook uses all chains.

3. residue ( i ), using the PDB residue number

4. residue ( j ), using the PDB residue number

5. residue ( k ), using the PDB residue number

6. the weighting scheme

The available weighting schemes are:

```text
unweighted
inv_d2
exp
```

Their meanings are:

| Scheme       | Edge weight                   |
| ------------ | ----------------------------- |
| `unweighted` | ( w_{ij} = 1 )                |
| `inv_d2`     | ( w_{ij} = 1/d_{ij}^{2} )     |
| `exp`        | ( w_{ij} = \exp(-d_{ij}/kT) ) |

For the exponential scheme, `kT` must be in Angstrom units. If you leave `kT` blank, the notebook uses the mean contact distance as a fallback.

### What it computes

After the user chooses residues ( i,j,k ), the notebook computes

[
R_{ijkl}
]

for every valid residue ( l ) in the protein, excluding ( i,j,k ).

In other words, ( i,j,k ) are fixed, and the notebook scans over all possible fourth residues ( l ).

### Output

The notebook prints a table like this:

```text
i       j       k       l          R_ijkl
--------------------------------------------------
330     327     372     298        0.902967
330     327     372     299        0.897498
...
```

It also saves the results to a text file with a name such as:

```text
5HED_chainA_exp_kT1_Rijkl_i330_j327_k372.txt
```

The output file contains the selected residues, the weighting scheme, diagnostic information, and the computed ( R_{ijkl} ) values.

---

## Required input files

You need protein structures in PDB format.

Each file should have the extension:

```text
.pdb
```

For the structures used in the manuscript, see the **Structures** section below.

---

## Python packages

If you run the notebooks in Google Colab, the required packages are handled inside the notebook environment.

The notebooks use common scientific Python packages, including:

```bash
pip install numpy scipy pandas biopython
```

The fourth-order notebook installs Biopython directly in Colab.

If you adapt the notebooks for local Python use, you may need to modify the file-upload and file-download parts of the code.

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

The fourth-order normalized invariant is

[
R_{ijkl}
========

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

* contact cutoff distance
* edge-weighting rule
* sequence-neighbor exclusion rule
* chain selection
* residue indexing convention

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

Download the PDB files from the RCSB PDB and upload them to the notebooks when prompted.

---

## Troubleshooting

### I clicked Copy, but I did not get a notebook file

The **Copy** button only copies text. It does not download the file.

Use **Download raw file**, **Download ZIP**, or **Open in Colab** instead.

### The notebook opens as strange text

A Jupyter notebook file, `.ipynb`, is stored internally as JSON. If you click **Raw**, GitHub may show that JSON text.

This is normal. To use the notebook, download the file and open it in Google Colab.

### Google Colab asks me to upload files

That is expected.

For `compute_chi_ijk.ipynb`, upload exactly two `.pdb` files.

For `compute_Rijklipynb.ipynb`, make sure the PDB file is available in the Colab session, then enter its file path when prompted.

### The fourth-order notebook says “File not found”

Make sure the PDB file has been uploaded to the same Colab session.

For example, if the uploaded file is named:

```text
5HED.pdb
```

then enter:

```text
5HED.pdb
```

when the notebook asks for the PDB file path.

### The notebook says a residue was not found

Check that you entered the PDB residue number, not the row number in the notebook output.

Also check that you selected the correct chain ID.

### The contact graph is disconnected

This can happen if the cutoff is too small or if the selected chain does not form a connected contact graph.

Try checking the chain ID and the PDB file. If necessary, adjust the cutoff in the notebook code.

---

## Citation and license

If you use this code, please cite the paper above.

Released under the MIT License. See [`LICENSE`](LICENSE).
