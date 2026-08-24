# Laplacian Minor Hierarchy for Many-Body Protein Communication

Code accompanying:

> Senguler Ciftci F. and Erman B. *The Geometry of Allostery: A Laplacian Minor Hierarchy for Many-Body Protein Communication.* (2026)

This repository computes higher-order effective-distance invariants of a protein contact network directly from PDB structures.

The hierarchy has three levels:

* the second-order effective resistance $R_{ij}$
* the third-order cooperation index $\chi_{ijk}$
* the fourth-order normalized invariant $R_{ijkl}$

**Sign convention.** Every index in this repository runs from 0 to 1, with

```text
0  uncorrelated  (independent communication pathways)
1  correlated    (fully overlapping pathways)
```

The same convention is used in all notebooks, all scripts, and the manuscript.

---

## What is in this repository?

### Colab notebooks

| File | What it does |
| ---- | ------------ |
| [`compute_chi_ijk.ipynb`](compute_chi_ijk.ipynb) | Third-order index $\chi_{ijk}$ for all ordered triples in each of two uploaded PDB files |
| [`compute_Rijkl.ipynb`](compute_Rijkl.ipynb) | Fourth-order invariant $R_{ijkl}$ for a fixed triple $i,j,k$ and every remaining residue $l$ |
| [`AlphaFold3_5HED_25_SECONDS_ensemble_UPLOAD_FIRST.ipynb`](AlphaFold3_5HED_25_SECONDS_ensemble_UPLOAD_FIRST.ipynb) | Runs AlphaFold 3 on 5HED, then applies the second-order analysis across a 25-model ensemble |

### Python modules

| File | What it does |
| ---- | ------------ |
| [`SECONDS.py`](SECONDS.py) | Second invariant. Compares the $\sqrt{\langle R_{ij}\rangle}$ residue profile of two structures |
| [`THIRDS.py`](THIRDS.py) | Third invariant. Computes $\chi^{(i)}_{j,k}$ for a fixed anchor $i$ and partner $j$, scanning $k$ |
| [`FOURTHS.py`](FOURTHS.py) | Fourth invariant. Computes $\chi^{(i)}_{j,k,x}$ for a fixed triple, scanning $x$ |

The three modules share the same PDB parser, weighted Laplacian, and Schur (Kron) ligand reduction. They run from the command line and prompt for their inputs:

```bash
python SECONDS.py
python THIRDS.py
python FOURTHS.py
```

`SECONDS.py` is also imported by the AlphaFold 3 notebook, so it must be present in the Colab working directory when that notebook runs.

---

## Before you start: a note about GitHub buttons

This page is hosted on GitHub. You do **not** need to know GitHub or Git to use the code.

GitHub has several buttons that can be confusing:

* **Copy** copies text to your clipboard. It does **not** download a runnable file.
* **Raw** shows the plain file contents in the browser. For a notebook, this looks like JSON text.
* **Download raw file** downloads the actual file to your computer.
* **Download ZIP** downloads the whole repository.

To run a notebook, use **Download raw file**, **Download ZIP**, or open it directly in Google Colab.

---

## Option 1 — Run a notebook in Google Colab

This is the easiest option for readers who do not use GitHub.

1. Open the notebook you want:

   * [`compute_chi_ijk.ipynb`](compute_chi_ijk.ipynb)
   * [`compute_Rijkl.ipynb`](compute_Rijkl.ipynb)
   * [`AlphaFold3_5HED_25_SECONDS_ensemble_UPLOAD_FIRST.ipynb`](AlphaFold3_5HED_25_SECONDS_ensemble_UPLOAD_FIRST.ipynb)

2. Click the **Open in Colab** button at the top of the notebook.

3. If you do not see that button:

   * click the download icon, usually labeled **Download raw file**
   * save the notebook to your computer
   * open [Google Colab](https://colab.research.google.com/)
   * choose **File → Upload notebook**
   * upload the `.ipynb` file

4. Run the cells from top to bottom.

Colab shows a warning when a notebook is loaded from GitHub. That warning appears for every GitHub-hosted notebook. Choose **Run anyway** to continue.

---

## Option 2 — Download the whole repository

1. Go to the main repository page.
2. Click the green **`<> Code`** button.
3. Choose **Download ZIP**.
4. Unzip the file on your computer.
5. Open a notebook in Colab using the steps in Option 1, or run a `.py` module locally.

---

## Option 3 — Use Git, if you already know Git

```bash
git clone https://github.com/fatmasenguler/laplacian-minor-hierarchy.git
cd laplacian-minor-hierarchy
```

---

## Notebook 1: `compute_chi_ijk.ipynb`

Computes the third-order cooperation index $\chi_{ijk}$.

### Input

Upload exactly two files with the extension `.pdb`.

### What it computes

For each structure, the notebook builds a $C_\alpha$ contact network and evaluates $\chi_{ijk}$ for every ordered triple of residues.

### Output

CSV tables. For large proteins the output is split into several files so that each stays within the row limit of spreadsheet software.

---

## Notebook 2: `compute_Rijkl.ipynb`

Computes the fourth-order invariant $R_{ijkl}$. It analyzes one structure at a time.

### Input

The notebook prompts for:

1. the PDB file path, for example `5HED.pdb`
2. the chain ID, for example `A` (blank uses all chains)
3. residue $i$, using the PDB residue number
4. residue $j$
5. residue $k$
6. the weighting scheme

The weighting schemes are:

| Scheme | Edge weight |
| ------ | ----------- |
| `unweighted` | $w_{ij} = 1$ |
| `inv_d2` | $w_{ij} = 1/d_{ij}^{2}$ |
| `exp` | $w_{ij} = \exp(-d_{ij}/kT)$ |

For the exponential scheme, `kT` must be in Angstrom units, because the exponent has to be dimensionless. If `kT` is left blank, the mean contact distance is used.

### What it computes

With $i,j,k$ fixed, the notebook scans every remaining residue $l$ and reports $R_{ijkl}$.

### Output

A printed table and a text file, for example:

```text
5HED_chainA_exp_kT1_Rijkl_i330_j327_k372.txt
```

The file records the selected residues, the weighting scheme, diagnostic information, and the computed values.

---

## Notebook 3: `AlphaFold3_5HED_25_SECONDS_ensemble_UPLOAD_FIRST.ipynb`

Predicts 5HED with AlphaFold 3 and applies the second-order analysis to the resulting ensemble.

### Requirements

Upload `5HED.pdb` in the first cell, and make sure [`SECONDS.py`](SECONDS.py) is present in the Colab working directory. The notebook imports it and calls its functions without changing their formulas.

### What it computes

Up to 25 AF3 models (5 seeds x 5 diffusion samples). For each model the peptide ligand in chain B is eliminated by Schur reduction, and the residue profile

$$R_{i,\mathrm{cumulative}} = \sqrt{\tfrac{1}{N}\sum_j R_{ij}}$$

is computed for chain A. The notebook then reports per-residue SD and CV, a mean $\pm$ 1.96 SD ensemble band, and model-wise NRMSE and correlation to the mean. The experimental reference is the same `5HED.pdb` uploaded in the first cell.

---

## The method in brief

The protein is represented as a $C_\alpha$ contact graph with weighted Laplacian $L$. When a ligand is present, its nodes are eliminated by Schur (Kron) reduction, which leaves an effective protein Laplacian carrying the ligand-mediated couplings.

From the Moore-Penrose pseudoinverse $K = L^{+}$, the pairwise effective distance is

$$R_{ab} = K_{aa} + K_{bb} - 2K_{ab}.$$

The overlap of two paths anchored at a reference residue $i$ is

$$K^{(i)}_{ab} = \tfrac{1}{2}\left(R_{ia} + R_{ib} - R_{ab}\right).$$

The third-order cooperation index is

$$\chi_{ijk} = \frac{\left(K^{(i)}_{jk}\right)^2}{R_{ij} R_{ik}}.$$

The fourth-order normalized invariant is

$$R_{ijkl} = 1 - \frac{\det G^{(i)}_{jkl}}{R_{ij} R_{ik} R_{il}},$$

where

$$G^{(i)}_{jkl} = \begin{pmatrix}
R_{ij} & K^{(i)}_{jk} & K^{(i)}_{jl} \\
K^{(i)}_{jk} & R_{ik} & K^{(i)}_{kl} \\
K^{(i)}_{jl} & K^{(i)}_{kl} & R_{il}
\end{pmatrix}.$$

Both indices lie in $[0,1]$ by the Cauchy-Schwarz and Hadamard bounds, and both follow the sign convention stated at the top of this page: 0 for independent pathways, 1 for fully overlapping ones.

Two boundary cases are worth noting. At third order, $\chi_{ijk} = 1$ is reached only in the limit $R_{jk} \to 0$, so distinct residues approach but do not attain the upper end. At fourth order, $R_{ijkl} = 1$ is attained whenever the three displacement vectors become coplanar, which does occur for distinct residues.

See the manuscript Methods and Appendix C for the full derivation.

---

## Anchoring

The command-line modules use an explicit reference residue. The normalizer contains only distances rooted at that anchor, so $\chi^{(i)}_{j,k}$ is deliberately not symmetric under $i \leftrightarrow j$: the exchange ratio is $R_{jk}/R_{ik}$. When scanning a running residue, the anchor stays fixed. Moving the anchor onto the running residue changes the reported profile.

---

## Python packages

In Google Colab the required packages are handled inside the notebook environment. Locally:

```bash
pip install numpy scipy pandas biopython
```

`SECONDS.py`, `THIRDS.py`, and `FOURTHS.py` need only NumPy. The fourth-order notebook installs Biopython itself.

---

## Structures

PSD-95 PDZ3 variants, all T-2F bound, from the [RCSB PDB](https://www.rcsb.org/):

| PDB ID | Variant |
| ------ | ------- |
| 5HED | Wild type |
| 5HF1 | G330T |
| 5HFC | H372A |
| 5HFF | G330T + H372A |

Download the PDB files from RCSB and upload them when prompted.

---

## Parameters used in the paper

To reproduce the manuscript figures, use the parameter choices given in the Methods section: contact cutoff distance, edge-weighting rule, sequence-neighbor exclusion rule, chain selection, and residue indexing convention. Use identical settings for every structure being compared.

---

## Troubleshooting

### Colab warns that the notebook is loaded from GitHub

This warning appears for every GitHub-hosted notebook. Review the code if you wish, then choose **Run anyway**.

### The AlphaFold 3 notebook says `SECONDS.py` was not found

Upload [`SECONDS.py`](SECONDS.py) to the Colab working directory, then rerun the analysis cell. The notebook looks for `SECONDS.py`, `SECONDS*.py`, `seconds.py`, or `seconds*.py`.

### I clicked Copy, but I did not get a notebook file

The **Copy** button only copies text. Use **Download raw file**, **Download ZIP**, or **Open in Colab**.

### The notebook opens as strange text

An `.ipynb` file is stored internally as JSON. Clicking **Raw** shows that JSON. Download the file and open it in Colab.

### Colab asks me to upload files

That is expected. `compute_chi_ijk.ipynb` needs exactly two `.pdb` files. `compute_Rijkl.ipynb` needs the PDB file present in the session; enter its path when prompted.

### The notebook says a residue was not found

Enter the PDB residue number, not the row number in the output. Check that the chain ID is correct.

### The contact graph is disconnected

The cutoff may be too small, or the selected chain may not form a connected contact graph. Check the chain ID and the PDB file, and raise the cutoff if necessary.

### A script refuses to run because residue numbers repeat

The modules reject duplicate residue numbers rather than silently keeping the last occurrence, because residue matching between two structures would be ambiguous. Restrict the chain selection or renumber the structure.

---

## Citation and license

If you use this code, please cite the paper above.

Released under the MIT License. See [`LICENSE`](LICENSE).
