# The Geometry of Allostery: Quantifying Pathway Organization in Protein Networks

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22172172.svg)](https://doi.org/10.5281/zenodo.22172172)

Code accompanying:

> Senguler Ciftci F. and Erman B. *The Geometry of Allostery: Quantifying Pathway Organization in Protein Networks.* (2026)

This repository computes higher-order effective-distance invariants of a protein contact network directly from PDB structures.

The hierarchy has three levels:

* the second-order effective distance $R_{ij}$ (Eq. 7)
* the third-order decoupling index $\chi_{ijk}$ (Eq. 15)
* the fourth-order decoupling index $\chi_{ijkl}$ (Eq. 19)

**Sign convention.** Both decoupling indices run from 0 to 1, with

```text
0  coupled     (complete pathway overlap, cooperative multi-node unit)
1  decoupled   (independent, topologically disjoint communication channels)
```

This is the convention of the manuscript (Section 2.1.3 and Appendix C.6). Note that it is the *decoupling* direction: a larger value means *less* correlation between pathways.

---

## Audio summary, appendices, and supplementary note

A spoken summary of the paper, the full appendices, and Supplementary Note S1 (AlphaFold 3 ensemble consistency of the 5HED effective-distance geometry) are archived on Zenodo:

**https://doi.org/10.5281/zenodo.22172172**

---

## What is in this repository?

### Colab notebooks

| File | What it does |
| ---- | ------------ |
| [`compute_chi_ijk.ipynb`](compute_chi_ijk.ipynb) | Third-order index $\chi_{ijk}$ for all ordered triples in each of two uploaded PDB files |
| [`compute_Rijkl.ipynb`](compute_Rijkl.ipynb) | Fourth-order index $\chi_{ijkl}$ for a fixed triple $i,j,k$ and every remaining residue $l$ |
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

## Module 1: `SECONDS.py`

Computes the second invariant, the pairwise effective distance $R_{ij}$, and reduces it to the residue profile of Eq. (10):

$$R_{i,\mathrm{cumul}} = \sqrt{\frac{1}{N}\sum_j R_{ij}}.$$

It compares two structures and reports the profile of each plus their difference. The average is taken over the residues present in both structures, so the two columns are on the same footing.

To reproduce Figure 2, run it once for each variant against the apo reference `1BFE.pdb`.

---

## Module 2: `THIRDS.py`

Computes the third-order decoupling index

$$\chi^{(i)}_{j,k} = 1 - \frac{\left(K^{(i)}_{jk}\right)^2}{R_{ij}R_{ik}}.$$

Residue $i$ is the anchor and $j$ the partner; both stay fixed while $k$ runs over the sequence. The script prompts for two PDB files, the chain IDs, $r_c$, $d_0$, and then $i$ and $j$.

To reproduce Figure 3, use $i = 330$ (the $\beta_2$–$\beta_3$ surface loop) and $j = 372$ (the binding pocket), with the structure pairs 1BFE→5HEB, 5HEB→5HEY, 5HEB→5HED, 5HEB→5HFB, 5HFC→5HFF, and 5HEB→5HF1.

At $k = i$ the index is undefined, because $R_{ii} = 0$ sits in the denominator. At $k = j$ it is exactly 0.

---

## Module 3: `FOURTHS.py`

Computes the fourth-order decoupling index

$$\chi^{(i)}_{j,k,l} = \frac{\det K_S}{R_{ij}R_{ik}R_{il}}, \qquad
K_S = \begin{pmatrix}
R_{ij} & K^{(i)}_{jk} & K^{(i)}_{jl} \\
K^{(i)}_{jk} & R_{ik} & K^{(i)}_{kl} \\
K^{(i)}_{jl} & K^{(i)}_{kl} & R_{il}
\end{pmatrix}.$$

The triple $i,j,k$ stays fixed while the fourth residue runs over the sequence. The script prompts for two PDB files, the chain IDs, $r_c$, $d_0$, and then $i$, $j$, $k$.

To reproduce Figure 4, fix the structural triad $i = 320$, $j = 330$, $k = 372$ and sweep the fourth node, using the pairs 5HEB→5HED, 5HEB→5HEY, 5HEB→5HF1, 5HFC→5HFF, and 5HEB→5HFB.

At $x = i$ the index is undefined. At $x = j$ and $x = k$ two rows of $K_S$ coincide, so $\det K_S = 0$ and the index is exactly 0.

---

## Notebook: `AlphaFold3_5HED_25_SECONDS_ensemble_UPLOAD_FIRST.ipynb`

Predicts 5HED with AlphaFold 3 and applies the second-order analysis to the resulting ensemble. This is the calculation behind Supplementary Note S1.

### Requirements

Upload `5HED.pdb` in the first cell, and make sure [`SECONDS.py`](SECONDS.py) is present in the Colab working directory. The notebook imports it and calls its functions without changing their formulas.

### What it computes

Up to 25 AF3 models (5 seeds × 5 diffusion samples). For each model the peptide ligand in chain B is eliminated by Schur reduction, and the residue profile $R_{i,\mathrm{cumul}}$ is computed for chain A. The notebook then reports per-residue SD and CV, a mean ± 1.96 SD ensemble band, and model-wise NRMSE and correlation to the mean. The experimental reference is the same `5HED.pdb` uploaded in the first cell.

---

## The method in brief

The protein is represented as a $C_\alpha$ contact graph with weighted Laplacian $L$, whose off-diagonal entries are $L_{ij} = -\exp(-d_{ij}/d_0)$ for $d_{ij} \le r_c$ and zero otherwise (Eq. 1). When a ligand is present, its nodes are eliminated by Schur (Kron) reduction, which leaves an effective protein Laplacian carrying the ligand-mediated couplings (Eq. 21):

$$L_\mathrm{eff} = L_{pp} - L_{pl}L_{ll}^{-1}L_{lp}.$$

From the Moore-Penrose pseudoinverse $K = L^{+}$, the pairwise effective distance is

$$R_{ab} = K_{aa} + K_{bb} - 2K_{ab}.$$

The overlap of two paths anchored at a reference residue $i$ is

$$K^{(i)}_{ab} = \tfrac{1}{2}\left(R_{ia} + R_{ib} - R_{ab}\right).$$

The third-order decoupling index is

$$\chi_{ijk} = \frac{R_{ijk}}{R_{ij}R_{ik}} = 1 - \frac{\left(K^{(i)}_{jk}\right)^2}{R_{ij}R_{ik}},$$

where $R_{ijk} = \det L(i,j,k)/\tau$ is the third Laplacian minor, equal to the $2 \times 2$ Gram determinant. The fourth-order decoupling index is

$$\chi_{ijkl} = \frac{R_{ijkl}}{R_{ij}R_{ik}R_{il}},$$

where $R_{ijkl} = \det L(i,j,k,l)/\tau = \det K_S$ is the $3 \times 3$ Gram determinant given above.

Note the distinction: $R_{ijk}$ and $R_{ijkl}$ are the *unnormalized* minors; $\chi_{ijk}$ and $\chi_{ijkl}$ are the normalized indices. Only the latter are bounded in $[0,1]$, by the Cauchy-Schwarz and Hadamard inequalities respectively.

Two boundary cases are worth noting. At third order, $\chi = 0$ requires the two displacement vectors anchored at $i$ to be collinear, which happens exactly when $k = j$. At fourth order, $\chi = 0$ whenever the three displacement vectors become linearly dependent, which occurs at $x = j$ and $x = k$ and can also occur for distinct residues. The upper end, $\chi = 1$, requires mutual orthogonality of the anchored displacement vectors.

See the manuscript Methods and Appendix C for the full derivation.

---

## Anchoring

The command-line modules use an explicit reference residue. The unnormalized minor $R_S = \det L(S)/\tau$ is symmetric under any permutation of the residues in $S$, but the normalizer contains only distances rooted at the anchor. Consequently $\chi^{(i)}_{j,k}$ is deliberately not symmetric under $i \leftrightarrow j$, and the exchange ratio is exactly

$$\frac{\chi^{(i)}_{j,k}}{\chi^{(j)}_{i,k}} = \frac{R_{jk}}{R_{ik}}.$$

The index *is* symmetric in the target indices: $\chi^{(i)}_{j,k} = \chi^{(i)}_{k,j}$, and at fourth order $\chi_{ijkl} = \chi_{ikjl} = \chi_{iljk}$.

When scanning a running residue, the anchor stays fixed. Moving the anchor onto the running residue changes the reported profile.

---

## Python packages

In Google Colab the required packages are handled inside the notebook environment. Locally:

```bash
pip install numpy scipy pandas biopython
```

`SECONDS.py`, `THIRDS.py`, and `FOURTHS.py` need only NumPy. The fourth-order notebook installs Biopython itself.

---

## Structures

PSD-95 PDZ3 variants from the [RCSB PDB](https://www.rcsb.org/). The manuscript analyzes all eight ligand-bound structures characterized by Raman et al. (2016), evaluated against the apo wild-type reference:

| PDB ID | State |
| ------ | ----- |
| 1BFE | Apo wild type (baseline reference) |
| 5HEB | Wild type, Class I (CRIPT) bound |
| 5HED | Wild type, Class II (T-2F) bound |
| 5HET | Wild type variant |
| 5HEY | G330T, Class I bound |
| 5HF1 | G330T, Class II bound |
| 5HFB | H372A |
| 5HFC | H372A, Class II bound |
| 5HFF | G330T + H372A, Class II bound |

Download the PDB files from RCSB and upload them when prompted.

---

## Parameters used in the paper

| Parameter | Value |
| --------- | ----- |
| Contact cutoff $r_c$ | 7.8 Å between $C_\alpha$ atoms |
| Length scale $d_0$ | 1.0 Å |
| Edge weight | $w_{ij} = \exp(-d_{ij}/d_0)$ |
| Protein chain | A ($C_\alpha$ atoms only) |
| Peptide ligand chain | B (all heavy atoms, eliminated by Schur reduction) |

These are the defaults in every module, so pressing Enter at each prompt reproduces the manuscript settings. Use identical settings for every structure being compared.

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

### The index comes out near 1 where I expected near 0

Check the sign convention at the top of this page. These are *decoupling* indices: 0 means overlapping pathways, 1 means independent ones. An earlier version of this repository used the complementary quantity $1 - \chi$.

---

## Citation and license

If you use this code, please cite the paper and the archived record:

> Senguler Ciftci F. and Erman B. *The Geometry of Allostery: Quantifying Pathway Organization in Protein Networks.* (2026)

> Audio summary, appendices, and Supplementary Note S1: https://doi.org/10.5281/zenodo.22172172

Released under the MIT License. See [`LICENSE`](LICENSE).
---

## Citation and license

If you use this code, please cite the paper above.

Released under the MIT License. See [`LICENSE`](LICENSE).
