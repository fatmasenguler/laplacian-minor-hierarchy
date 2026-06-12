# Laplacian Minor Hierarchy for Many-Body Protein Communication

Code accompanying:

> Senguler Ciftci F and Erman B. *The Geometry of Allostery: A Laplacian Minor
> Hierarchy for Many-Body Protein Communication.* (2026)

The scripts compute the higher-order effective-distance invariants of a protein
contact network directly from a PDB structure: the third-order cooperation
index ($\chi_{ijk}$) and the fourth-order normalized invariant ($\chi_{ijkl}$).

---

## How to get the code

The Python scripts are in the [`src/`](src) folder. You do **not** need to know
Git to use them. Pick whichever option suits you.

### Option 1 — Download everything as a ZIP (easiest)

1. Go to the repository's main page (the page you are reading now).
2. Click the green **`<> Code`** button near the top right.
3. Choose **Download ZIP**.
4. Unzip the file on your computer. The scripts are inside the `src/` folder.

### Option 2 — Download a single script

1. Open the `src/` folder and click the file you want
   (for example `compute_Rijkl.py`).
2. On the file page, click the **download icon** at the top right of the file
   view (tooltip: *Download raw file*). If you don't see it, click **Raw** and
   then use your browser's **File → Save As**.
3. Save it with a `.py` extension and run it like any other Python script.

### Option 3 — Clone with Git (if you already use Git)

```bash
git clone https://github.com/fatmasenguler/laplacian-minor-hierarchy.git
```

---

## What's in this repository

| File | Computes | How you run it |
|------|----------|----------------|
| [`src/compute_chi_ijk.py`](src/compute_chi_ijk.py) | Third-order cooperation index $\chi_{ijk}$ for all ordered triples of two PDB files | Google Colab (upload structures in the browser) |
| [`src/compute_Rijkl.py`](src/compute_Rijkl.py)     | Fourth-order normalized invariant for a fixed triple $i,j,k$ and every remaining residue $l$ | Command line (answers a few prompts) |

---

## Installing the requirements

Once you have the files, install the Python packages they need:

```bash
pip install -r requirements.txt
```

That installs `numpy`, `scipy`, `pandas`, and `biopython`.
`compute_chi_ijk.py` also needs the Google Colab runtime (`google.colab`),
because it reads structures through the Colab file-upload widget.

---

## Running the scripts

### Fourth-order invariant — `compute_Rijkl.py` (command line)

```bash
python src/compute_Rijkl.py
```

It asks for the PDB file path, chain ID, the residues $i$, $j$, $k$, and a
weighting scheme, then prints the invariant for every other residue $l$ and
saves the table to a `.txt` file.

### Third-order index — `compute_chi_ijk.py` (Colab)

Open the file in [Google Colab](https://colab.research.google.com/), run all
cells, and upload exactly two PDB files when prompted. It writes one CSV table
per protein (split into Excel-safe chunks), a summary, and a ZIP for download.

---

## The method in brief

The protein is a $C_\alpha$ contact graph with weighted Laplacian $L$. From its
Moore–Penrose pseudoinverse $K = L^{+}$, the pairwise effective distance is

$$R_{ab} = K_{aa} + K_{bb} - 2K_{ab}.$$

The overlap of two paths anchored at a reference residue $i$ is

$$K^{(i)}_{ab} = \tfrac{1}{2}\left(R_{ia} + R_{ib} - R_{ab}\right).$$

The third-order cooperation index is

$$\chi_{ijk} = 1 - \frac{\left(K^{(i)}_{jk}\right)^{2}}{R_{ij}\,R_{ik}},$$

and the fourth-order normalized invariant is the determinant of the $3\times 3$
Gram matrix of overlaps, divided by the product of the three branch distances:

$$\chi_{ijkl} = \frac{\det G^{(i)}_{jkl}}{R_{ij}\,R_{ik}\,R_{il}},
\qquad
G^{(i)}_{jkl} =
\begin{pmatrix}
R_{ij} & K^{(i)}_{jk} & K^{(i)}_{jl}\\[2pt]
K^{(i)}_{jk} & R_{ik} & K^{(i)}_{kl}\\[2pt]
K^{(i)}_{jl} & K^{(i)}_{kl} & R_{il}
\end{pmatrix}.$$

Both indices lie in $[0,1]$. See the manuscript Methods and Appendix C for the
full derivation.

---

## Parameters used in the paper

To reproduce the manuscript figures, set both scripts to match the Methods:

- **Cutoff:** 7.8 Å
- **Edge weight:** $\exp(-d_{ij}/d_0)$ with $d_0 = 1$
- **Sequence-neighbor exclusion:** use the same rule in both scripts

> The current defaults differ slightly between the two scripts (the $\chi$ script
> uses an 8.0 Å cutoff and $d_\text{mean}$ scaling; the $R_{ijkl}$ script uses
> 7.8 Å with a mean-distance fallback). Align these before regenerating figures.

---

## Structures

PSD-95 PDZ3 variants, all T-2F bound, from the [RCSB PDB](https://www.rcsb.org/):

| PDB | Variant |
|-----|---------|
| 5HED | Wild type |
| 5HF1 | G330T |
| 5HFC | H372A |
| 5HFF | G330T + H372A |

---

## Citation and license

If you use this code, please cite the paper above. Released under the MIT
License — see [`LICENSE`](LICENSE).
