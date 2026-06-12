# laplacian-minor-hierarchy
# Laplacian Minor Hierarchy for Many-Body Protein Communication

Code accompanying:

> Senguler Ciftci F and Erman B. *The Geometry of Allostery: A Laplacian Minor
> Hierarchy for Many-Body Protein Communication.* (2026)

The scripts compute the higher-order effective-distance invariants of a protein
contact network from a PDB structure: the third-order cooperation index
(chi_ijk) and the fourth-order normalized invariant (chi_ijkl).

## Contents

| File | Computes | Interface |
|------|----------|-----------|
| `src/compute_chi_ijk.py` | Third-order cooperation index chi_ijk for all ordered triples of two PDB files | Google Colab (file upload) |
| `src/compute_Rijkl.py`   | Fourth-order normalized invariant for fixed i, j, k and all l | Command line (interactive) |

## Method

The protein is represented as a Calpha contact graph with weighted Laplacian L.
From the Moore-Penrose pseudoinverse K = L^+, the pairwise effective distance is

    R_ab = K_aa + K_bb - 2 K_ab.

The reference-anchored overlap is

    K_ab^(i) = (1/2) (R_ia + R_ib - R_ab),

the third-order cooperation index is

    chi_ijk = 1 - (K_jk^(i))^2 / (R_ij R_ik),

and the fourth-order normalized invariant is the determinant of the 3x3 Gram
matrix of overlaps divided by R_ij R_ik R_il. Both indices lie in [0, 1].
See the manuscript Methods and Appendix C for the full derivation.

## Installation

```bash
pip install -r requirements.txt
```

`compute_chi_ijk.py` additionally requires the Google Colab runtime
(`google.colab`); it reads structures through the Colab file-upload widget.

## Usage

### Fourth-order invariant (command line)

```bash
python src/compute_Rijkl.py
```

The script prompts for the PDB file path, chain ID, residues i / j / k, and a
weighting scheme, then reports the invariant for every remaining residue l and
writes a `.txt` table.

### Third-order index (Colab)

Open `src/compute_chi_ijk.py` in Google Colab, run all cells, and upload exactly
two PDB files when prompted. The script writes per-protein CSV tables (split into
Excel-safe chunks) and a summary, then zips them for download.

## Parameters used in the paper

To reproduce the manuscript figures, set the following so both scripts agree
with the Methods section:

- **Cutoff:** 7.8 A
- **Edge weight:** `exp(-d_ij / d0)` with **d0 = 1**
- **Sequence-neighbor exclusion:** apply the same rule in both scripts

> Note: the current defaults differ between the two scripts (chi script: 8.0 A
> cutoff and `d_mean` scaling; Rijkl script: 7.8 A and a mean-distance fallback).
> Align these before regenerating figures.

## Structures

PSD-95 PDZ3 variants, all T-2F bound, from the RCSB PDB:

| PDB | Variant |
|-----|---------|
| 5HED | Wild type |
| 5HF1 | G330T |
| 5HFC | H372A |
| 5HFF | G330T + H372A |

## Citation

If you use this code, please cite the paper above.

## License

See `LICENSE`.
