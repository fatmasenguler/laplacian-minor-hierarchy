#!/usr/bin/env python3
"""
compare_linearized_cumulative_resistance.py
===========================================

Compares the second Laplacian invariant, effective resistance R_ij,
for two proteins through a residue-level summary profile:
  1. Computes R_ij accurately for each protein.
  2. Matches the two proteins by residue number.
  3. Averages R_ij only over residues present in BOTH structures.
  4. Reports sqrt(mean_common R_ij) as a derived dimensionless profile.
  5. Exports 4 columns:
     Residue Index | sqrt(mean R) (Protein 1) | sqrt(mean R) (Protein 2) | Difference (P2 - P1)

Note: R_ij is the second invariant. The reported profile is averaged over
the SAME common residue set in both structures. With w_ij = exp(-d_ij/d0),
both quantities are dimensionless.

Ligand handling follows effective_laplacian.py: peptide ligands are read
from ligand chains (default B), while non-water/non-ion HETATM ligands
are read as heavy-atom nodes and eliminated by Schur reduction.

Dependencies: numpy only
"""

import os
import sys
import warnings
import numpy as np

# ======================================================================
# 1. Structure parsing
# ======================================================================

WATER = {"HOH", "WAT", "DOD", "H2O"}
IONS = {"NA", "CL", "K", "MG", "CA", "ZN", "SO4", "PO4", "GOL", "EDO", "ACT"}

def read_complex(pdb_path, ligand_resnames=None, protein_chains=None,
                 ligand_chains=None, skip_ions=True, model=1):
    """
    Read protein CA atoms and ligand heavy atoms from a PDB file.

    Ligands appear in PDB files in two different ways, so there are two
    ways to select them:

        small molecule / cofactor / drug -> HETATM -> ligand_resnames
        peptide ligand                   -> ATOM   -> ligand_chains

    A peptide ligand cannot be selected by residue name, because its
    residues carry ordinary amino-acid names indistinguishable from the
    protein's own. It is selected by chain instead, and ALL of its heavy
    atoms become nodes, not just CA. The protein is represented by CA
    atoms only, so the two node sets have different granularity by
    design: the ligand is usually small, and its internal geometry is
    what mediates the new pathways.

    Parameters
    ----------
    pdb_path : str
        Path to a PDB-format file.
    ligand_resnames : list of str, optional
        HETATM residue names to treat as ligand. If None, every HETATM
        other than water and, when skip_ions is True, the additives in
        IONS. An explicit name here overrides skip_ions.
    protein_chains : list of str, optional
        Chain IDs forming the protein. If None, every chain that is not
        listed in ligand_chains.
    ligand_chains : list of str, optional
        Chain IDs holding a peptide ligand.
    skip_ions : bool
        Exclude the entries of IONS when ligand_resnames is not given.
    model : int
        Which MODEL to read from a multi-model file such as an NMR
        ensemble. Files without MODEL records are unaffected. Without
        this filter every atom of such a file would be read once per
        model, silently multiplying the node set.

    Notes
    -----
    Only altloc ' ' or 'A' is kept, so a residue modelled in two
    conformations contributes a single node. Hydrogen and deuterium are
    excluded; when the element column is blank the atom name is used as
    a fallback.

    Returns
    -------
    protein_xyz : ndarray, shape (n_p, 3)
    protein_ids : list of str
        Labels of the form "chain:resnum:resname". The residue field may
        carry an insertion code, e.g. "A:100A:ALA"; use
        parse_residue_numbers to convert them safely.
    ligand_xyz : ndarray, shape (n_l, 3)
    ligand_ids : list of str
        Labels of the form "chain:resnum:resname:atomname".
    """
    ligand_chains = set(ligand_chains or [])
    p_xyz, p_ids, l_xyz, l_ids = [], [], [], []
    current_model, seen_model_record = 1, False

    with open(pdb_path) as fh:
        for line in fh:
            rec = line[:6].strip()

            if rec == "MODEL":
                seen_model_record = True
                try:
                    current_model = int(line[10:14])
                except ValueError:
                    current_model += 1
                continue
            if rec == "ENDMDL":
                continue
            if rec not in ("ATOM", "HETATM"):
                continue
            if seen_model_record and current_model != model:
                continue
            if line[16] not in (" ", "A"):                    # altloc
                continue

            name = line[12:16].strip()
            resname = line[17:20].strip()
            chain = line[21].strip()
            resseq = line[22:26].strip()
            icode = line[26].strip()
            element = line[76:78].strip().upper()
            is_hydrogen = (element in ("H", "D")
                           or (not element and name[:1] in ("H", "D")))

            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            label = f"{chain}:{resseq}{icode}:{resname}"

            if rec == "ATOM" and chain in ligand_chains:
                # peptide ligand: every heavy atom is a node
                if is_hydrogen:
                    continue
                l_xyz.append(xyz)
                l_ids.append(f"{label}:{name}")

            elif rec == "ATOM":
                if protein_chains and chain not in protein_chains:
                    continue
                if name == "CA":
                    p_xyz.append(xyz)
                    p_ids.append(label)

            else:                                             # HETATM
                if is_hydrogen or resname in WATER:
                    continue
                if ligand_resnames is not None:
                    if resname not in ligand_resnames:
                        continue
                elif skip_ions and resname in IONS:
                    continue
                l_xyz.append(xyz)
                l_ids.append(f"{label}:{name}")

    p_arr = np.asarray(p_xyz, float).reshape(-1, 3)
    l_arr = np.asarray(l_xyz, float).reshape(-1, 3)
    if len(p_arr) == 0:
        warnings.warn(f"{pdb_path}: no protein CA atoms selected. "
                      f"Check the protein_chains argument.")
    return p_arr, p_ids, l_arr, l_ids


def parse_residue_numbers(ids):
    """Convert labels into integer numbers and names."""
    nums, names = [], []
    for s in ids:
        field = s.split(":")[1]
        digits = "".join(c for c in field if c.isdigit() or c == "-")
        nums.append(int(digits))
        names.append(s.split(":")[2])
    nums = np.asarray(nums, int)
    return nums, np.asarray(names)


# ======================================================================
# 2. Laplacian & Schur reduction
# ======================================================================

def build_laplacian(coords, rc=7.8, d0=1.0):
    """Distance-dependent Laplacian of a point set."""
    if rc <= 0:
        raise ValueError("rc must be positive")
    if d0 <= 0:
        raise ValueError("d0 must be positive")

    coords = np.asarray(coords, float)
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError(f"coords must have shape (n, 3), got {coords.shape}")
    if len(coords) == 0:
        raise ValueError("cannot build a Laplacian from an empty node set")

    diff = coords[:, None, :] - coords[None, :, :]
    d = np.sqrt((diff ** 2).sum(-1))

    off = -np.exp(-d / d0)
    off[d > rc] = 0.0
    np.fill_diagonal(off, 0.0)

    L = off.copy()
    np.fill_diagonal(L, -off.sum(axis=1))
    return L


def schur_reduce(L, n_protein, cond_warn=1e12):
    """Kron reduction: eliminate the trailing block of L."""
    Hpp = L[:n_protein, :n_protein]
    Hpl = L[:n_protein, n_protein:]
    Hlp = L[n_protein:, :n_protein]
    Hll = L[n_protein:, n_protein:]

    if Hll.shape[0] == 0:
        return Hpp.copy()

    cond = np.linalg.cond(Hll)
    if not np.isfinite(cond) or cond > cond_warn:
        warnings.warn(f"H_ll is ill-conditioned (cond = {cond:.2e}). "
                      f"Falling back to least-squares solve.")
        X = np.linalg.lstsq(Hll, Hlp, rcond=None)[0]
    else:
        X = np.linalg.solve(Hll, Hlp)

    return Hpp - Hpl @ X


def build_effective_laplacian(pdb_path, ligand_resnames=None,
                              protein_chains=("A",), ligand_chains=("B",),
                              rc=7.8, d0=1.0, model=1, verbose=False):
    """Build the protein effective Laplacian including ligand influence.

    Default PDZ convention:
        protein chain  = A
        peptide ligand = B (all heavy ATOM records)

    Small-molecule/cofactor/drug ligands written as HETATM records are also
    read by the same parser.  ``ligand_resnames`` can be supplied to restrict
    HETATM selection to specific residue names.
    """
    p_xyz, p_ids, l_xyz, l_ids = read_complex(
        pdb_path,
        ligand_resnames=ligand_resnames,
        protein_chains=list(protein_chains) if protein_chains is not None else None,
        ligand_chains=list(ligand_chains) if ligand_chains is not None else None,
        model=model,
    )

    n_p, n_l = len(p_xyz), len(l_xyz)
    if verbose:
        print(f"{pdb_path}: {n_p} residues (CA), {n_l} ligand heavy atoms")
        if n_l:
            # Show what was actually recognized as ligand.  This makes an
            # incorrect chain/residue selection immediately visible.
            ligand_residues = []
            seen = set()
            for s in l_ids:
                fields = s.split(":")
                key = ":".join(fields[:3])
                if key not in seen:
                    seen.add(key)
                    ligand_residues.append(key)
            preview = ", ".join(ligand_residues[:12])
            if len(ligand_residues) > 12:
                preview += ", ..."
            print(f"  ligand residues recognized   : {preview}")

    if n_l == 0:
        if verbose:
            print("  no ligand found: returning plain protein Laplacian")
        return build_laplacian(p_xyz, rc, d0), p_ids

    # Protein nodes first, ligand heavy atoms last: required by schur_reduce.
    coords = np.vstack([p_xyz, l_xyz])
    L_full = build_laplacian(coords, rc, d0)
    L_eff = schur_reduce(L_full, n_p)
    L_eff = (L_eff + L_eff.T) / 2

    if verbose:
        # Number of protein-protein couplings that arise only after ligand
        # elimination.  A positive value confirms ligand-mediated coupling.
        off_eff = np.abs(L_eff - np.diag(np.diag(L_eff))) > 1e-10
        direct = L_full[:n_p, :n_p]
        off_dir = np.abs(direct - np.diag(np.diag(direct))) > 1e-10
        created = int((off_eff & ~off_dir).sum() // 2)
        print(f"  ligand-mediated edges created: {created}")

    return L_eff, p_ids


# ======================================================================
# 3. Invariants & Resistance Calculation
# ======================================================================

def _require_connected_laplacian(L):
    """Raise if L is not numerically a connected Laplacian."""
    L = np.asarray(L, float)
    if L.ndim != 2 or L.shape[0] != L.shape[1]:
        raise ValueError("L must be a square matrix")
    if L.shape[0] < 2:
        raise ValueError("L must contain at least two nodes")

    Ls = (L + L.T) / 2
    eig = np.linalg.eigvalsh(Ls)
    scale = max(float(abs(eig[-1])), 1.0)
    tol = 1e-12 * scale
    if eig[1] <= tol:
        raise ValueError(
            "The protein network is disconnected at the chosen rc/d0. "
            "Effective resistance between components is not finite."
        )


def effective_distances(L):
    """
    Compute the second invariant, pairwise effective resistance R_ij.

    For a connected Laplacian,
        L^+ = (L + J)^(-1) - J,  J = 11^T/N,
    and R_ij = L^+_ii + L^+_jj - 2 L^+_ij.
    """
    _require_connected_laplacian(L)
    N = L.shape[0]
    J = np.ones((N, N)) / N
    K = np.linalg.solve(L + J, np.eye(N)) - J
    
    diag = np.diag(K)
    R = diag[:, None] + diag[None, :] - 2 * K
    
    # Clip tiny negative numerical artifact values (e.g. -1e-16) to 0.0
    return np.clip(R, 0.0, None)


def compute_sqrt_mean_R_from_R(R, column_indices):
    """Derived residue profile using the same common residue set in both structures.

    For each residue i:
        profile_i = sqrt((1/N_common) * sum_{j in common} R_ij)

    R is computed on the complete network of each structure. Residues present
    only in one structure may therefore still affect the network geometry, but
    they do not enter the averaging sum for only one of the two structures.
    """
    column_indices = np.asarray(column_indices, dtype=int)
    if column_indices.size == 0:
        raise ValueError("column_indices must contain at least one common residue")

    cum_R = np.mean(R[:, column_indices], axis=1)
    return np.sqrt(cum_R)

# ======================================================================
# 4. Interactive Input Functions
# ======================================================================

def _prompt_pdb_file(label):
    """Prompt explicitly for one PDB file and require that it exists."""
    while True:
        pdb_file = input(f"\nEnter PDB filename for {label} (e.g., 5HED.pdb): ").strip()
        if not pdb_file:
            print("Please enter a filename.")
            continue
        if not os.path.exists(pdb_file):
            print(f"File '{pdb_file}' not found in current working directory.")
            print(f"Current directory: {os.getcwd()}")
            continue
        return pdb_file


def get_inputs():
    """Prompt separately for Protein 1, Protein 2, cutoff, and d0."""
    print("\n" + "="*60)
    print("SECOND INVARIANT: EFFECTIVE RESISTANCE COMPARISON")
    print("="*60)

    # Restored explicit interactive prompts for the two structures.
    pdb1 = _prompt_pdb_file("Protein 1")
    pdb2 = _prompt_pdb_file("Protein 2")

    print("\n" + "-"*60)
    print("PARAMETERS")
    print("-"*60)

    while True:
        try:
            rc = float(input("Enter cutoff distance in Å (default 7.8): ").strip() or "7.8")
            if rc <= 0:
                print("Cutoff must be positive.")
                continue
            break
        except ValueError:
            print("Please enter a valid number.")

    while True:
        try:
            d0 = float(input("Enter temperature/length scale d0 (default 1.0): ").strip() or "1.0")
            if d0 <= 0:
                print("d0 must be positive.")
                continue
            break
        except ValueError:
            print("Please enter a valid number.")

    return pdb1, pdb2, rc, d0


# ======================================================================
# 5. Main Execution
# ======================================================================

def main():
    np.seterr(all='ignore')
    
    pdb1, pdb2, rc, d0 = get_inputs()
    
    base1 = os.path.splitext(pdb1)[0]
    base2 = os.path.splitext(pdb2)[0]
    output_file = f"{base1}_vs_{base2}_sqrt_mean_R.txt"
    
    print(f"\nProcessing proteins...")
    print(f"  Protein 1: {pdb1}")
    print(f"  Protein 2: {pdb2}")
    print(f"  rc = {rc} Å, d0 = {d0} Å")
    
    try:
        print("\nBuilding Laplacian for Protein 1...")
        L1, ids1 = build_effective_laplacian(pdb1, rc=rc, d0=d0, verbose=True)
        
        print("\nBuilding Laplacian for Protein 2...")
        L2, ids2 = build_effective_laplacian(pdb2, rc=rc, d0=d0, verbose=True)
    except Exception as e:
        print(f"\nError building Laplacians: {e}")
        return

    resnums1, _ = parse_residue_numbers(ids1)
    resnums2, _ = parse_residue_numbers(ids2)

    # Residue numbers must be unique within each selected chain.
    if len(set(resnums1.tolist())) != len(resnums1):
        raise ValueError("Protein 1 has duplicate residue numbers; cannot match safely.")
    if len(set(resnums2.tolist())) != len(resnums2):
        raise ValueError("Protein 2 has duplicate residue numbers; cannot match safely.")

    # Match residues BEFORE constructing the residue-level summary profile.
    # R is still computed on each COMPLETE reduced protein network.
    idx1 = {int(r): i for i, r in enumerate(resnums1)}
    idx2 = {int(r): i for i, r in enumerate(resnums2)}
    common_res = sorted(set(idx1) & set(idx2))
    if not common_res:
        raise ValueError("The two proteins have no residue numbers in common.")

    common_idx1 = np.asarray([idx1[r] for r in common_res], dtype=int)
    common_idx2 = np.asarray([idx2[r] for r in common_res], dtype=int)

    print("\nCalculating effective resistance matrices...")
    R1 = effective_distances(L1)
    R2 = effective_distances(L2)

    print(f"Calculating residue profiles over {len(common_res)} common residues...")
    profile1 = compute_sqrt_mean_R_from_R(R1, common_idx1)
    profile2 = compute_sqrt_mean_R_from_R(R2, common_idx2)

    only1 = sorted(set(idx1) - set(idx2))
    only2 = sorted(set(idx2) - set(idx1))
    if only1 or only2:
        print(f"\nResidues matched by number: {len(common_res)} common residues.")
        if only1:
            print(f"  Present only in Protein 1: {only1}")
        if only2:
            print(f"  Present only in Protein 2: {only2}")

    with open(output_file, 'w') as f:
        f.write("# Effective-resistance residue profile\n")
        f.write("# Second invariant: R_ij\n")
        f.write("# Reported profile: sqrt((1/N_common) * sum_{j in common} R_ij), dimensionless\n")
        f.write(f"# Protein 1: {pdb1}\n")
        f.write(f"# Protein 2: {pdb2}\n")
        f.write(f"# Parameters: rc = {rc} Å, d0 = {d0} Å\n")
        f.write("# Difference = Protein 2 - Protein 1\n")
        f.write("# Both profiles are averaged over the same common residue set.\n")
        f.write("# --------------------------------------------------------\n")
        f.write(f"# {'Res_ID':<8} {'P1_sqrtMeanR':<16} {'P2_sqrtMeanR':<16} {'Diff_P2-P1':<16}\n")

        for res_id in common_res:
            p1_val = profile1[idx1[res_id]]
            p2_val = profile2[idx2[res_id]]
            d_val = p2_val - p1_val
            f.write(f"{res_id:<8d} {p1_val:<16.6f} {p2_val:<16.6f} {d_val:<16.6f}\n")

    print("\n" + "="*60)
    print("CALCULATION COMPLETE")
    print("="*60)
    print(f"Results saved to: {output_file}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nProcess canceled by user.")
        sys.exit(0)