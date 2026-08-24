#!/usr/bin/env python3
"""
calculate_chi_abci.py
=====================

Interactive script to calculate the fourth-order cooperation index chi_{a,b,c,i}
for two proteins and their difference (Protein 2 - Protein 1).

Output format (4 columns):
  Residue_ID | chi_a,b,c,i (Protein 1) | chi_a,b,c,i (Protein 2) | Diff (P2 - P1)

Dependencies: numpy only
"""

import os
import sys
import warnings
import numpy as np

# ======================================================================
# 1. Structure Parsing
# ======================================================================

WATER = {"HOH", "WAT", "DOD", "H2O"}
IONS = {"NA", "CL", "K", "MG", "CA", "ZN", "SO4", "PO4", "GOL", "EDO", "ACT"}

def read_complex(pdb_path, ligand_resnames=None, protein_chains=None,
                 ligand_chains=None, skip_ions=True, model=1):
    """Read protein CA atoms and ligand heavy atoms from a PDB file."""
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
            if line[16] not in (" ", "A"):
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

            else:
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
        warnings.warn(f"{pdb_path}: no protein CA atoms selected.")
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
# 2. Laplacian & Schur Reduction
# ======================================================================

def build_laplacian(coords, rc=7.8, d0=1.0):
    """Distance-dependent Laplacian of a point set."""
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


def build_effective_laplacian(pdb_path, rc=7.8, d0=1.0, model=1,
                              protein_chains=("A",), ligand_chains=("B",),
                              verbose=False):
    """Build effective Laplacian for the protein."""
    p_xyz, p_ids, l_xyz, l_ids = read_complex(
        pdb_path,
        protein_chains=protein_chains,
        ligand_chains=ligand_chains,
        model=model
    )

    n_p, n_l = len(p_xyz), len(l_xyz)
    if verbose:
        print(f"{pdb_path}: {n_p} residues (CA), {n_l} ligand heavy atoms")
        if n_l:
            seen, ligand_residues = set(), []
            for s in l_ids:
                key = ":".join(s.split(":")[:3])
                if key not in seen:
                    seen.add(key)
                    ligand_residues.append(key)
            print(f"  ligand residues recognized: {', '.join(ligand_residues)}")

    if n_l == 0:
        if verbose:
            print("  no ligand found: returning plain protein Laplacian")
        return build_laplacian(p_xyz, rc, d0), p_ids

    coords = np.vstack([p_xyz, l_xyz])
    L_full = build_laplacian(coords, rc, d0)
    L_eff = schur_reduce(L_full, n_p)
    L_eff = (L_eff + L_eff.T) / 2

    return L_eff, p_ids


# ======================================================================
# 3. Effective Distances & Higher-Order Cooperation Index
# ======================================================================

def _require_connected_laplacian(L):
    """Raise if L is not numerically a connected Laplacian."""
    L = np.asarray(L, float)
    if L.ndim != 2 or L.shape[0] != L.shape[1]:
        raise ValueError("L must be a square matrix")
    if L.shape[0] < 2:
        raise ValueError("L must contain at least two nodes")

    eig = np.linalg.eigvalsh((L + L.T) / 2)
    tol = 1e-12 * max(float(abs(eig[-1])), 1.0)
    if eig[1] <= tol:
        raise ValueError(
            "The protein network is disconnected at the chosen rc/d0. "
            "Effective resistance between components is not finite."
        )


def effective_distances(L):
    """Compute pairwise effective resistance R_ij using regularized inverse."""
    # FIX: check connectivity before inversion
    _require_connected_laplacian(L)

    N = L.shape[0]
    J = np.ones((N, N)) / N
    K = np.linalg.inv(L + J) - J
    diag = np.diag(K)
    R = diag[:, None] + diag[None, :] - 2 * K
    return np.clip(R, 0.0, None)


def path_overlap(R, i, a, b):
    """Inner product of communication directions relative to reference node i."""
    return 0.5 * (R[i, a] + R[i, b] - R[a, b])


def Gram_matrix(R, subset):
    """
    Gram matrix K_S of Eq. C.5, rooted at the reference residue
    i = subset[0]:  (K_S)_pq = (R_ip + R_iq - R_pq) / 2.
    """
    i, rest = subset[0], list(subset[1:])
    m = len(rest)
    K = np.empty((m, m), float)
    for p in range(m):
        for q in range(m):
            K[p, q] = path_overlap(R, i, rest[p], rest[q])
    return K


def higher_order_chi(R, subset):
    """
    Normalized higher-order cooperation index (Eq. C.23).

        chi_S^(i) = det K_S / prod_{s in rest} R_{i,s}

    subset[0] is the REFERENCE residue i and stays FIXED; the remaining
    entries are the other members of S. For the fourth-order scan of
    Figure 3 pass [i, j, k, x] with i = 330, j = 372, k = 400 fixed and
    x running over the sequence -- NOT [x, i, j, k].

    det K_S is permutation-symmetric in S, so only the denominator
    depends on which residue is the reference. Moving the reference
    onto the running residue therefore changes the reported profile.
    """
    i, rest = subset[0], list(subset[1:])
    K = Gram_matrix(R, subset)
    R_S = float(np.linalg.det(K))
    denom = float(np.prod([R[i, s] for s in rest]))

    if denom <= 0 or np.isnan(R_S):
        return np.nan

    chi = R_S / denom

    # Clip tiny numerical excursions outside [0,1]
    if np.isfinite(chi):
        if chi < -1e-8 or chi > 1 + 1e-8:
            warnings.warn(f"chi={chi:.3e} outside [0,1]; clipping.")
        chi = np.clip(chi, 0.0, 1.0)

    return chi


# ======================================================================
# 4. Interactive Input Processing
# ======================================================================

def get_inputs():
    """Prompt user for proteins, cutoff, d0, and chain IDs."""
    print("\n" + "="*60)
    print("FOURTH-ORDER COOPERATION INDEX (chi_a,b,c,i) CALCULATOR")
    print("="*60)

    proteins = []
    for i in [1, 2]:
        while True:
            pdb_file = input(f"\nEnter PDB filename for protein {i} (e.g., 5ghd.pdb): ").strip()
            if not pdb_file:
                print("Please enter a filename.")
                continue
            if not os.path.exists(pdb_file):
                print(f"File '{pdb_file}' not found in current directory.")
                print(f"Current directory: {os.getcwd()}")
                continue
            break
        proteins.append(pdb_file)

    print("\n" + "-"*60)
    print("CHAIN SELECTION")
    print("-"*60)
    prot_chains = input("Enter protein chain ID(s) (comma-separated, default 'A'): ").strip()
    prot_chains = tuple(c.strip() for c in prot_chains.split(",")) if prot_chains else ("A",)
    lig_chains = input("Enter ligand chain ID(s) (comma-separated, default 'B'): ").strip()
    lig_chains = tuple(c.strip() for c in lig_chains.split(",")) if lig_chains else ("B",)

    print("\n" + "-"*60)
    print("PARAMETERS")
    print("-"*60)

    # Cutoff distance
    while True:
        try:
            rc = float(input("Enter cutoff distance in Å (default 7.8): ").strip() or "7.8")
            if rc <= 0:
                print("Cutoff must be positive.")
                continue
            break
        except ValueError:
            print("Please enter a valid number.")

    # Temperature/length scale d0
    while True:
        try:
            d0 = float(input("Enter temperature/length scale d0 (default 1.0): ").strip() or "1.0")
            if d0 <= 0:
                print("d0 must be positive.")
                continue
            break
        except ValueError:
            print("Please enter a valid number.")

    return proteins[0], proteins[1], prot_chains, lig_chains, rc, d0


def prompt_reference_residues(common_res):  # a = anchor i, b = j, c = k
    """Prompt user for reference residues a, b, and c from available common set."""
    print("\n" + "-"*60)
    print("REFERENCE RESIDUES (a, b, and c)")
    print("-"*60)
    print(f"Available common residue range: {min(common_res)} to {max(common_res)} (Total: {len(common_res)})")

    ref_res = []
    for label in ['a', 'b', 'c']:
        while True:
            try:
                val = input(f"Enter reference residue {label}: ").strip()
                res_num = int(val)
                if res_num not in common_res:
                    print(f"Residue {res_num} is not present in both proteins. Please select a common residue.")
                    continue
                if res_num in ref_res:
                    print(f"Residue {label} must be distinct from previous reference residues ({ref_res}).")
                    continue
                ref_res.append(res_num)
                break
            except ValueError:
                print("Please enter a valid integer residue number.")

    return ref_res[0], ref_res[1], ref_res[2]


# ======================================================================
# 5. Main Processing
# ======================================================================

def main():
    # np.seterr left at defaults: numerical warnings should be visible.

    pdb1, pdb2, prot_chains, lig_chains, rc, d0 = get_inputs()

    # Build Laplacians
    try:
        print("\nBuilding Laplacian for Protein 1...")
        L1, ids1 = build_effective_laplacian(pdb1, rc=rc, d0=d0,
                                             protein_chains=prot_chains,
                                             ligand_chains=lig_chains,
                                             verbose=True)

        print("\nBuilding Laplacian for Protein 2...")
        L2, ids2 = build_effective_laplacian(pdb2, rc=rc, d0=d0,
                                             protein_chains=prot_chains,
                                             ligand_chains=lig_chains,
                                             verbose=True)
    except Exception as e:
        print(f"\nError building Laplacians: {e}")
        return

    # Parse Residues
    resnums1, resnames1 = parse_residue_numbers(ids1)
    resnums2, resnames2 = parse_residue_numbers(ids2)

    # Create mapping from residue number -> matrix index.
    # Duplicate numbers would silently keep only the last occurrence,
    # so refuse them instead.
    def _index_map(resnums, tag):
        seen = {}
        for idx, r in enumerate(resnums):
            r = int(r)
            if r in seen:
                raise ValueError(
                    f"{tag}: residue number {r} appears more than once in "
                    f"the selected chains. Residue matching would be "
                    f"ambiguous. Restrict the chain selection or renumber."
                )
            seen[r] = idx
        return seen

    res2idx1 = _index_map(resnums1, "Protein 1")
    res2idx2 = _index_map(resnums2, "Protein 2")

    common_res = sorted(set(res2idx1) & set(res2idx2))

    if not common_res:
        print("\nError: No common residue indices found between the two proteins.")
        return

    print(f"\nResidue count comparison:")
    print(f"  Protein 1 total residues: {len(resnums1)}")
    print(f"  Protein 2 total residues: {len(resnums2)}")
    print(f"  Common residues count:    {len(common_res)}")

    # Matching by number alone would happily align two structures with
    # offset numbering. Report where the residue names disagree.
    mismatches = [(r, resnames1[res2idx1[r]], resnames2[res2idx2[r]])
                  for r in common_res
                  if resnames1[res2idx1[r]] != resnames2[res2idx2[r]]]
    if mismatches:
        print(f"  residue-name mismatches at {len(mismatches)} positions "
              f"(expected only at mutated sites):")
        for r, n1, n2 in mismatches[:10]:
            print(f"    {r}: {n1} -> {n2}")
        if len(mismatches) > 10:
            print(f"    ... and {len(mismatches) - 10} more")
        if len(mismatches) > 0.1 * len(common_res):
            warnings.warn("More than 10% of common positions disagree in "
                          "residue name. The two structures are probably "
                          "numbered differently.")

    # Prompt for reference residues a, b, and c
    res_a, res_b, res_c = prompt_reference_residues(common_res)

    print("\nCalculating effective resistances...")
    R1 = effective_distances(L1)
    R2 = effective_distances(L2)

    # Fixed reference quartet members. res_a is the anchor i.
    a1, b1, c1 = res2idx1[res_a], res2idx1[res_b], res2idx1[res_c]
    a2, b2, c2 = res2idx2[res_a], res2idx2[res_b], res2idx2[res_c]

    # Output filename setup
    base1 = os.path.splitext(os.path.basename(pdb1))[0]
    base2 = os.path.splitext(os.path.basename(pdb2))[0]
    output_file = f"{base1}_vs_{base2}_chi_i{res_a}_j{res_b}_k{res_c}.txt"

    print(f"\nCalculating chi^({res_a})_({res_b},{res_c},x) "
          f"and writing results to: {output_file}")
    with open(output_file, 'w') as f:
        f.write("# Fourth-order cooperation index, Laplacian minor hierarchy\n")
        f.write(f"# chi^(i)_(j,k,x) = R_ijkx / (R_ij R_ik R_ix)\n")
        f.write(f"# i = {res_a} (reference, FIXED), j = {res_b}, "
                f"k = {res_c}, x = running residue\n")
        f.write("# Anchored at i: the reference does not move with x.\n")
        f.write(f"# Protein 1: {pdb1}\n")
        f.write(f"# Protein 2: {pdb2}\n")
        f.write(f"# Protein chains: {','.join(prot_chains)}   Ligand chains: {','.join(lig_chains)}\n")
        f.write(f"# Parameters: rc = {rc} Å, d0 = {d0} Å\n")
        f.write(f"# Difference = Protein 2 - Protein 1\n")
        f.write(f"# x = {res_a} is undefined (R_ii = 0); "
                f"x = {res_b} and x = {res_c} are exactly 0.\n")
        f.write(f"# -------------------------------------------------------------------\n")
        f.write(f"# {'Res_x':<8} {'P1_chi':<18} {'P2_chi':<18} {'Diff_P2-P1':<18}\n")

        for r in common_res:
            x1 = res2idx1[r]
            x2 = res2idx2[r]

            # Reference residue i first: it stays fixed for every x.
            val1 = higher_order_chi(R1, [a1, b1, c1, x1])
            val2 = higher_order_chi(R2, [a2, b2, c2, x2])

            s1 = f"{val1:18.6f}" if np.isfinite(val1) else f"{'NaN':>18}"
            s2 = f"{val2:18.6f}" if np.isfinite(val2) else f"{'NaN':>18}"

            if not (np.isfinite(val1) and np.isfinite(val2)):
                sd = f"{'NaN':>18}"
            else:
                sd = f"{(val2 - val1):18.6f}"

            f.write(f"  {r:<8d} {s1} {s2} {sd}\n")

    print("\n" + "="*60)
    print("CALCULATION COMPLETE")
    print("="*60)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nProcess canceled by user.")
        sys.exit(0)