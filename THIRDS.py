#!/usr/bin/env python3
"""
THIRDS_ligand_patched.py
========================

Third-order cooperation index of the Laplacian minor hierarchy,

    chi^(i)_{j,k} = (K^(i)_jk)^2 / (R_ij R_ik)

with K^(i)_jk = (R_ij + R_ik - R_jk) / 2   (Eq. B.8 / C.5 / C.19).

Convention:

    chi = 0  uncorrelated
    chi = 1  correlated

Residue i is the REFERENCE (anchor). The normalizer contains only
distances rooted at i, so chi is deliberately NOT symmetric under
i <-> j; the exchange ratio is chi^(i)_{j,k} / chi^(j)_{i,k} = R_jk / R_ik.
This follows Eq. C.23, where the index is written chi_S^(i) precisely
because the denominator depends on the reference choice.

Compares two structures and reports the difference (Protein 2 - Protein 1)
over residues present in both.

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
            # FIX 4b: PDB hydrogen names may start with a digit ("1HB", "2HG1").
            # Strip leading digits before the fallback element test.
            bare = name.lstrip("0123456789")
            is_hydrogen = (element in ("H", "D")
                           or (not element and bare[:1] in ("H", "D")))

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
    """Convert labels into integer numbers and residue names."""
    nums, names = [], []
    for s in ids:
        field = s.split(":")[1]
        digits = "".join(c for c in field if c.isdigit() or c == "-")
        nums.append(int(digits))
        names.append(s.split(":")[2])
    return np.asarray(nums, int), np.asarray(names)


# ======================================================================
# 2. Laplacian & Schur Reduction
# ======================================================================

def build_laplacian(coords, rc=7.8, d0=1.0):
    """Distance-dependent Laplacian of a point set (Eq. 1)."""
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
    """Kron reduction: eliminate the trailing (ligand) block of L."""
    Hpp = L[:n_protein, :n_protein]
    Hpl = L[:n_protein, n_protein:]
    Hlp = L[n_protein:, :n_protein]
    Hll = L[n_protein:, n_protein:]

    if Hll.shape[0] == 0:
        return Hpp.copy()

    cond = np.linalg.cond(Hll)
    if not np.isfinite(cond) or cond > cond_warn:
        warnings.warn(f"H_ll is ill-conditioned (cond = {cond:.2e}). "
                      f"Falling back to least-squares solve. Check for "
                      f"ligand atoms isolated at the chosen rc.")
        X = np.linalg.lstsq(Hll, Hlp, rcond=None)[0]
    else:
        X = np.linalg.solve(Hll, Hlp)

    return Hpp - Hpl @ X


def build_effective_laplacian(pdb_path, ligand_resnames=None,
                              protein_chains=("A",), ligand_chains=("B",),
                              rc=7.8, d0=1.0, model=1, verbose=False):
    """Effective protein Laplacian H_pp - H_pl H_ll^-1 H_lp.

    FIX 4: chain selection is now an argument, and the ligand actually
    recognized is printed, so a wrong chain choice is visible.
    """
    p_xyz, p_ids, l_xyz, l_ids = read_complex(
        pdb_path,
        ligand_resnames=ligand_resnames,
        protein_chains=list(protein_chains) if protein_chains else None,
        ligand_chains=list(ligand_chains) if ligand_chains else None,
        model=model,
    )

    n_p, n_l = len(p_xyz), len(l_xyz)
    if verbose:
        print(f"{pdb_path}: {n_p} residues (CA), {n_l} ligand heavy atoms")
        print(f"  protein chains: {list(protein_chains)}   "
              f"ligand chains: {list(ligand_chains)}")
        if n_l:
            seen, ligand_residues = set(), []
            for s in l_ids:
                key = ":".join(s.split(":")[:3])
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

    coords = np.vstack([p_xyz, l_xyz])
    L_full = build_laplacian(coords, rc, d0)
    L_eff = schur_reduce(L_full, n_p)
    L_eff = (L_eff + L_eff.T) / 2

    if verbose:
        off_eff = np.abs(L_eff - np.diag(np.diag(L_eff))) > 1e-10
        direct = L_full[:n_p, :n_p]
        off_dir = np.abs(direct - np.diag(np.diag(direct))) > 1e-10
        print(f"  ligand-mediated edges created: "
              f"{int((off_eff & ~off_dir).sum() // 2)}")

    return L_eff, p_ids


# ======================================================================
# 3. Effective Distances & Cooperation Index
# ======================================================================

def _require_connected_laplacian(L):
    """FIX 2: raise if L is not numerically a connected Laplacian."""
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
    """Pairwise effective resistance R_ij = K_ii + K_jj - 2 K_ij (Eq. 5)."""
    _require_connected_laplacian(L)
    N = L.shape[0]
    J = np.ones((N, N)) / N
    K = np.linalg.solve(L + J, np.eye(N)) - J

    diag = np.diag(K)
    R = diag[:, None] + diag[None, :] - 2 * K

    # FIX 2b: report, rather than silently absorb, a genuine failure.
    worst = float(R.min())
    if worst < -1e-8 * float(R.max()):
        warnings.warn(f"Effective resistance came out negative "
                      f"(min = {worst:.3e}); the inversion is unreliable.")
    return np.clip(R, 0.0, None)


def cooperation_index_vector(R, ref_idx, j_idx):
    """chi^(i)_{j,k} for fixed reference i = ref_idx and partner j = j_idx,
    evaluated for every running residue k.

        K^(i)_jk = (R_ij + R_ik - R_jk) / 2
        chi      = (K^(i)_jk)^2 / (R_ij R_ik)

    chi = 0 is uncorrelated, chi = 1 is correlated.

    FIX 3: only k = i is undefined (R_ii = 0 in the denominator).
    At k = j the value is exactly 1 - the two displacement vectors
    coincide, so the overlap is complete - and it is returned.
    """
    if ref_idx == j_idx:
        raise ValueError("reference residue i and partner residue j "
                         "must be distinct")

    R_ij = R[ref_idx, j_idx]
    R_ik = R[ref_idx, :]
    R_jk = R[j_idx, :]

    K_ijk = 0.5 * (R_ij + R_ik - R_jk)

    with np.errstate(divide="ignore", invalid="ignore"):
        chi = (K_ijk ** 2) / (R_ij * R_ik)

    chi = np.asarray(chi, float)
    chi[ref_idx] = np.nan          # R_ii = 0: genuinely undefined
    chi[j_idx] = 1.0               # exact value, not a singularity

    # Gram/Hadamard bounds (Eq. C.28). Clip only float dust.
    if chi.size:
        finite = np.isfinite(chi)
        stray = finite & ((chi < -1e-8) | (chi > 1 + 1e-8))
        if stray.any():
            warnings.warn(f"{int(stray.sum())} chi values fell outside [0,1] "
                          f"beyond tolerance; check network conditioning.")
        chi[finite] = np.clip(chi[finite], 0.0, 1.0)
    return chi


# ======================================================================
# 4. Interactive Input
# ======================================================================

def _prompt_pdb_file(label):
    while True:
        pdb_file = input(f"\nEnter PDB filename for {label} "
                         f"(e.g., 5HED.pdb): ").strip()
        if not pdb_file:
            print("Please enter a filename.")
            continue
        if not os.path.exists(pdb_file):
            print(f"File '{pdb_file}' not found.")
            print(f"Current directory: {os.getcwd()}")
            continue
        return pdb_file


def _prompt_positive(prompt, default):
    while True:
        try:
            val = float(input(prompt).strip() or str(default))
            if val <= 0:
                print("Value must be positive.")
                continue
            return val
        except ValueError:
            print("Please enter a valid number.")


def _prompt_chains(role, default):
    raw = input(f"Chain ID(s) for {role} "
                f"[default {','.join(default) or 'none'}]: ").strip()
    if not raw:
        return default
    if raw.lower() in ("none", "-"):
        return ()
    return tuple(c.strip() for c in raw.split(",") if c.strip())


def get_inputs():
    print("\n" + "=" * 62)
    print("THIRD-ORDER COOPERATION INDEX  chi^(i)_{j,k}")
    print("=" * 62)

    pdb1 = _prompt_pdb_file("Protein 1")
    pdb2 = _prompt_pdb_file("Protein 2")

    print("\n" + "-" * 62)
    print("CHAIN SELECTION")
    print("-" * 62)
    protein_chains = _prompt_chains("the protein", ("A",))
    ligand_chains = _prompt_chains("the peptide ligand", ("B",))

    print("\n" + "-" * 62)
    print("PARAMETERS")
    print("-" * 62)
    rc = _prompt_positive("Enter cutoff distance in A (default 7.8): ", 7.8)
    d0 = _prompt_positive("Enter effective temperature d0 (default 1.0): ", 1.0)

    return pdb1, pdb2, protein_chains, ligand_chains, rc, d0


def prompt_reference_residues(common_res):
    """Reference residue i (the anchor) and partner residue j."""
    print("\n" + "-" * 62)
    print("REFERENCE RESIDUE i (anchor) AND PARTNER RESIDUE j")
    print("-" * 62)
    print(f"Common residues: {min(common_res)} to {max(common_res)} "
          f"(total {len(common_res)})")
    print("Note: chi is anchored at i. Swapping i and j changes the result.")

    chosen = []
    for label, role in (("i", "reference/anchor"), ("j", "partner")):
        while True:
            try:
                num = int(input(f"Enter residue {label} ({role}): ").strip())
            except ValueError:
                print("Please enter an integer residue number.")
                continue
            if num not in common_res:
                print(f"Residue {num} is not present in both structures.")
                continue
            if chosen and num == chosen[0]:
                print(f"Residue j must differ from residue i ({chosen[0]}).")
                continue
            chosen.append(num)
            break
    return chosen[0], chosen[1]


# ======================================================================
# 5. Main
# ======================================================================

def _index_map(resnums, tag):
    """FIX 1: refuse duplicate residue numbers instead of silently
    keeping the last occurrence."""
    seen = {}
    for idx, r in enumerate(resnums):
        r = int(r)
        if r in seen:
            raise ValueError(
                f"{tag}: residue number {r} appears more than once in the "
                f"selected chains. Residue matching would be ambiguous. "
                f"Restrict the chain selection or renumber the structure."
            )
        seen[r] = idx
    return seen


def main():
    pdb1, pdb2, protein_chains, ligand_chains, rc, d0 = get_inputs()

    print("\nBuilding Laplacian for Protein 1...")
    L1, ids1 = build_effective_laplacian(
        pdb1, protein_chains=protein_chains, ligand_chains=ligand_chains,
        rc=rc, d0=d0, verbose=True)

    print("\nBuilding Laplacian for Protein 2...")
    L2, ids2 = build_effective_laplacian(
        pdb2, protein_chains=protein_chains, ligand_chains=ligand_chains,
        rc=rc, d0=d0, verbose=True)

    resnums1, resnames1 = parse_residue_numbers(ids1)
    resnums2, resnames2 = parse_residue_numbers(ids2)

    idx1 = _index_map(resnums1, "Protein 1")
    idx2 = _index_map(resnums2, "Protein 2")

    common_res = sorted(set(idx1) & set(idx2))
    if not common_res:
        raise ValueError("The two structures share no residue numbers.")

    # FIX 5: matching by number alone would happily align two different
    # depositions with offset numbering. Report the mismatches.
    mismatches = [(r, resnames1[idx1[r]], resnames2[idx2[r]])
                  for r in common_res
                  if resnames1[idx1[r]] != resnames2[idx2[r]]]

    print(f"\nResidue counts: P1 = {len(resnums1)}, P2 = {len(resnums2)}, "
          f"common = {len(common_res)}")
    if mismatches:
        print(f"  residue-name mismatches at {len(mismatches)} common "
              f"positions (expected only at mutated sites):")
        for r, n1, n2 in mismatches[:10]:
            print(f"    {r}: {n1} -> {n2}")
        if len(mismatches) > 10:
            print(f"    ... and {len(mismatches) - 10} more")
        if len(mismatches) > 0.1 * len(common_res):
            warnings.warn("More than 10% of common positions disagree in "
                          "residue name. The two structures are probably "
                          "numbered differently.")

    res_i, res_j = prompt_reference_residues(common_res)

    print("\nCalculating effective resistances...")
    R1 = effective_distances(L1)
    R2 = effective_distances(L2)

    chi1 = cooperation_index_vector(R1, idx1[res_i], idx1[res_j])
    chi2 = cooperation_index_vector(R2, idx2[res_i], idx2[res_j])

    base1 = os.path.splitext(os.path.basename(pdb1))[0]
    base2 = os.path.splitext(os.path.basename(pdb2))[0]
    output_file = f"{base1}_vs_{base2}_chi_i{res_i}_j{res_j}.txt"

    with open(output_file, "w") as f:
        f.write("# Third-order cooperation index, Laplacian minor hierarchy\n")
        f.write(f"# chi^(i)_(j,k) = (K^(i)_jk)^2 / (R_ij R_ik), "
                f"i = {res_i} (reference), j = {res_j} (partner)\n")
        f.write("# chi = 0 uncorrelated, chi = 1 correlated.\n")
        f.write("# Anchored at i: chi is NOT symmetric under i <-> j.\n")
        f.write(f"# Protein 1: {pdb1}\n")
        f.write(f"# Protein 2: {pdb2}\n")
        f.write(f"# Protein chains: {','.join(protein_chains) or 'all'}   "
                f"Ligand chains: {','.join(ligand_chains) or 'none'}\n")
        f.write(f"# Parameters: rc = {rc} A, d0 = {d0}\n")
        f.write("# Difference = Protein 2 - Protein 1\n")
        f.write(f"# k = {res_i} is undefined (R_ii = 0); "
                f"k = {res_j} is exactly 1.\n")
        f.write("# ---------------------------------------------------------"
                "----------\n")
        f.write(f"# {'Res_k':<8} {'P1_chi':<18} {'P2_chi':<18} "
                f"{'Diff_P2-P1':<18}\n")

        for r in common_res:
            v1, v2 = chi1[idx1[r]], chi2[idx2[r]]
            s1 = f"{v1:18.6f}" if np.isfinite(v1) else f"{'NaN':>18}"
            s2 = f"{v2:18.6f}" if np.isfinite(v2) else f"{'NaN':>18}"
            sd = (f"{(v2 - v1):18.6f}"
                  if np.isfinite(v1) and np.isfinite(v2) else f"{'NaN':>18}")
            f.write(f"  {r:<8d} {s1} {s2} {sd}\n")

    print("\n" + "=" * 62)
    print("CALCULATION COMPLETE")
    print("=" * 62)
    print(f"Results saved to: {output_file}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nProcess canceled by user.")
        sys.exit(0)
    except Exception as exc:
        print(f"\nError: {exc}")
        sys.exit(1)
