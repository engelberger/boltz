#!/usr/bin/env python3
"""Align Boltz prediction structures using the Kabsch algorithm.

Given a predictions directory produced by Boltz (containing mmcif files), this
script reads all predicted structures, aligns them to the first structure in
Cα space and writes aligned copies into `<pred_dir>/aligned`.

Dependencies: gemmi, numpy
"""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import gemmi as gm


def parse_args():
    parser = argparse.ArgumentParser(description="Align Boltz predictions via Kabsch")
    parser.add_argument("pred_dir", type=Path, help="Path to Boltz predictions directory")
    parser.add_argument("--glob", default="*.cif", help="Glob pattern for CIF files inside pred_dir")
    parser.add_argument("--out", type=Path, default=None, help="Output directory (default: pred_dir/aligned)")
    return parser.parse_args()


def get_ca_coords(struct: gm.Structure) -> np.ndarray:
    """Return array of CA coordinates (N,3) for the first model."""
    ca_coords = []
    for model in struct:
        for chain in model:
            for res in chain:
                for atom in res:
                    if atom.name == "CA":
                        ca_coords.append(atom.pos.tolist())
        break  # only first model
    return np.asarray(ca_coords, dtype=np.float32)


def kabsch(P: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Compute optimal rotation matrix (Kabsch)."""
    C = P.T @ Q
    V, _, Wt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1.0, 1.0, d])
    return (V @ D @ Wt).astype(np.float32, copy=False)


def align_structure(mobile: gm.Structure, ref_ca: np.ndarray) -> gm.Structure:
    mob_ca = get_ca_coords(mobile)
    if mob_ca.shape != ref_ca.shape:
        raise ValueError("Residue count mismatch between structures")

    # centre
    mob_centroid = mob_ca.mean(axis=0)
    ref_centroid = ref_ca.mean(axis=0)
    P = mob_ca - mob_centroid
    Q = ref_ca - ref_centroid
    R = kabsch(P, Q)

    # apply rotation + translation to all atoms
    for model in mobile:
        for chain in model:
            for res in chain:
                for atom in res:
                    pos = atom.pos.tolist()
                    new = (np.dot((pos - mob_centroid), R) + ref_centroid).astype(np.float32)
                    atom.pos = gm.Position(*new)
        break
    return mobile


def main():
    args = parse_args()
    pred_dir: Path = args.pred_dir
    cif_files = sorted(pred_dir.glob(args.glob))
    if not cif_files:
        raise RuntimeError(f"No CIF files found in {pred_dir}")

    out_dir = args.out or pred_dir / "aligned"
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Found {len(cif_files)} structures. Using first as reference.")
    ref_struct = gm.read_structure(str(cif_files[0]))
    ref_ca = get_ca_coords(ref_struct)

    # Save reference copy as is
    ref_struct.write_mmCIF(str(out_dir / cif_files[0].name))

    # Align remaining
    for cif_path in cif_files[1:]:
        struct = gm.read_structure(str(cif_path))
        struct = align_structure(struct, ref_ca)
        struct.write_mmCIF(str(out_dir / cif_path.name))
        print(f"Aligned {cif_path.name} -> {out_dir/cif_path.name}")

    print("Done. Aligned files written to", out_dir)


if __name__ == "__main__":
    main() 