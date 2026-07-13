#!/usr/bin/env python3

from __future__ import annotations

import argparse
import os

import numpy as np
from scipy.io import mmwrite
from scipy.sparse import coo_matrix

from band_permanent_computer import exact_values
from registry_utils import MATRICES_DIR, read_matrix_market_metadata, repo_relative, upsert_registry_entry

OUTPUT_DIR = os.path.join(MATRICES_DIR, "Banded")


def generate_band_matrix(
    n: int,
    k: int,
    weighted: bool = False,
    weight_low: float = 1.0,
    weight_high: float = 2.0,
    seed: int | None = None,
) -> coo_matrix:
    perm_rng = np.random.default_rng(seed)
    weight_rng = np.random.default_rng(None if seed is None else seed + 1)
    rows = []
    cols = []
    data = []
    dtype = np.longdouble if weighted else int

    for i in range(n):
        positions = {(i + delta) % n for delta in range(-k, k + 1)}
        for j in sorted(positions):
            rows.append(i)
            cols.append(j)
            if weighted:
                data.append(dtype(weight_rng.uniform(weight_low, weight_high)))
            else:
                data.append(1)

    matrix = coo_matrix((data, (rows, cols)), shape=(n, n), dtype=dtype).tocsr()
    row_perm = perm_rng.permutation(n)
    col_perm = perm_rng.permutation(n)
    matrix = matrix[row_perm, :][:, col_perm].tocoo()
    return matrix


def exact_band_permanent(n: int, k: int) -> int:
    if n <= 0:
        raise ValueError("n must be > 0")
    if k < 0:
        raise ValueError("k must be >= 0")
    return int(exact_values(k, [n])[n])


def update_registry(matrix_path: str, permanent: int, n: int, k: int) -> None:
    metadata = read_matrix_market_metadata(matrix_path)
    upsert_registry_entry(
        "Banded",
        repo_relative(matrix_path),
        {
            "Permanent": permanent,
            **metadata,
            "Matrix type": "cyclic band",
            "Parameters": f"n={n}, k={k}",
            "Construction": "A_{n,k}(i,j)=1 iff min(|pi(i)-sigma(j)|, n-|pi(i)-sigma(j)|) <= k for permutations pi and sigma.",
        },
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate cyclic band matrices and record exact permanents"
    )
    parser.add_argument("n", type=int, help="Matrix dimension")
    parser.add_argument("k", type=int, help="Band width")
    parser.add_argument(
        "--weighted",
        action="store_true",
        help="Generate weighted band matrices with real-valued entries",
    )
    parser.add_argument(
        "--weight-low",
        type=float,
        default=1.0,
        help="Lower bound for weighted entries",
    )
    parser.add_argument(
        "--weight-high",
        type=float,
        default=2.0,
        help="Upper bound for weighted entries",
    )
    parser.add_argument(
        "--output-dir",
        default=OUTPUT_DIR,
        help="Directory to write generated matrices",
    )
    parser.add_argument(
        "--filename",
        default=None,
        help="Optional matrix filename (without path)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Unused. Present only for CLI compatibility.",
    )

    args = parser.parse_args()

    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)
    matrix = generate_band_matrix(
        args.n,
        args.k,
        weighted=args.weighted,
        weight_low=args.weight_low,
        weight_high=args.weight_high,
        seed=args.seed,
    )
    permanent = None if args.weighted else exact_band_permanent(args.n, args.k)

    filename = args.filename or (
        f"band_n{args.n}_k{args.k}_weighted.mtx" if args.weighted else f"band_n{args.n}_k{args.k}.mtx"
    )
    filepath = os.path.abspath(os.path.join(output_dir, filename))
    mmwrite(filepath, matrix)

    if not args.weighted:
        update_registry(filepath, permanent, args.n, args.k)
    print(f"Generated {filepath}")
    if args.weighted:
        print(f"Weighted: True")
        print(f"Weight range: [{args.weight_low}, {args.weight_high}]")
    else:
        print(f"Permanent: {permanent}")


if __name__ == "__main__":
    main()
