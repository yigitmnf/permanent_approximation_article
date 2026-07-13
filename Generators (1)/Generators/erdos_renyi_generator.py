#!/usr/bin/env python3
"""
Erdős-Rényi random matrix generator.

Generates random n×n matrices where each entry is non-zero with probability p = k/n,
resulting in expected k non-zeros per row/column.
"""

import numpy as np
from scipy.sparse import random as sparse_random
from scipy.sparse import coo_matrix
from scipy.io import mmwrite
import argparse
import os

from registry_utils import MATRICES_DIR


def generate_erdos_renyi_matrix(n, k, seed=None, weighted=False, weight_low=1.0, weight_high=2.0):
    """
    Generate Erdős-Rényi random matrix with expected k non-zeros per row.
    Matrix values are binary (1) by default and include a full diagonal for
    perfect matching. In weighted mode, the non-zero entries are sampled
    uniformly from [weight_low, weight_high].

    Args:
        n: matrix dimension
        k: expected non-zeros per row/column
        seed: random seed for reproducibility
        weighted: if True, generate real-valued weights
        weight_low: lower bound for weighted entries
        weight_high: upper bound for weighted entries

    Returns:
        Sparse COO matrix
    """
    p = k / n  # Probability of each entry being non-zero
    rng = np.random.default_rng(seed)

    if weighted:
        def data_rvs(nnz):
            return rng.uniform(weight_low, weight_high, size=nnz)
        matrix = sparse_random(
            n,
            n,
            density=p,
            format='coo',
            random_state=seed,
            data_rvs=data_rvs,
        )
    else:
        # Generate random sparse matrix with density p
        matrix = sparse_random(n, n, density=p, format='coo', random_state=seed)
        # Make values binary (1)
        matrix.data = np.ones_like(matrix.data, dtype=int)

    # Add full diagonal to ensure perfect matching
    # First, remove any existing diagonal entries to avoid duplicates
    off_diag_mask = matrix.row != matrix.col
    matrix.row = matrix.row[off_diag_mask]
    matrix.col = matrix.col[off_diag_mask]
    matrix.data = matrix.data[off_diag_mask]
    
    # Add diagonal entries
    diag_rows = np.arange(n)
    diag_cols = np.arange(n)
    if weighted:
        diag_data = rng.uniform(weight_low, weight_high, size=n)
    else:
        diag_data = np.ones(n, dtype=int)
    
    # Combine off-diagonal and diagonal
    all_rows = np.concatenate([matrix.row, diag_rows])
    all_cols = np.concatenate([matrix.col, diag_cols])
    all_data = np.concatenate([matrix.data, diag_data])

    matrix = coo_matrix((all_data, (all_rows, all_cols)), shape=(n, n))
    
    return matrix


def main():
    parser = argparse.ArgumentParser(
        description="Generate Erdős-Rényi random sparse matrices"
    )
    parser.add_argument("n", type=int, help="Matrix dimension")
    parser.add_argument(
        "k", type=int, help="Expected non-zeros per row/column"
    )
    parser.add_argument(
        "--output-dir",
        default=os.path.join(MATRICES_DIR, "ErdosRenyi"),
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
        help="Random seed for reproducibility",
    )
    parser.add_argument(
        "--weighted",
        action="store_true",
        help="Generate real-valued entries instead of binary 0/1 entries",
    )
    parser.add_argument(
        "--weight-low",
        type=float,
        default=0.0,
        help="Lower bound for weighted entries",
    )
    parser.add_argument(
        "--weight-high",
        type=float,
        default=10.0,
        help="Upper bound for weighted entries",
    )

    args = parser.parse_args()

    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    if args.weight_low > args.weight_high:
        raise ValueError("weight-low must be <= weight-high")

    matrix = generate_erdos_renyi_matrix(
        args.n,
        args.k,
        seed=args.seed,
        weighted=args.weighted,
        weight_low=args.weight_low,
        weight_high=args.weight_high,
    )

    if args.filename:
        filename = args.filename
    elif args.weighted:
        filename = f"erdos_renyi_n{args.n}_k{args.k}_weighted.mtx"
    else:
        filename = f"erdos_renyi_n{args.n}_k{args.k}.mtx"
    filepath = os.path.abspath(os.path.join(output_dir, filename))
    mmwrite(filepath, matrix)

    # Print statistics
    actual_nnz = matrix.nnz
    actual_nnz_per_row = actual_nnz / args.n if args.n > 0 else 0
    probability = args.k / args.n

    print(f"Generated {filepath}")
    print(f"Dimension: {args.n}x{args.n}")
    print(f"Expected nnz per row: {args.k}")
    print(f"Probability p = k/n: {probability:.6f}")
    print(f"Weighted: {args.weighted}")
    if args.weighted:
        print(f"Weight range: [{args.weight_low}, {args.weight_high}]")
    print(f"Total non-zeros: {actual_nnz}")
    print(f"Actual nnz per row: {actual_nnz_per_row:.2f}")


if __name__ == "__main__":
    main()
