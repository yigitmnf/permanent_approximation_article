#!/usr/bin/env python3
"""
Block Diagonal Matrix Generator.

Generates n×n matrices composed of blocks where k is chosen uniformly
from {5, 6, 7, 8, 9}. Block sizes are sampled randomly and the whole draw is
restarted until the sizes sum exactly to n. The diagonal is always preserved.
Off-diagonal entries can be removed independently with probability p =
reduction_ratio, and the final matrix is permuted by independent row and column
permutations.
"""

import numpy as np
from scipy.sparse import block_diag, coo_matrix
from scipy.io import mmwrite
import argparse
import os
from collections import Counter

from registry_utils import MATRICES_DIR, read_matrix_market_metadata, repo_relative, upsert_registry_entry

BLOCK_SIZES = (5, 6, 7, 8, 9)

def format_ratio(value):
    if float(value).is_integer():
        return str(int(value))
    text = f"{value}".rstrip("0").rstrip(".")
    return text


def ratio_tag(value):
    return format_ratio(value).replace(".", "p")


def exact_block_permanent(block):
    """Compute the permanent of a small dense block via Ryser's formula."""
    n = block.shape[0]
    if n == 0:
        return 1

    row_sums = [np.longdouble(0)] * n
    total = np.longdouble(0)
    previous_gray = 0

    for subset_index in range(1, 1 << n):
        gray = subset_index ^ (subset_index >> 1)
        changed = gray ^ previous_gray
        changed_bit = changed.bit_length() - 1
        delta = np.longdouble(1 if (gray & changed) else -1)

        for row in range(n):
            row_sums[row] += delta * np.longdouble(block[row, changed_bit])

        product = np.longdouble(1)
        for row_sum in row_sums:
            product *= row_sum
            if product == 0:
                break

        if (gray.bit_count() - n) % 2:
            total -= product
        else:
            total += product

        previous_gray = gray

    return total


def generate_block(k, reduction_ratio, rng, weighted=False, weight_low=1.0, weight_high=2.0):
    """Generate one block with preserved diagonal and optional off-diagonal thinning."""
    dtype = np.longdouble if weighted else int
    if weighted:
        block = rng.uniform(weight_low, weight_high, size=(k, k)).astype(np.longdouble)
    else:
        block = np.ones((k, k), dtype=dtype)
    if k > 1 and reduction_ratio > 0:
        off_diagonal = ~np.eye(k, dtype=bool)
        removal_mask = rng.random((k, k)) < reduction_ratio
        block[off_diagonal & removal_mask] = 0
    permanent = exact_block_permanent(block)
    return block, permanent


def sample_block_sizes(n, rng, max_attempts=100000):
    """
    Sample block sizes by repeated independent draws from BLOCK_SIZES.
    Restart until the draw sums exactly to n.
    """
    for _ in range(max_attempts):
        total = 0
        block_sizes = []
        while total < n:
            k = int(rng.choice(BLOCK_SIZES))
            block_sizes.append(k)
            total += k
        if total == n:
            return block_sizes
    raise ValueError(
        f"Could not sample block sizes summing to n={n} after {max_attempts} attempts"
    )


def generate_block_diagonal_matrix(
    n,
    reduction_ratio=0.0,
    seed=None,
    weighted=False,
    weight_low=1.0,
    weight_high=2.0,
):
    """
    Generate a block diagonal matrix with blocks where k ∈ {5, 6, 7, 8, 9}.
    The block diagonal is preserved and off-diagonal entries are removed with
    probability `reduction_ratio`. In weighted mode, block entries are sampled
    uniformly from [weight_low, weight_high].

    Args:
        n: matrix dimension (will be adjusted to fit blocks)
        reduction_ratio: probability of removing an off-diagonal entry
        seed: random seed for reproducibility
        weighted: if True, generate real-valued blocks
        weight_low: lower bound for weighted entries
        weight_high: upper bound for weighted entries

    Returns:
        Sparse COO matrix, list of block sizes, and list of exact block permanents
    """
    if reduction_ratio < 0 or reduction_ratio > 1:
        raise ValueError("reduction_ratio must be between 0 and 1")
    if weight_low > weight_high:
        raise ValueError("weight_low must be <= weight_high")

    size_rng = np.random.default_rng(seed)
    block_rng_seed = None if seed is None else seed + 1
    block_rng = np.random.default_rng(block_rng_seed)
    perm_rng_seed = None if seed is None else seed + 2
    perm_rng = np.random.default_rng(perm_rng_seed)

    blocks = []
    block_sizes = sample_block_sizes(n, size_rng)
    block_permanents = []

    for k in block_sizes:
        block, block_permanent = generate_block(
            k,
            reduction_ratio,
            block_rng,
            weighted=weighted,
            weight_low=weight_low,
            weight_high=weight_high,
        )
        blocks.append(block)
        block_permanents.append(block_permanent)

    # Create block diagonal matrix
    if blocks:
        matrix = block_diag(blocks, format='csr')
        matrix.eliminate_zeros()
        row_perm = perm_rng.permutation(n)
        col_perm = perm_rng.permutation(n)
        matrix = matrix[row_perm, :][:, col_perm]
        matrix = matrix.tocoo()
    else:
        # Fallback for very small n
        matrix = coo_matrix((n, n))

    return matrix, block_sizes, block_permanents


def calculate_permanent(block_permanents):
    """
    Calculate the permanent of a block diagonal matrix.
    For a block diagonal matrix, permanent = product of the block permanents.
    """
    permanent = np.longdouble(1)
    for block_permanent in block_permanents:
        permanent *= block_permanent
    return permanent


def update_registry(matrix_path, permanent, n, block_sizes, reduction_ratio, weighted, weight_low, weight_high):
    metadata = read_matrix_market_metadata(matrix_path)
    counts = Counter(block_sizes)
    count_summary = ", ".join(
        f"{size}:{counts[size]}" for size in sorted(counts)
    )
    ratio_text = format_ratio(reduction_ratio)
    weight_range_text = f"[{format_ratio(weight_low)}, {format_ratio(weight_high)}]"
    matrix_type = (
        ("full weighted block diagonal" if reduction_ratio == 0 else "reduced weighted block diagonal")
        if weighted
        else ("full block diagonal" if reduction_ratio == 0 else "reduced block diagonal")
    )
    construction = (
        (
            "Matrix is built block diagonally from weighted blocks with entries "
            f"sampled uniformly from {weight_range_text} and then independently "
            "permuted in rows and columns. The exact permanent is the product of "
            "the exact block permanents."
            if weighted
            else "Matrix is built block diagonally from the sampled blocks and then "
            "independently permuted in rows and columns. Permanent equals the product "
            "of factorials of the block sizes."
        )
        if reduction_ratio == 0
        else (
            "Diagonal entries are always kept. Off-diagonal entries inside each "
            f"block are removed independently with probability {ratio_text}. "
            "The block-diagonal matrix is then independently permuted in rows and "
            "columns. "
            + (
                "The exact permanent is the product of the exact block permanents."
                if weighted
                else "The exact permanent is the product of the exact block permanents."
            )
        )
    )
    upsert_registry_entry(
        "BlockDiagonal",
        repo_relative(matrix_path),
        {
            "Permanent": permanent,
            **metadata,
            "Matrix type": matrix_type,
            "Parameters": (
                f"n={n}, block_size_counts={{{count_summary}}}, "
                f"reduction_ratio={ratio_text}, weighted={weighted}"
                + (f", weight_range={weight_range_text}" if weighted else "")
            ),
            "Construction": construction,
        },
    )


def main():
    parser = argparse.ArgumentParser(
        description="Generate block diagonal matrices with block sizes in {5,6,7,8,9}"
    )
    parser.add_argument("n", type=int, help="Matrix dimension")
    parser.add_argument(
        "--output-dir",
        default=os.path.join(MATRICES_DIR, "BlockDiagonal"),
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
        "--reduction-ratio",
        type=float,
        default=0.0,
        help="Probability of removing each off-diagonal block entry",
    )
    parser.add_argument(
        "--weighted",
        action="store_true",
        help="Generate weighted blocks with real-valued entries",
    )
    parser.add_argument(
        "--weight-low",
        type=float,
        default=1.0,
        help="Lower bound for weighted block entries",
    )
    parser.add_argument(
        "--weight-high",
        type=float,
        default=2.0,
        help="Upper bound for weighted block entries",
    )

    args = parser.parse_args()

    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    if args.weight_low > args.weight_high:
        raise ValueError("weight-low must be <= weight-high")

    matrix, block_sizes, block_permanents = generate_block_diagonal_matrix(
        args.n,
        reduction_ratio=args.reduction_ratio,
        seed=args.seed,
        weighted=args.weighted,
        weight_low=args.weight_low,
        weight_high=args.weight_high,
    )
    permanent = calculate_permanent(block_permanents)

    if args.filename:
        filename = args.filename
    elif args.reduction_ratio == 0:
        filename = (
            f"block_diagonal_n{args.n}_weighted.mtx"
            if args.weighted
            else f"block_diagonal_n{args.n}.mtx"
        )
    else:
        filename = (
            f"block_diagonal_n{args.n}_p{ratio_tag(args.reduction_ratio)}_weighted.mtx"
            if args.weighted
            else f"block_diagonal_n{args.n}_p{ratio_tag(args.reduction_ratio)}.mtx"
        )
    filepath = os.path.abspath(os.path.join(output_dir, filename))
    mmwrite(filepath, matrix)
    try:
        update_registry(
            filepath,
            permanent,
            args.n,
            block_sizes,
            args.reduction_ratio,
            args.weighted,
            args.weight_low,
            args.weight_high,
        )
    except ValueError as exc:
        print(f"Skipping registry update: {exc}")

    # Print statistics
    print(f"Generated {filepath}")
    print(f"Dimension: {args.n}x{args.n}")
    print(f"Block sizes: {block_sizes}")
    print(f"Number of blocks: {len(block_sizes)}")
    print(f"Reduction ratio: {format_ratio(args.reduction_ratio)}")
    print(f"Weighted: {args.weighted}")
    if args.weighted:
        print(f"Weight range: [{format_ratio(args.weight_low)}, {format_ratio(args.weight_high)}]")
    print(f"Permanent: {permanent}")
    print(f"Total non-zeros: {matrix.nnz}")


if __name__ == "__main__":
    main()
