#!/usr/bin/env python3

import os
import random
import argparse
from collections import Counter
from decimal import Decimal, localcontext

from scipy.io import mmread
from scipy.sparse import block_diag
import numpy as np

from registry_utils import (
    MATRICES_DIR,
    parse_registry,
    read_matrix_market_metadata,
    repo_relative,
    upsert_registry_entry,
)

TINY_MTX_DIR = os.path.join(MATRICES_DIR, "TinyOriginal")
GENERATED_DIR = os.path.join(MATRICES_DIR, "FromTinyToLarge")

# Ensure the default output directory exists.
os.makedirs(GENERATED_DIR, exist_ok=True)


def load_tiny_permanents():
    tiny_entries = parse_registry().get("TinyOriginal", {})
    permanents = {}
    for relative_path, fields in tiny_entries.items():
        permanent_text = fields.get("Permanent")
        if permanent_text is None:
            continue
        permanent_value = Decimal(permanent_text)
        if permanent_value.is_zero():
            continue
        matrix_name = os.path.splitext(os.path.basename(relative_path))[0]
        permanents[matrix_name] = permanent_text
    return permanents

def load_tiny_matrices():
    matrices = {}
    tiny_permanents = load_tiny_permanents()
    tiny_files = [f for f in os.listdir(TINY_MTX_DIR) if f.endswith('.mtx')]
    for fname in tiny_files:
        name = fname[:-4]  # remove .mtx
        permanent_text = tiny_permanents.get(name)
        if permanent_text is None:
            continue
        mtx_path = os.path.join(TINY_MTX_DIR, fname)

        # Load matrix
        try:
            mat = mmread(mtx_path).tocsr().astype(np.int64)
            size = mat.shape[0]
            matrices[name] = {
                'size': size,
                'perm_text': permanent_text,
                'matrix': mat,
            }
        except Exception as e:
            print(f"Error loading {fname}: {e}")
    
    return matrices

def multiply_permanents(permanent_texts):
    decimals = [Decimal(text) for text in permanent_texts]
    precision = max(50, sum(len(value.as_tuple().digits) for value in decimals) + 10)
    with localcontext() as ctx:
        ctx.prec = precision
        product = Decimal("1")
        for value in decimals:
            product *= value
    product_text = format(product, "f")
    if "." in product_text:
        product_text = product_text.rstrip("0").rstrip(".")
    return product_text


def generate_block_matrix(target_N, tiny_matrices, seed=None):
    rng_py = random.Random(seed)
    rng_np = np.random.default_rng(seed)

    names = list(tiny_matrices.keys())
    selected = []
    total_size = 0
    
    while total_size < target_N:
        remaining = target_N - total_size
        feasible_names = [name for name in names if tiny_matrices[name]['size'] <= remaining]
        if not feasible_names:
            break
        name = rng_py.choice(feasible_names)
        info = tiny_matrices[name]
        selected.append(name)
        total_size += info['size']

    if not selected:
        raise ValueError(
            f"Could not assemble any TinyOriginal blocks within target size {target_N}"
        )
    
    # Build block diagonal without padding
    blocks = [tiny_matrices[name]['matrix'] for name in selected]
    large_matrix = block_diag(blocks, format='csr')
    
    # Compute permanent product
    perm_product = multiply_permanents(
        [tiny_matrices[name]['perm_text'] for name in selected]
    )
    
    # Now, permute rows and columns
    N = large_matrix.shape[0]
    row_perm = rng_np.permutation(N)
    col_perm = rng_np.permutation(N)
    large_matrix = large_matrix[row_perm, :][:, col_perm]
    
    return large_matrix, perm_product, selected

def save_matrix(matrix, filename):
    from scipy.io import mmwrite
    mmwrite(filename, matrix)

def update_registry(matrix_path, target_n, permanent, selected):
    metadata = read_matrix_market_metadata(matrix_path)
    counts = Counter(selected)
    composed_from = ", ".join(
        f"{name} x{counts[name]}" for name in sorted(counts)
    )
    upsert_registry_entry(
        "FromTinyToLarge",
        repo_relative(matrix_path),
        {
            "Permanent": permanent,
            **metadata,
            "Matrix type": "block diagonal from TinyOriginal",
            "Parameters": f"target_n={target_n}",
            "Composed from": composed_from,
            "Construction": "Block diagonal composition of TinyOriginal matrices with independent row and column permutations.",
        },
    )

def main():
    parser = argparse.ArgumentParser(description="Generate larger matrices from tiny ones with known permanents")
    parser.add_argument("target_N", type=int, help="Target size of the matrix (e.g., 100, 250)")
    parser.add_argument("--seed", type=int, default=None, help="Random seed")
    parser.add_argument("--num", type=int, default=1, help="Number of matrices to generate")
    parser.add_argument(
        "--output-dir",
        default=GENERATED_DIR,
        help="Directory to write generated matrices",
    )
    
    args = parser.parse_args()
    
    tiny_matrices = load_tiny_matrices()
    if not tiny_matrices:
        print("No tiny matrices loaded")
        return
    
    for i in range(args.num):
        iteration_seed = None if args.seed is None else args.seed + i
        matrix, perm, selected = generate_block_matrix(
            args.target_N,
            tiny_matrices,
            seed=iteration_seed,
        )
        actual_N = matrix.shape[0]
        matrix_name = f"generated_{args.target_N}_{actual_N}_{i+1}.mtx"
        output_dir = os.path.abspath(args.output_dir)
        os.makedirs(output_dir, exist_ok=True)
        filepath = os.path.join(output_dir, matrix_name)
        save_matrix(matrix, filepath)
        update_registry(filepath, args.target_N, perm, selected)
        print(f"Generated {matrix_name} with size {actual_N}x{actual_N} and permanent {perm}")

if __name__ == "__main__":
    main()
