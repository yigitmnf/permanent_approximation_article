from __future__ import annotations

from collections import defaultdict
from typing import Dict, Iterable, List
import numpy as np


def build_transitions(k: int) -> List[List[int]]:
    """Return the state graph for the cyclic bandwidth-k permanent problem.

    A state is a pair (L, R) of k-bit masks:
      - L marks unmatched columns on the left side of the cut
      - R marks already-matched columns on the right side of the cut

    The number of perfect matchings / permanent is trace(T^n), where T is the
    adjacency matrix of this state graph.
    """
    m = 1 << (2 * k)
    mask = (1 << k) - 1
    trans: List[List[int]] = [[] for _ in range(m)]

    for s in range(m):
        L = s & mask
        R = s >> k
        Lbits = [(L >> t) & 1 for t in range(k)]
        Rbits = [(R >> t) & 1 for t in range(k)]

        for delta in range(-k, k + 1):
            if Lbits[0] == 1 and delta != -k:
                continue
            if delta < 0:
                t = delta + k
                if Lbits[t] != 1:
                    continue
            elif delta < k:
                t = delta
                if Rbits[t] != 0:
                    continue

            Lp = [0] * k
            Rp = [0] * k

            for t in range(k - 1):
                bit = Lbits[t + 1]
                if delta == -k + (t + 1):
                    bit = 0
                Lp[t] = bit
            Lp[k - 1] = 1 if (Rbits[0] == 0 and delta != 0) else 0

            for t in range(k - 1):
                bit = Rbits[t + 1]
                if delta == t + 1:
                    bit = 1
                Rp[t] = bit
            Rp[k - 1] = 1 if delta == k else 0

            ns = sum(b << t for t, b in enumerate(Lp)) | (
                sum(b << t for t, b in enumerate(Rp)) << k
            )
            trans[s].append(ns)
    return trans


def build_blocks(k: int) -> Dict[int, List[List[int]]]:
    """Split states by invariant popcount(L)-popcount(R) and build predecessor lists."""
    trans = build_transitions(k)
    mask = (1 << k) - 1
    blocks: Dict[int, List[int]] = defaultdict(list)

    for s in range(len(trans)):
        d = (s & mask).bit_count() - (s >> k).bit_count()
        blocks[d].append(s)

    out: Dict[int, List[List[int]]] = {}
    for d, states in blocks.items():
        idx = {s: i for i, s in enumerate(states)}
        preds: List[List[int]] = [[] for _ in range(len(states))]
        for s in states:
            i = idx[s]
            for t in trans[s]:
                j = idx[t]
                preds[j].append(i)
        out[d] = preds
    return out


def exact_values(k: int, targets: Iterable[int]) -> Dict[int, int]:
    """Exact integer permanents for the cyclic band matrix of width k.

    This is practical for k=3,4,5 in the n-range requested by the user.
    It can also be run for k=6, but will be substantially slower.
    """
    targets = sorted(set(targets))
    if not targets:
        return {}
    N = targets[-1]
    target_set = set(targets)

    blocks = build_blocks(k)
    totals = {t: 0 for t in targets}

    # Symmetry: block d and block -d have the same trace sequence.
    ds = sorted(blocks)
    used = set()
    for d in ds:
        if d in used:
            continue
        preds = blocks[d]
        mult = 1 if d == 0 else 2
        used.add(d)
        used.add(-d)

        r = len(preds)
        P = np.eye(r, dtype=object)

        for n in range(1, N + 1):
            new = np.empty_like(P)
            for j, pr in enumerate(preds):
                if len(pr) == 1:
                    new[:, j] = P[:, pr[0]]
                else:
                    col = np.zeros(r, dtype=object)
                    for idx in pr:
                        col += P[:, idx]
                    new[:, j] = col
            P = new
            if n in target_set:
                totals[n] += mult * sum(P[i, i] for i in range(r))
    return totals


if __name__ == "__main__":
    targets = [500, 600, 700, 800, 900, 1000]
    for k in [3, 4, 5]:
        vals = exact_values(k, targets)
        print(f"k={k}")
        for n in targets:
            print(n, vals[n])
        print()
