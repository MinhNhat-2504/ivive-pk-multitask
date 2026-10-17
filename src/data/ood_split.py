"""Out-of-domain splits, and evidence that the benchmark's split is not one.

A scaffold split sorts molecules by Bemis-Murcko framework, which separates
test from train by *identity* of scaffold but not by *distance*: two different
scaffolds can be one atom apart. Whether such a split actually probes
extrapolation is an empirical question, answered here by measuring how similar
each test molecule is to its nearest training neighbour.

The alternative built here clusters molecules by fingerprint similarity and
holds out whole clusters, so a test molecule has no close relative in training
by construction. Both splits are then characterised by the same statistic, and
the models of Section 7 are compared across both.
"""

from __future__ import annotations

import numpy as np
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.ML.Cluster import Butina

RDLogger.DisableLog("rdApp.*")

_generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def fingerprints(smiles: list[str]) -> tuple[list, np.ndarray]:
    """Morgan bit vectors, plus a mask of which inputs parsed."""
    bits, ok = [], np.zeros(len(smiles), dtype=bool)
    for i, smi in enumerate(smiles):
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            bits.append(_generator.GetFingerprint(mol))
            ok[i] = True
    return bits, ok


def nearest_neighbour_similarity(query: list, reference: list) -> np.ndarray:
    """Max Tanimoto from each query molecule to anything in the reference set.

    This is the statistic that says whether a split demands extrapolation. A
    test molecule with a 0.9-similar neighbour in training is being
    interpolated to, whatever the scaffold bookkeeping says.
    """
    return np.array(
        [max(DataStructs.BulkTanimotoSimilarity(fp, reference)) for fp in query]
    )


def butina_clusters(bits: list, cutoff: float = 0.65) -> list[list[int]]:
    """Butina clustering on Tanimoto distance.

    `cutoff` is a distance, so 0.65 groups molecules sharing at least 0.35
    similarity -- loose enough that clusters correspond to chemical series
    rather than near-duplicates.
    """
    n = len(bits)
    distances: list[float] = []
    for i in range(1, n):
        similarities = DataStructs.BulkTanimotoSimilarity(bits[i], bits[:i])
        distances.extend(1.0 - s for s in similarities)
    clusters = Butina.ClusterData(distances, n, cutoff, isDistData=True)
    return [list(c) for c in clusters]


def cluster_split(
    smiles: list[str], test_size: float = 0.2, seed: int = 0, cutoff: float = 0.65
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Hold out whole clusters, so no test molecule has a training sibling.

    Returns boolean train and test masks over the input list, and the mask of
    molecules that parsed. Clusters are shuffled rather than taken largest- or
    smallest-first: taking the smallest would stack the test set with
    singletons and overstate the difficulty, taking the largest would hand it
    whole well-populated series and understate it.
    """
    bits, ok = fingerprints(smiles)
    index_of = np.flatnonzero(ok)
    clusters = butina_clusters(bits, cutoff)

    rng = np.random.default_rng(seed)
    order = rng.permutation(len(clusters))
    target = int(round(test_size * len(bits)))

    test_local: list[int] = []
    for c in order:
        if len(test_local) >= target:
            break
        test_local.extend(clusters[c])

    test_mask = np.zeros(len(smiles), dtype=bool)
    test_mask[index_of[np.array(test_local, dtype=int)]] = True
    train_mask = ok & ~test_mask
    return train_mask, test_mask, ok


def describe_split(
    smiles: list[str], train_mask: np.ndarray, test_mask: np.ndarray
) -> dict:
    """Nearest-neighbour similarity from test to train, as split difficulty."""
    train_bits, train_ok = fingerprints([s for s, m in zip(smiles, train_mask) if m])
    test_bits, test_ok = fingerprints([s for s, m in zip(smiles, test_mask) if m])
    if not train_bits or not test_bits:
        return {}
    similarity = nearest_neighbour_similarity(test_bits, train_bits)
    return {
        "n_train": int(train_ok.sum()),
        "n_test": int(test_ok.sum()),
        "nn_sim_mean": float(similarity.mean()),
        "nn_sim_median": float(np.median(similarity)),
        "frac_nn_above_0.5": float((similarity > 0.5).mean()),
        "frac_nn_above_0.7": float((similarity > 0.7).mean()),
    }
