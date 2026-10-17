"""Molecular featurisation: RDKit descriptors + Morgan counts.

Gradient boosting on descriptors is the backbone of this project rather than a
graph network, because boosting currently leads the graph baselines on exactly
these PK endpoints (e.g. Caco2: CaliciBoost 0.256 vs MiniMol 0.350 MAE) and
because it trains on a CPU in seconds.

Descriptors are also what makes the IVIVE constraint interpretable later: logP,
TPSA and the acid/base counts are the same quantities that drive protein binding
and hepatic extraction in the mechanistic model, so the constrained model's
behaviour can be read against pharmacological expectation.
"""

from __future__ import annotations

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, rdFingerprintGenerator, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

MORGAN_RADIUS = 2
MORGAN_BITS = 1024

_morgan = rdFingerprintGenerator.GetMorganGenerator(
    radius=MORGAN_RADIUS, fpSize=MORGAN_BITS
)

#: Physicochemical descriptors, each a (name, callable) pair.
_DESCRIPTORS: list[tuple[str, object]] = [
    ("mw", Descriptors.MolWt),
    ("logp", Crippen.MolLogP),
    ("mr", Crippen.MolMR),
    ("tpsa", rdMolDescriptors.CalcTPSA),
    ("hbd", rdMolDescriptors.CalcNumHBD),
    ("hba", rdMolDescriptors.CalcNumHBA),
    ("rotb", rdMolDescriptors.CalcNumRotatableBonds),
    ("rings", rdMolDescriptors.CalcNumRings),
    ("aromatic_rings", rdMolDescriptors.CalcNumAromaticRings),
    ("heavy_atoms", lambda m: m.GetNumHeavyAtoms()),
    ("fsp3", rdMolDescriptors.CalcFractionCSP3),
    ("heteroatoms", rdMolDescriptors.CalcNumHeteroatoms),
    ("stereocentres", rdMolDescriptors.CalcNumAtomStereoCenters),
    ("amide_bonds", rdMolDescriptors.CalcNumAmideBonds),
    ("labute_asa", rdMolDescriptors.CalcLabuteASA),
    ("balaban_j", Descriptors.BalabanJ),
    ("bertz_ct", Descriptors.BertzCT),
    ("hall_kier_alpha", Descriptors.HallKierAlpha),
    ("kappa1", Descriptors.Kappa1),
    ("kappa2", Descriptors.Kappa2),
    ("kappa3", Descriptors.Kappa3),
    ("chi0v", Descriptors.Chi0v),
    ("chi1v", Descriptors.Chi1v),
    ("chi2v", Descriptors.Chi2v),
    ("max_partial_charge", Descriptors.MaxPartialCharge),
    ("min_partial_charge", Descriptors.MinPartialCharge),
    ("num_valence_electrons", Descriptors.NumValenceElectrons),
    ("qed", Descriptors.qed),
]

DESCRIPTOR_NAMES: list[str] = [name for name, _ in _DESCRIPTORS]
FEATURE_NAMES: list[str] = DESCRIPTOR_NAMES + [f"morgan_{i}" for i in range(MORGAN_BITS)]


def featurise_one(smiles: str) -> np.ndarray | None:
    """Return the feature vector for one SMILES, or None if it does not parse."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    values = np.empty(len(_DESCRIPTORS), dtype=np.float64)
    for i, (_, fn) in enumerate(_DESCRIPTORS):
        try:
            values[i] = fn(mol)
        except Exception:
            # A few descriptors (the partial-charge pair especially) fail on
            # unusual valences. LightGBM handles NaN natively, so leaving the
            # gap is better than substituting a value that is not measured.
            values[i] = np.nan

    counts = _morgan.GetCountFingerprintAsNumPy(mol).astype(np.float64)
    return np.concatenate([values, counts])


def featurise(smiles: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """Featurise many SMILES.

    Returns the feature matrix and a boolean mask of which inputs parsed, so the
    caller can drop the matching labels.
    """
    rows: list[np.ndarray] = []
    ok = np.zeros(len(smiles), dtype=bool)
    for i, smi in enumerate(smiles):
        vec = featurise_one(smi)
        if vec is not None:
            rows.append(vec)
            ok[i] = True
    matrix = np.vstack(rows) if rows else np.empty((0, len(FEATURE_NAMES)))
    # Infinities arise from a handful of topological indices on degenerate
    # graphs; map them to NaN so they take the same missing-value path.
    return np.where(np.isfinite(matrix), matrix, np.nan), ok
