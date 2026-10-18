"""Multitask gradient boosting with the IVIVE constraint applied to predictions.

`notes/01` showed the constraint cannot be supervised: only 20 molecules carry
the labels needed to close the chain. `notes/06` then showed which coupling is
actually worth enforcing -- the hepatocyte/microsome pair, where 395 uncensored
molecules carry both labels, rather than the full chain where 16 do.

So the constraint acts on predictions. Every molecule gets a prediction from
every head whether or not it carries that label, the two clearance heads are
required to agree up to a learned affine correction between the assay systems,
and the requirement therefore applies to all 3 898 molecules in the union --
including the 2 791 that carry a single label and would otherwise contribute
nothing to any task but their own.

Learning that correction rather than fixing it is the point of departure from
PEMAL, which uses constants fitted on mouse data. The two constants this data
calls for were measured in `notes/06`: a 3.8-fold offset between the systems,
and a 2.4-fold under-prediction from microsomes.

Algorithm, one round:

    1. fit each head on its own labels
    2. predict every head for every molecule
    3. project each molecule's pair of clearance predictions onto the
       constraint line, in closed form, weighted by how much each head is
       trusted -- a head holding a real label is held fixed
    4. refit on labels where they exist and projections where they do not,
       with projections down-weighted
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
from lightgbm import LGBMRegressor

from src.models.ivive import HEPATOCYTE_SCALE, MICROSOME_SCALE

log = logging.getLogger(__name__)

#: Default weight standing for "this value is a measurement". Not zero, because
#: a label is itself noisy. The original justification quoted a within-molecule
#: spread, which `notes/15` retracted; the value is kept as a default and swept
#: in `src/experiments/weight_sweep.py` rather than defended from a number that
#: no longer exists.
LABEL_WEIGHT = 0.05
#: Weight for a value the model merely predicted, so the projection moves it far
#: more readily than it moves a measurement.
PREDICTION_WEIGHT = 1.0
#: Below this many usable paired molecules the coupling is not fitted and the
#: model falls back to pure physiological scaling. `notes/16` found this rule,
#: not the choice of coupling, is what the 25%-label ablation row measures.
MIN_COUPLING_PAIRS = 20


@dataclass
class AffineCoupling:
    """log10 CL_mic_scaled = intercept + slope * log10 CL_hep_scaled."""

    intercept: float = 0.0
    slope: float = 1.0
    n_fitted: int = 0
    residual_sd: float = field(default=np.nan)

    def apply(self, log_hep: np.ndarray) -> np.ndarray:
        return self.intercept + self.slope * log_hep


def deming_fit(x: np.ndarray, y: np.ndarray, variance_ratio: float = 1.0) -> tuple[float, float]:
    """Errors-in-variables line fit.

    Ordinary least squares is the wrong estimator here. It assumes the
    predictor is measured without error, and attenuates the slope towards zero
    by roughly the predictor's reliability when that assumption fails. These
    are two noisy assays of the same quantity -- `notes/04` puts the hepatocyte
    ICC near 0.57 -- so OLS returns a slope near 0.4 where the relation is much
    closer to proportional, and the projection built on it would drag
    predictions towards the mean.

    `variance_ratio` is sigma_y^2 / sigma_x^2, left at 1 because both are the
    same class of assay run in the same laboratory; at 1 this is orthogonal
    regression.
    """
    x_centred, y_centred = x - x.mean(), y - y.mean()
    sxx = float(x_centred @ x_centred)
    syy = float(y_centred @ y_centred)
    sxy = float(x_centred @ y_centred)
    if abs(sxy) < 1e-12:
        return 0.0, float(y.mean())
    discriminant = (syy - variance_ratio * sxx) ** 2 + 4.0 * variance_ratio * sxy**2
    slope = (syy - variance_ratio * sxx + np.sqrt(discriminant)) / (2.0 * sxy)
    return float(slope), float(y.mean() - slope * x.mean())


def fit_coupling(
    clint_hep: np.ndarray,
    clint_mic: np.ndarray,
    censored: np.ndarray,
    mode: str = "deming",
    min_pairs: int = MIN_COUPLING_PAIRS,
) -> AffineCoupling:
    """Fit the between-system correction on molecules measured in both systems.

    Censored molecules are excluded. `notes/06` measured what including them
    does: the rank agreement between the two assays reads 0.642 with them and
    0.412 without, because most of the apparent agreement is the two assays
    both reporting the same limit flag rather than both measuring the compound.
    """
    usable = ~censored & np.isfinite(clint_hep) & np.isfinite(clint_mic)
    log_hep = np.log10(np.clip(clint_hep[usable] * HEPATOCYTE_SCALE, 1e-6, None))
    log_mic = np.log10(np.clip(clint_mic[usable] * MICROSOME_SCALE, 1e-6, None))

    if mode == "identity":
        # Ablation standing in for the textbook assumption, and for PEMAL's use
        # of coefficients fixed in advance: physiological scaling alone is taken
        # to make the two systems agree, with nothing learned from the data.
        residuals = log_mic - log_hep
        return AffineCoupling(
            intercept=0.0,
            slope=1.0,
            n_fitted=int(usable.sum()),
            residual_sd=float(residuals.std(ddof=1)) if usable.sum() > 1 else np.nan,
        )

    if usable.sum() < min_pairs:
        log.warning("only %d usable pairs (< %d); falling back to identity",
                    int(usable.sum()), min_pairs)
        return AffineCoupling(n_fitted=int(usable.sum()))

    if mode == "ols":
        slope, intercept = np.polyfit(log_hep, log_mic, deg=1)
    elif mode == "deming":
        slope, intercept = deming_fit(log_hep, log_mic)
    else:
        raise ValueError(f"unknown coupling mode: {mode}")
    residuals = log_mic - (intercept + slope * log_hep)
    return AffineCoupling(
        intercept=float(intercept),
        slope=float(slope),
        n_fitted=int(usable.sum()),
        residual_sd=float(residuals.std(ddof=1)),
    )


def project_pair(
    log_hep: np.ndarray,
    log_mic: np.ndarray,
    coupling: AffineCoupling,
    weight_hep: np.ndarray,
    weight_mic: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Nearest point on the constraint line, under the given weights.

    Minimises (u - log_hep)^2 / w_h + (v - log_mic)^2 / w_m subject to
    v = intercept + slope * u, which has the closed form below. A head carrying
    a measurement gets a small weight and so barely moves; a purely predicted
    head absorbs most of the correction.
    """
    a, b = coupling.intercept, coupling.slope
    denominator = 1.0 / weight_hep + b**2 / weight_mic
    numerator = log_hep / weight_hep + b * (log_mic - a) / weight_mic
    projected_hep = numerator / denominator
    return projected_hep, coupling.apply(projected_hep)


class ConstrainedMultitask:
    """Two clearance heads coupled by a learned IVIVE correction.

    The remaining three heads (ppbr, vdss, half_life) are fitted independently
    for now and enter the chain only at validation time; `notes/06` found the
    full chain too weakly supported (n=16) to drive training.
    """

    def __init__(
        self,
        n_rounds: int = 2,
        pseudo_weight: float = 0.3,
        coupling_mode: str = "deming",
        label_weight: float = LABEL_WEIGHT,
        prediction_weight: float = PREDICTION_WEIGHT,
        min_coupling_pairs: int = MIN_COUPLING_PAIRS,
        random_state: int = 0,
        **lgbm_params,
    ):
        self.n_rounds = n_rounds
        self.pseudo_weight = pseudo_weight
        self.coupling_mode = coupling_mode
        self.label_weight = label_weight
        self.prediction_weight = prediction_weight
        self.min_coupling_pairs = min_coupling_pairs
        self.random_state = random_state
        self.lgbm_params = lgbm_params
        self.heads_: dict[str, LGBMRegressor] = {}
        self.coupling_: AffineCoupling | None = None

    def fit(
        self,
        x: np.ndarray,
        y_hep: np.ndarray,
        y_mic: np.ndarray,
        censored_hep: np.ndarray,
        censored_mic: np.ndarray,
    ) -> "ConstrainedMultitask":
        """Labels are arrays over all rows with NaN where the label is absent."""
        has_hep, has_mic = np.isfinite(y_hep), np.isfinite(y_mic)
        both = has_hep & has_mic
        self.coupling_ = fit_coupling(
            y_hep[both],
            y_mic[both],
            (censored_hep | censored_mic)[both],
            mode=self.coupling_mode,
            min_pairs=self.min_coupling_pairs,
        )
        log.info(
            "coupling: log10 CLmic = %.3f + %.3f * log10 CLhep  (n=%d, residual sd %.3f)",
            self.coupling_.intercept,
            self.coupling_.slope,
            self.coupling_.n_fitted,
            self.coupling_.residual_sd,
        )

        # Work in log10 of the physiologically scaled quantity, the space the
        # constraint is linear in.
        log_hep = np.where(has_hep, np.log10(np.clip(y_hep * HEPATOCYTE_SCALE, 1e-6, None)), np.nan)
        log_mic = np.where(has_mic, np.log10(np.clip(y_mic * MICROSOME_SCALE, 1e-6, None)), np.nan)

        targets = {"hep": log_hep.copy(), "mic": log_mic.copy()}
        weights = {"hep": np.where(has_hep, 1.0, 0.0), "mic": np.where(has_mic, 1.0, 0.0)}

        for round_index in range(self.n_rounds):
            for name in ("hep", "mic"):
                labelled = weights[name] > 0
                self.heads_[name] = LGBMRegressor(
                    random_state=self.random_state, **self.lgbm_params
                ).fit(x[labelled], targets[name][labelled], sample_weight=weights[name][labelled])

            if round_index == self.n_rounds - 1:
                break

            pred_hep = self.heads_["hep"].predict(x)
            pred_mic = self.heads_["mic"].predict(x)
            # A measured head is trusted; a predicted one is not.
            w_hep = np.where(has_hep, self.label_weight, self.prediction_weight)
            w_mic = np.where(has_mic, self.label_weight, self.prediction_weight)
            proj_hep, proj_mic = project_pair(
                np.where(has_hep, log_hep, pred_hep),
                np.where(has_mic, log_mic, pred_mic),
                self.coupling_,
                w_hep,
                w_mic,
            )

            # Labels stay; everything else becomes a down-weighted pseudo-target.
            targets["hep"] = np.where(has_hep, log_hep, proj_hep)
            targets["mic"] = np.where(has_mic, log_mic, proj_mic)
            weights["hep"] = np.where(has_hep, 1.0, self.pseudo_weight)
            weights["mic"] = np.where(has_mic, 1.0, self.pseudo_weight)

        return self

    def predict(self, x: np.ndarray, head: str) -> np.ndarray:
        """Prediction on the measured in vitro scale (uL/min per cells or mg)."""
        scale = {"hep": HEPATOCYTE_SCALE, "mic": MICROSOME_SCALE}[head]
        return 10.0 ** self.heads_[head].predict(x) / scale
