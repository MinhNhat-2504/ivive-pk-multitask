"""Tobit objective for LightGBM, for labels that are detection-limit flags.

Between a quarter and a third of the AZ clearance labels sit exactly on an assay
limit (see `notes/04`). Ordinary least squares treats "3.0" as a measurement of
three, when the datum only asserts "less than three"; the model is then pulled
towards a value the experiment never observed.

The Tobit likelihood states what the experiment actually said. For an
uncensored point the contribution is Gaussian; for a left-censored point it is
the probability that the true value fell below the limit, and symmetrically
above for a right-censored one:

    uncensored      -log phi((y - p) / sigma) / sigma
    left-censored   -log Phi((c - p) / sigma)
    right-censored  -log Phi((p - c) / sigma)

Gradients and Hessians are supplied analytically. The Hessian of the censored
terms is lambda(z) * (lambda(z) + z) with lambda the inverse Mills ratio, which
lies strictly in (0, 1), so the objective stays convex in the score and
LightGBM's Newton step is well behaved.
"""

from __future__ import annotations

import numpy as np
from lightgbm import LGBMRegressor
from scipy.stats import norm

#: Floor on the Hessian. The censored terms approach zero curvature far from the
#: limit, and LightGBM divides by the summed Hessian when it computes leaf
#: values, so an unfloored term can produce enormous leaves.
MIN_HESSIAN = 1e-6


def inverse_mills(z: np.ndarray) -> np.ndarray:
    """phi(z) / Phi(z), evaluated so it stays finite in the far-left tail."""
    # norm.logpdf - norm.logcdf avoids 0/0 where both underflow; for z < -30 the
    # ratio tends to |z|, which is the asymptotic form.
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.exp(norm.logpdf(z) - norm.logcdf(z))
    return np.where(np.isfinite(ratio), ratio, np.abs(z))


def tobit_objective(
    left: np.ndarray, right: np.ndarray, sigma: float
):
    """Build a LightGBM objective closure.

    `left` and `right` are boolean masks over the training rows, aligned with
    the order LightGBM passes them back in. A row flagged in neither mask is an
    ordinary observation.
    """
    observed = ~(left | right)

    def objective(y_true: np.ndarray, y_pred: np.ndarray):
        grad = np.zeros_like(y_pred)
        hess = np.zeros_like(y_pred)

        # Ordinary Gaussian term.
        grad[observed] = (y_pred[observed] - y_true[observed]) / sigma**2
        hess[observed] = 1.0 / sigma**2

        # Left-censored: the true value lies below the limit, so only a
        # prediction that climbs above it is penalised.
        if left.any():
            z = (y_true[left] - y_pred[left]) / sigma
            lam = inverse_mills(z)
            grad[left] = lam / sigma
            hess[left] = lam * (lam + z) / sigma**2

        # Right-censored: mirror image.
        if right.any():
            z = (y_pred[right] - y_true[right]) / sigma
            lam = inverse_mills(z)
            grad[right] = -lam / sigma
            hess[right] = lam * (lam + z) / sigma**2

        return grad, np.maximum(hess, MIN_HESSIAN)

    return objective


def estimate_sigma(residuals: np.ndarray) -> float:
    """Scale of the Gaussian term, from the uncensored residuals of a first fit.

    A robust scale is used rather than the standard deviation: the uncensored
    residuals of a clearance model have heavy tails, and a few large ones would
    otherwise inflate sigma and flatten the censored gradients towards zero.
    """
    mad = float(np.median(np.abs(residuals - np.median(residuals))))
    return max(1.4826 * mad, 1e-3)


class CensoredLGBM:
    """LightGBM with a Tobit objective, fitted in two stages.

    The first stage is an ordinary fit whose uncensored residuals set sigma; the
    second refits with the Tobit objective. Estimating sigma jointly with the
    trees is possible but makes the comparison against the plain baseline less
    clean, since the two models would then differ in more than the objective.
    """

    def __init__(self, left_limit: float, right_limit: float, **params):
        self.left_limit = left_limit
        self.right_limit = right_limit
        self.params = params
        self.model_: LGBMRegressor | None = None
        self.init_score_: float = 0.0
        self.sigma_: float = 1.0

    def fit(self, x: np.ndarray, y: np.ndarray) -> "CensoredLGBM":
        left = np.isclose(y, self.left_limit)
        right = np.isclose(y, self.right_limit)
        observed = ~(left | right)

        warm = LGBMRegressor(**self.params).fit(x[observed], y[observed])
        self.sigma_ = estimate_sigma(y[observed] - warm.predict(x[observed]))

        # With a custom objective LightGBM does not start from the label mean,
        # so the offset is supplied explicitly and added back at predict time.
        self.init_score_ = float(np.mean(y[observed]))

        self.model_ = LGBMRegressor(
            objective=tobit_objective(left, right, self.sigma_), **self.params
        )
        self.model_.fit(x, y, init_score=np.full(y.size, self.init_score_))
        return self

    def predict(self, x: np.ndarray) -> np.ndarray:
        if self.model_ is None:
            raise RuntimeError("fit before predict")
        return self.model_.predict(x) + self.init_score_
