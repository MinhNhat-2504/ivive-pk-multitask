"""The IVIVE relation that couples the five endpoints.

Chain, for a hepatically eliminated compound:

    CL_int,in vivo = CL_int,in vitro * scaling factor      (physiological scaling)
    CL_h           = Q_h * fu_b * CL_int / (Q_h + fu_b * CL_int)   (well-stirred)
    t_1/2          = ln2 * V_ss / CL_h                     (one-compartment)

Plasma protein binding enters as fu, an intermediate of the chain, rather than
as an independent prediction target -- which is the respect in which this
differs from PEMAL, where the coupling is a pair of constants fitted on mouse
data.

Every constant below is a population mean. Real compounds deviate, and the
deviation is not noise to be averaged away: extrahepatic and renal elimination,
transporter-mediated uptake, and blood-to-plasma partitioning all break the
relation for identifiable classes of compound. `validate_relation` measures how
badly, on the molecules where all the labels needed to check it exist.
"""

from __future__ import annotations

import numpy as np

# -- Human physiological scaling constants -------------------------------------
# Barter et al. (2007) consensus values; see notes for the citation to add.
HEPATOCELLULARITY = 99.0  # 10^6 cells per gram of liver
MICROSOMAL_PROTEIN = 32.0  # mg microsomal protein per gram of liver
LIVER_WEIGHT = 25.7  # gram of liver per kg body weight
HEPATIC_BLOOD_FLOW = 20.7  # mL/min/kg body weight

#: uL/min/10^6 cells -> mL/min/kg body weight
HEPATOCYTE_SCALE = HEPATOCELLULARITY * LIVER_WEIGHT / 1000.0
#: uL/min/mg protein -> mL/min/kg body weight
MICROSOME_SCALE = MICROSOMAL_PROTEIN * LIVER_WEIGHT / 1000.0

#: The smallest unbound fraction the binding assay can express. Its declared
#: range stops at 99.95% bound, so a prediction beyond that is clipped here
#: rather than at some arbitrary value: anything further is outside what the
#: measurement could have reported. Using a tighter clip inflates fold errors
#: for models that saturate -- clipping at 99.99% instead turns a ridge model's
#: median fold error from 20 into 30 -- and the inflation is a property of the
#: clip, not of the model.
ASSAY_FU_FLOOR = 5e-4

LN2 = float(np.log(2.0))


def unbound_fraction(ppbr_percent: np.ndarray, floor: float = 1e-3) -> np.ndarray:
    """fu from percent bound.

    Floored because highly bound compounds are reported as 99.x% and the
    well-stirred model divides by fu downstream; an fu of exactly zero would
    make hepatic clearance zero and the predicted half-life infinite.
    """
    return np.clip(1.0 - np.asarray(ppbr_percent, dtype=float) / 100.0, floor, 1.0)


def scale_intrinsic_clearance(clint_in_vitro: np.ndarray, system: str) -> np.ndarray:
    """In vitro intrinsic clearance to whole-body units (mL/min/kg)."""
    factor = {"hepatocyte": HEPATOCYTE_SCALE, "microsome": MICROSOME_SCALE}[system]
    return np.asarray(clint_in_vitro, dtype=float) * factor


def well_stirred_clearance(clint_scaled: np.ndarray, fu: np.ndarray) -> np.ndarray:
    """Hepatic clearance (mL/min/kg) from scaled intrinsic clearance and fu.

    Blood-to-plasma ratio is taken as 1, so fu in blood equals fu in plasma.
    That assumption is wrong for compounds that partition into erythrocytes and
    is one of the known sources of the deviation measured in validate_relation.
    """
    numerator = HEPATIC_BLOOD_FLOW * fu * clint_scaled
    return numerator / (HEPATIC_BLOOD_FLOW + fu * clint_scaled)


def half_life_from(vdss_l_per_kg: np.ndarray, clearance_ml_min_kg: np.ndarray) -> np.ndarray:
    """Terminal half-life in hours from V_ss (L/kg) and clearance (mL/min/kg)."""
    clearance_l_per_h_kg = np.asarray(clearance_ml_min_kg, dtype=float) * 60.0 / 1000.0
    return LN2 * np.asarray(vdss_l_per_kg, dtype=float) / clearance_l_per_h_kg


def predict_half_life(
    clint_in_vitro: np.ndarray,
    ppbr_percent: np.ndarray,
    vdss_l_per_kg: np.ndarray,
    system: str = "hepatocyte",
) -> np.ndarray:
    """Run the whole chain: in vitro clearance + binding + volume -> half-life."""
    fu = unbound_fraction(ppbr_percent)
    clint = scale_intrinsic_clearance(clint_in_vitro, system)
    return half_life_from(vdss_l_per_kg, well_stirred_clearance(clint, fu))


def consistency_residual(
    clint_in_vitro: np.ndarray,
    ppbr_percent: np.ndarray,
    vdss_l_per_kg: np.ndarray,
    half_life_h: np.ndarray,
    system: str = "hepatocyte",
) -> np.ndarray:
    """log10 ratio of observed half-life to the half-life the chain implies.

    Expressed in log10 so that a value of 1.0 means "wrong by a factor of ten",
    which is the unit the IVIVE literature reports accuracy in.
    """
    predicted = predict_half_life(clint_in_vitro, ppbr_percent, vdss_l_per_kg, system)
    return np.log10(np.asarray(half_life_h, dtype=float) / predicted)
