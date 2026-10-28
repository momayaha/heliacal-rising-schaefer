#!/usr/bin/env python3
"""Schaefer-style naked-eye stellar visibility calculations.

This module implements the *visibility-physics* layer intended for the
Chiang Saen heliacal-rising revision.  It deliberately does NOT calculate
solar/stellar ephemerides or search across dates; those tasks belong in
separate geometry/search modules.

Primary sources implemented here
--------------------------------
Schaefer, B. E. (1986), "Atmospheric Extinction Effects on Stellar
Alignments", Archaeoastronomy 10 / JHA 17, S32-S42.

Schaefer, B. E. (1987), "Heliacal Rise Phenomena", Archaeoastronomy 11 /
JHA 18, S19-S33.

Schaefer, B. E. (1993), "Astronomy and the Limits of Vision", Vistas in
Astronomy 36, 311-361.

Model scope
-----------
For heliacal work, Schaefer (1987) models the twilight background at the
star's position and combines it with atmospheric extinction and a human-eye
point-source threshold.  This file follows that architecture:

    geometry -> twilight background -> air mass/extinction
             -> human-eye threshold -> visible / not visible

All angles accepted by the public API are in DEGREES.  Schaefer (1987)
Eq. (1), however, uses angular quantities in RADIANS internally.

The project baseline uses m_z = 6.0 and k_v in {0.20, 0.30} mag/airmass.
The functions accept other k_v values so that published Schaefer cases can
be reproduced during validation.

Important limitation
--------------------
The empirical twilight relation of Schaefer (1987), Eq. (1), was optimized
for approximately:

    -18 deg < Sun altitude < -6 deg
    stellar zenith distance > 60 deg  (star altitude < 30 deg)
    Sun-star azimuth difference < 90 deg

The code therefore reports whether a calculation lies inside that
calibration domain.  It does not silently claim validity outside it.

No third-party packages are required.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Tuple


# ---------------------------------------------------------------------------
# Project defaults / published constants
# ---------------------------------------------------------------------------

PROJECT_KV_SCENARIOS: Tuple[float, float] = (0.20, 0.30)
DEFAULT_MZ = 6.0

# Schaefer (1987), Eq. (2): the standard-sky reference background.
# Units in the paper are millimicroLamberts (mμL), numerically equivalent
# to nanoLamberts (nL).
STANDARD_DARK_SKY_NL = 118.0
STANDARD_KV = 0.20

# Hecht point-source threshold constants as reproduced by
# Schaefer (1986), Eq. (4), and Schaefer (1993), Eq. (34).
# Branch boundary: log10(background[nL]) = 3.17.
HECHT_LOG_B_BREAK = 3.17
HECHT_LOW_K = 10.0 ** (-1.90)
HECHT_LOW_C = 10.0 ** (-9.80)
HECHT_HIGH_K = 10.0 ** (-5.90)
HECHT_HIGH_C = 10.0 ** (-8.35)

# Visual-magnitude / illuminance zero point used by Schaefer
# (1986 Eq. 5; 1993 Eq. 10).  Illuminance is in foot-candles.
MAG_ZEROPOINT = -16.57


@dataclass(frozen=True)
class VisibilityResult:
    """Diagnostic output for one star at one instant."""

    visible: bool
    star_magnitude_v: float
    limiting_magnitude_v: float
    visibility_margin_mag: float
    # Positive margin means the star is brighter than the threshold.

    background_nl: float
    standard_twilight_background_nl: float
    threshold_illuminance_fc: float
    star_illuminance_above_atmosphere_fc: float
    star_illuminance_at_eye_fc: float

    airmass: float
    extinction_mag: float

    sun_altitude_deg: float
    star_altitude_deg: float
    star_zenith_distance_deg: float
    sun_star_azimuth_difference_deg: float
    kv: float
    mz: float

    within_schaefer1987_calibration_domain: bool


# ---------------------------------------------------------------------------
# Human-eye threshold (Hecht 1947 as used by Schaefer)
# ---------------------------------------------------------------------------


def hecht_constants(background_nl: float) -> Tuple[float, float]:
    """Return (K, C) for the Hecht point-source threshold relation.

    Source: Schaefer (1986), Eq. (4); Schaefer (1993), Eq. (34).

    Parameters
    ----------
    background_nl : float
        Effective background luminance in nanoLamberts (nL), equivalent
        numerically to the mμL unit used in Schaefer (1986/1987).
    """
    if background_nl <= 0.0:
        raise ValueError("background_nl must be positive")

    if math.log10(background_nl) < HECHT_LOG_B_BREAK:
        return HECHT_LOW_K, HECHT_LOW_C
    return HECHT_HIGH_K, HECHT_HIGH_C


def threshold_illuminance(background_nl: float) -> float:
    """Minimum detectable point-source illuminance at the eye [foot-candle].

    Implements Hecht's relation as reproduced by Schaefer (1986), Eq. (4),
    and Schaefer (1993), Eq. (34):

        i_th = C * (1 + sqrt(K * b))^2

    where b is the effective background in nL (mμL in the original papers).
    """
    K, C = hecht_constants(background_nl)
    return C * (1.0 + math.sqrt(K * background_nl)) ** 2


# ---------------------------------------------------------------------------
# Air mass and stellar extinction
# ---------------------------------------------------------------------------


def rozenberg_airmass_from_zenith_distance(zenith_distance_deg: float) -> float:
    """Rozenberg (1966) air mass used by Schaefer near the horizon.

    Schaefer (1986), Eq. (6), and Schaefer (1993), Eq. (2b):

        X = 1 / [cos(Z) + 0.025 exp(-11 cos(Z))]

    Z is the apparent zenith distance.

    Notes
    -----
    The formula gives X ~ 40 at the geometric horizon (Z = 90 deg),
    avoiding the divergence of the simple sec(Z) approximation.
    """
    if not (0.0 <= zenith_distance_deg <= 90.0):
        raise ValueError(
            "zenith_distance_deg must be between 0 and 90 degrees "
            "for this above-horizon visibility calculation"
        )

    z = math.radians(zenith_distance_deg)
    cos_z = math.cos(z)
    denominator = cos_z + 0.025 * math.exp(-11.0 * cos_z)
    return 1.0 / denominator


def rozenberg_airmass_from_altitude(star_altitude_deg: float) -> float:
    """Convenience wrapper returning air mass from apparent altitude."""
    if not (0.0 <= star_altitude_deg <= 90.0):
        raise ValueError("star_altitude_deg must be between 0 and 90 degrees")
    return rozenberg_airmass_from_zenith_distance(90.0 - star_altitude_deg)


def extinction_magnitude(kv: float, airmass: float) -> float:
    """Simple visual extinction Δm = k_v X [magnitudes].

    This is the scalar-k form used in Schaefer (1986), Eq. (5), and in the
    Schaefer (1987) heliacal model.  Schaefer (1993) also discusses a more
    detailed component-by-component near-horizon treatment; that refinement
    is intentionally kept separate from this 1987-compatible core.
    """
    if kv < 0.0:
        raise ValueError("kv must be non-negative")
    if airmass <= 0.0:
        raise ValueError("airmass must be positive")
    return kv * airmass


# ---------------------------------------------------------------------------
# Schaefer (1987) heliacal twilight background
# ---------------------------------------------------------------------------


def standard_twilight_background(
    sun_altitude_deg: float,
    star_zenith_distance_deg: float,
    sun_star_azimuth_difference_deg: float,
) -> float:
    """Standard-sky twilight luminance at the star position [nL].

    Implements Schaefer (1987), Eq. (1):

        log10(b_ST) = 4.75 - theta*Z/3
                      + h0*(12 + 8.21*Z) + 2.86*Z

    h0 is the Sun's altitude (negative below the horizon), Z is the stellar
    zenith distance, and theta is the Sun-star azimuth difference.  Schaefer
    states that these angular quantities are measured in RADIANS.

    The returned unit is nL, numerically equivalent to the mμL used in the
    original paper.
    """
    if not (-90.0 <= sun_altitude_deg <= 0.0):
        raise ValueError("sun_altitude_deg must be between -90 and 0 degrees")
    if not (0.0 <= star_zenith_distance_deg <= 90.0):
        raise ValueError("star_zenith_distance_deg must be between 0 and 90 degrees")
    if not (0.0 <= sun_star_azimuth_difference_deg <= 180.0):
        raise ValueError(
            "sun_star_azimuth_difference_deg must be between 0 and 180 degrees"
        )

    h0 = math.radians(sun_altitude_deg)
    z = math.radians(star_zenith_distance_deg)
    theta = math.radians(sun_star_azimuth_difference_deg)

    log10_b_st = 4.75 - (theta * z / 3.0) + h0 * (12.0 + 8.21 * z) + 2.86 * z
    return 10.0 ** log10_b_st


def _dark_sky_background_from_mz(mz: float, kv: float) -> float:
    """Dark-sky background implied by (m_z, k_v) in Schaefer (1987) Eq. 2.

    Schaefer (1987) defines

        Y = -0.2 [m_z + 16.57 + k_v + 2.5 log10(C)]
        b_dark = (10^Y - 1)^2 / K

    where C and K are the Hecht constants.  Because the Hecht constants have
    two branches, this helper checks which branch is self-consistent with the
    resulting background.
    """
    if kv < 0.0:
        raise ValueError("kv must be non-negative")

    candidates = []
    for K, C, branch_name in (
        (HECHT_LOW_K, HECHT_LOW_C, "low"),
        (HECHT_HIGH_K, HECHT_HIGH_C, "high"),
    ):
        Y = -0.2 * (mz + 16.57 + kv + 2.5 * math.log10(C))
        b_dark = (10.0 ** Y - 1.0) ** 2 / K
        if b_dark <= 0.0:
            continue
        is_low = math.log10(b_dark) < HECHT_LOG_B_BREAK
        if (branch_name == "low" and is_low) or (branch_name == "high" and not is_low):
            candidates.append(b_dark)

    if len(candidates) == 1:
        return candidates[0]
    if len(candidates) > 1:
        # This should be rare; choose the smaller-background solution, which
        # corresponds to the dark-adapted branch relevant to naked-eye night sky.
        return min(candidates)

    raise RuntimeError(
        "No self-consistent Hecht branch found for the supplied mz and kv. "
        "Check whether the parameters are physically appropriate."
    )


def adjusted_twilight_background(
    standard_background_nl: float,
    kv: float,
    mz: float = DEFAULT_MZ,
) -> float:
    """Adjust Schaefer's standard twilight sky to non-standard conditions.

    Implements Schaefer (1987), Eq. (2):

        b = (10^Y - 1)^2 / K
            + (k_v / 0.20) * (b_ST - 118)

    with

        Y = -0.2 [m_z + 16.57 + k_v + 2.5 log10(C)].

    The first term is the dark-sky background implied by m_z and k_v.

    Implementation note
    -------------------
    The twilight fit b_ST is not allowed to drive the total luminance below
    the dark-sky level represented by m_z.  For the standard sky this means
    b_ST is floored at 118 nL.  Without this physical floor, faint-star cases
    can acquire an artificial visibility window against a background darker
    than the model's own specified dark sky.
    """
    if standard_background_nl <= 0.0:
        raise ValueError("standard_background_nl must be positive")

    b_dark = _dark_sky_background_from_mz(mz, kv)

    # Dark-sky floor.  Eq. (1) is an empirical *twilight* fit and, near the
    # dark end of its useful range, its extrapolated b_ST can fall below the
    # standard dark-sky background (118 nL) that Schaefer explicitly adopts.
    # A physical sky background cannot become darker than the dark-sky level
    # specified by m_z merely because the twilight term has become small.
    # We therefore floor the standard-sky term at 118 nL before applying the
    # non-standard-sky correction.  This also reproduces the faint-star
    # transition in Schaefer (1987), Table 4 (m=3--5), where first visibility
    # occurs near the onset of twilight rather than in an artificially
    # sub-dark background.
    b_st_effective = max(standard_background_nl, STANDARD_DARK_SKY_NL)
    b = b_dark + (kv / STANDARD_KV) * (
        b_st_effective - STANDARD_DARK_SKY_NL
    )

    if b <= 0.0:
        raise ValueError(
            "Adjusted twilight background is non-positive. "
            "Check the supplied m_z and k_v values."
        )
    return b


# ---------------------------------------------------------------------------
# Magnitudes, illuminance, and visibility decision
# ---------------------------------------------------------------------------


def illuminance_above_atmosphere(star_magnitude_v: float) -> float:
    """Visual illuminance above the atmosphere [foot-candle].

    Rearrangement of Schaefer (1993), Eq. (10):

        m_v = -16.57 - 2.5 log10(I_0)
    """
    return 10.0 ** ((MAG_ZEROPOINT - star_magnitude_v) / 2.5)


def illuminance_at_eye(
    star_magnitude_v: float,
    kv: float,
    airmass: float,
) -> float:
    """Apparent stellar illuminance after atmospheric extinction [fc]."""
    i0 = illuminance_above_atmosphere(star_magnitude_v)
    dm = extinction_magnitude(kv, airmass)
    return i0 * 10.0 ** (-0.4 * dm)


def limiting_magnitude(
    background_nl: float,
    kv: float,
    airmass: float,
) -> float:
    """Faintest above-atmosphere V magnitude detectable at this geometry.

    Schaefer (1986), Eq. (5):

        m_lim = -16.57 - 2.5 log10(i_th) - k_v X

    where i_th is the Hecht threshold illuminance at the eye.
    """
    i_th = threshold_illuminance(background_nl)
    return MAG_ZEROPOINT - 2.5 * math.log10(i_th) - kv * airmass


def in_schaefer1987_calibration_domain(
    sun_altitude_deg: float,
    star_zenith_distance_deg: float,
    sun_star_azimuth_difference_deg: float,
) -> bool:
    """Whether geometry lies in the optimized range quoted by Schaefer 1987."""
    return (
        -18.0 < sun_altitude_deg < -6.0
        and star_zenith_distance_deg > 60.0
        and sun_star_azimuth_difference_deg < 90.0
    )


def evaluate_visibility(
    *,
    star_magnitude_v: float,
    sun_altitude_deg: float,
    star_altitude_deg: float,
    sun_star_azimuth_difference_deg: float,
    kv: float,
    mz: float = DEFAULT_MZ,
) -> VisibilityResult:
    """Evaluate naked-eye visibility for one star at one instant.

    Parameters
    ----------
    star_magnitude_v : float
        Above-atmosphere visual magnitude of the star.
    sun_altitude_deg : float
        Apparent solar altitude; negative below the horizon.
    star_altitude_deg : float
        Apparent stellar altitude above an ideal horizon.
    sun_star_azimuth_difference_deg : float
        Absolute azimuth difference between Sun and star, 0..180 degrees.
    kv : float
        Visual extinction coefficient [mag/airmass].  Chiang Saen project
        scenarios are 0.20 and 0.30.
    mz : float, default 6.0
        Dark-sky zenith limiting magnitude.

    Returns
    -------
    VisibilityResult
        Full diagnostics.  The star is classified as visible when
        star_magnitude_v <= limiting_magnitude_v, equivalently when the
        visibility margin is >= 0 mag.
    """
    if not (0.0 <= star_altitude_deg <= 90.0):
        raise ValueError("star_altitude_deg must be between 0 and 90 degrees")

    z_deg = 90.0 - star_altitude_deg
    x = rozenberg_airmass_from_zenith_distance(z_deg)

    b_st = standard_twilight_background(
        sun_altitude_deg=sun_altitude_deg,
        star_zenith_distance_deg=z_deg,
        sun_star_azimuth_difference_deg=sun_star_azimuth_difference_deg,
    )
    b = adjusted_twilight_background(b_st, kv=kv, mz=mz)

    i_th = threshold_illuminance(b)
    i0 = illuminance_above_atmosphere(star_magnitude_v)
    i_eye = illuminance_at_eye(star_magnitude_v, kv=kv, airmass=x)

    m_lim = limiting_magnitude(b, kv=kv, airmass=x)
    margin = m_lim - star_magnitude_v

    return VisibilityResult(
        visible=(margin >= 0.0),
        star_magnitude_v=star_magnitude_v,
        limiting_magnitude_v=m_lim,
        visibility_margin_mag=margin,
        background_nl=b,
        standard_twilight_background_nl=b_st,
        threshold_illuminance_fc=i_th,
        star_illuminance_above_atmosphere_fc=i0,
        star_illuminance_at_eye_fc=i_eye,
        airmass=x,
        extinction_mag=kv * x,
        sun_altitude_deg=sun_altitude_deg,
        star_altitude_deg=star_altitude_deg,
        star_zenith_distance_deg=z_deg,
        sun_star_azimuth_difference_deg=sun_star_azimuth_difference_deg,
        kv=kv,
        mz=mz,
        within_schaefer1987_calibration_domain=in_schaefer1987_calibration_domain(
            sun_altitude_deg,
            z_deg,
            sun_star_azimuth_difference_deg,
        ),
    )


# ---------------------------------------------------------------------------
# Lightweight internal checks (not a substitute for published-case validation)
# ---------------------------------------------------------------------------


def _self_test() -> None:
    """Run basic algebra/unit sanity checks."""

    # Rozenberg air mass should be ~40 at the horizon and ~1 at zenith.
    assert abs(rozenberg_airmass_from_altitude(0.0) - 40.0) < 1e-9
    assert abs(rozenberg_airmass_from_altitude(90.0) - 1.0) < 1e-3

    # Schaefer's standard m_z=6, k_v=0.20 should imply ~118 nL dark sky.
    b_dark = _dark_sky_background_from_mz(6.0, 0.20)
    assert abs(b_dark - 118.0) < 0.2

    # Eq. (2) should reproduce b_ST to rounding under standard conditions.
    for b_st in (200.0, 1000.0, 10000.0):
        b = adjusted_twilight_background(b_st, kv=0.20, mz=6.0)
        assert abs(b - b_st) < 0.2

    # Magnitude/illuminance inversion check.
    for m in (-1.46, 0.0, 2.86, 5.0):
        i0 = illuminance_above_atmosphere(m)
        recovered = MAG_ZEROPOINT - 2.5 * math.log10(i0)
        assert abs(recovered - m) < 1e-12


if __name__ == "__main__":
    _self_test()
    print("schaefer_visibility.py: internal checks passed")
    print(f"Project k_v scenarios: {PROJECT_KV_SCENARIOS}")
