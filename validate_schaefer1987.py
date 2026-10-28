#!/usr/bin/env python3
"""Numerical cross-check against Schaefer (1987), Table 3.

Purpose
-------
This script validates the *visibility* implementation in
``schaefer_visibility.py`` against the published heliacal-rising results for
Sirius in Bradley E. Schaefer (1987), "Heliacal Rise Phenomena",
Archaeoastronomy 11 / Journal for the History of Astronomy 18, S19-S33.

The paper's Table 3 gives, for Sirius at epoch A.D. 1950 and latitude +30 deg,
with m_z = 6.0:

    k_v = 0.20 -> Delta = 133.8 days
    k_v = 0.30 -> Delta = 136.3 days
    k_v = 0.40 -> Delta = 138.5 days

where Delta is the number of days elapsed since the vernal equinox.

The Sirius parameters printed by Schaefer (1987, p. S25) are:

    m_V = -1.46
    alpha = 101.09 deg
    delta = -16.69 deg
    epoch = A.D. 1950

Method
------
For each trial Delta, the script:

1. Computes an approximate solar right ascension and declination from a
   low-order solar-orbit model anchored so that solar ecliptic longitude is
   exactly 0 deg at Delta = 0 (the vernal equinox).
2. Scans the pre-sunrise sky in solar hour angle.
3. At every trial instant, derives Sun/star altitude and azimuth and calls
   ``schaefer_visibility.evaluate_visibility``.
4. Records the maximum visibility margin achieved during that dawn.
5. Solves for the Delta at which the maximum margin first reaches zero.

This is an *independent numerical cross-check*, not a byte-for-byte recreation
of Schaefer's original program.  Schaefer states that his solar position was
obtained from equations in the 1987 Nautical Almanac, but those equations are
not reproduced in the paper.  Accordingly, agreement within about half a day
is an appropriate target for this stand-alone validator.  The scientific
analysis code may later use a higher-precision ephemeris layer while retaining
the same visibility engine.

A note on Schaefer (1987), Eq. (2)
---------------------------------
The scanned 1987 paper prints the first term of Eq. (2) with ``(10^Y + 1)^2``.
However, solving the Hecht threshold relation algebraically requires
``(10^Y - 1)^2``; Schaefer (1986), Eq. (3), explicitly prints the minus sign,
and only the minus sign reproduces Schaefer's stated standard dark-sky value
of 118 m-microLambert for m_z=6.0, k_v=0.20.  The visibility module therefore
uses the algebraically consistent minus-sign form.

No third-party packages are required.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import math
from pathlib import Path
from typing import Iterable, Optional

import schaefer_visibility as sv


# ---------------------------------------------------------------------------
# Published Schaefer (1987) Sirius validation case
# ---------------------------------------------------------------------------

SIRIUS_MV = -1.46
SIRIUS_RA_DEG = 101.09
SIRIUS_DEC_DEG = -16.69
OBSERVER_LATITUDE_DEG = 30.0
MZ = 6.0

PUBLISHED_TABLE3 = {
    0.20: 133.8,
    0.30: 136.3,
    0.40: 138.5,
}

# A half-day is a deliberately transparent tolerance for this validator,
# because the 1987 paper refers to Nautical Almanac solar equations that are
# not reproduced in the article, whereas this script uses the compact solar
# approximation below.
DEFAULT_TOLERANCE_DAYS = 0.50


@dataclass(frozen=True)
class DawnMaximum:
    """Best visibility point found during one trial dawn."""

    delta_days: float
    visibility_margin_mag: float
    sun_hour_angle_deg: float
    sun_altitude_deg: float
    sun_azimuth_deg: float
    star_hour_angle_deg: float
    star_altitude_deg: float
    star_azimuth_deg: float
    sun_star_azimuth_difference_deg: float
    limiting_magnitude_v: float
    airmass: float
    background_nl: float
    calibration_domain: bool


@dataclass(frozen=True)
class ValidationRow:
    kv: float
    published_delta_days: float
    calculated_delta_days: float
    residual_days: float
    absolute_residual_days: float
    passed: bool
    dawn: DawnMaximum


# ---------------------------------------------------------------------------
# Basic angle / spherical-astronomy helpers
# ---------------------------------------------------------------------------


def wrap180(angle_deg: float) -> float:
    """Wrap an angle to [-180, +180) degrees."""
    return (angle_deg + 180.0) % 360.0 - 180.0


def angular_difference_deg(a_deg: float, b_deg: float) -> float:
    """Smallest absolute angular separation of two azimuths [0, 180]."""
    return abs(wrap180(a_deg - b_deg))


def altitude_azimuth(
    latitude_deg: float,
    declination_deg: float,
    hour_angle_deg: float,
) -> tuple[float, float]:
    """Return geometric altitude and azimuth in degrees.

    Azimuth is measured clockwise from geographic north (N=0, E=90,
    S=180, W=270), matching the convention used in the Chiang Saen study.
    """
    phi = math.radians(latitude_deg)
    dec = math.radians(declination_deg)
    H = math.radians(hour_angle_deg)

    sin_h = (
        math.sin(phi) * math.sin(dec)
        + math.cos(phi) * math.cos(dec) * math.cos(H)
    )
    sin_h = max(-1.0, min(1.0, sin_h))
    altitude_deg = math.degrees(math.asin(sin_h))

    # atan2 form; the +180 shift converts to north=0, east=90 convention.
    azimuth_deg = (
        math.degrees(
            math.atan2(
                math.sin(H),
                math.cos(H) * math.sin(phi) - math.tan(dec) * math.cos(phi),
            )
        )
        + 180.0
    ) % 360.0

    return altitude_deg, azimuth_deg


# ---------------------------------------------------------------------------
# Compact solar ephemeris for the Table-3 validation only
# ---------------------------------------------------------------------------


def _solar_equation_of_center_deg(mean_anomaly_deg: float) -> float:
    """Low-order solar equation of center, adequate for this cross-check."""
    M = math.radians(mean_anomaly_deg)
    return (
        1.9148 * math.sin(M)
        + 0.0200 * math.sin(2.0 * M)
        + 0.0003 * math.sin(3.0 * M)
    )


def solar_radec_from_delta(delta_days: float) -> tuple[float, float]:
    """Approximate solar RA/Dec for Delta days after the vernal equinox.

    The orbital model is appropriate to the A.D. 1950 validation epoch.
    The mean anomaly at the 1950 vernal equinox is approximately 76 deg.
    We anchor ecliptic longitude to exactly 0 deg at Delta=0, so the
    validation is insensitive to the absolute calendar date of the equinox.
    """
    mean_daily_motion = 360.0 / 365.2422
    mean_anomaly_at_equinox = 76.0  # deg; appropriate to ca. A.D. 1950
    obliquity_deg = 23.445

    M0 = mean_anomaly_at_equinox
    longitude_zero = M0 + _solar_equation_of_center_deg(M0)

    M = M0 + mean_daily_motion * delta_days
    longitude_deg = (
        M + _solar_equation_of_center_deg(M) - longitude_zero
    ) % 360.0

    lam = math.radians(longitude_deg)
    eps = math.radians(obliquity_deg)

    ra_deg = math.degrees(
        math.atan2(math.cos(eps) * math.sin(lam), math.cos(lam))
    ) % 360.0
    dec_deg = math.degrees(math.asin(math.sin(eps) * math.sin(lam)))

    return ra_deg, dec_deg


# ---------------------------------------------------------------------------
# Dawn visibility search
# ---------------------------------------------------------------------------


def dawn_maximum_visibility(
    delta_days: float,
    kv: float,
    *,
    hour_angle_step_deg: float = 0.05,
) -> DawnMaximum:
    """Find the maximum visibility margin during dawn for one Delta.

    We scan the morning sky by solar hour angle.  Only instants with the Sun
    below the horizon and Sirius above the horizon are evaluated.  The broad
    solar-altitude range includes Schaefer's calibrated twilight domain and
    permits bright-star solutions close to sunrise, as discussed in the 1987
    paper.
    """
    if hour_angle_step_deg <= 0.0:
        raise ValueError("hour_angle_step_deg must be positive")

    sun_ra_deg, sun_dec_deg = solar_radec_from_delta(delta_days)

    best: Optional[DawnMaximum] = None

    # Dawn occurs for negative solar hour angle.  -150..0 deg safely brackets
    # the relevant pre-sunrise interval for latitude +30 deg.
    H_sun = -150.0
    while H_sun <= 0.0 + 1e-12:
        sun_alt_deg, sun_az_deg = altitude_azimuth(
            OBSERVER_LATITUDE_DEG, sun_dec_deg, H_sun
        )

        # Schaefer's heliacal application concerns the dawn twilight.  The
        # visibility routine itself flags whether a point lies inside the
        # optimized -18..-6 deg solar-altitude range.
        if -25.0 <= sun_alt_deg <= 0.0:
            H_star = wrap180(H_sun + sun_ra_deg - SIRIUS_RA_DEG)
            star_alt_deg, star_az_deg = altitude_azimuth(
                OBSERVER_LATITUDE_DEG, SIRIUS_DEC_DEG, H_star
            )

            if star_alt_deg >= 0.0:
                theta_deg = angular_difference_deg(star_az_deg, sun_az_deg)

                try:
                    result = sv.evaluate_visibility(
                        star_magnitude_v=SIRIUS_MV,
                        sun_altitude_deg=sun_alt_deg,
                        star_altitude_deg=star_alt_deg,
                        sun_star_azimuth_difference_deg=theta_deg,
                        kv=kv,
                        mz=MZ,
                    )
                except ValueError:
                    # Eq. (2) can yield a non-physical adjusted background
                    # when extrapolated too far outside the twilight range;
                    # those instants are not admissible model points.
                    H_sun += hour_angle_step_deg
                    continue

                candidate = DawnMaximum(
                    delta_days=delta_days,
                    visibility_margin_mag=result.visibility_margin_mag,
                    sun_hour_angle_deg=H_sun,
                    sun_altitude_deg=sun_alt_deg,
                    sun_azimuth_deg=sun_az_deg,
                    star_hour_angle_deg=H_star,
                    star_altitude_deg=star_alt_deg,
                    star_azimuth_deg=star_az_deg,
                    sun_star_azimuth_difference_deg=theta_deg,
                    limiting_magnitude_v=result.limiting_magnitude_v,
                    airmass=result.airmass,
                    background_nl=result.background_nl,
                    calibration_domain=result.within_schaefer1987_calibration_domain,
                )

                if (
                    best is None
                    or candidate.visibility_margin_mag > best.visibility_margin_mag
                ):
                    best = candidate

        H_sun += hour_angle_step_deg

    if best is None:
        raise RuntimeError(
            f"No valid dawn geometry found at Delta={delta_days:.3f} d"
        )

    return best


def first_visible_delta(
    kv: float,
    *,
    delta_lo: float = 126.0,
    delta_hi: float = 150.0,
    delta_tolerance: float = 1.0e-4,
    hour_angle_step_deg: float = 0.05,
) -> DawnMaximum:
    """Solve for the first Delta at which maximum dawn margin reaches zero."""
    lo = float(delta_lo)
    hi = float(delta_hi)

    lo_dawn = dawn_maximum_visibility(
        lo, kv, hour_angle_step_deg=hour_angle_step_deg
    )
    hi_dawn = dawn_maximum_visibility(
        hi, kv, hour_angle_step_deg=hour_angle_step_deg
    )

    if lo_dawn.visibility_margin_mag >= 0.0:
        raise ValueError(
            "Lower Delta bound is already visible; decrease delta_lo."
        )
    if hi_dawn.visibility_margin_mag < 0.0:
        raise ValueError(
            "Upper Delta bound is still invisible; increase delta_hi."
        )

    # Bisection on the maximum dawn visibility margin.
    while (hi - lo) > delta_tolerance:
        mid = 0.5 * (lo + hi)
        mid_dawn = dawn_maximum_visibility(
            mid, kv, hour_angle_step_deg=hour_angle_step_deg
        )
        if mid_dawn.visibility_margin_mag >= 0.0:
            hi = mid
        else:
            lo = mid

    return dawn_maximum_visibility(
        0.5 * (lo + hi), kv, hour_angle_step_deg=hour_angle_step_deg
    )


# ---------------------------------------------------------------------------
# Validation report / CSV
# ---------------------------------------------------------------------------


def run_validation(
    *,
    tolerance_days: float = DEFAULT_TOLERANCE_DAYS,
    hour_angle_step_deg: float = 0.05,
) -> list[ValidationRow]:
    """Run the three extinction cases printed in Schaefer (1987), Table 3."""
    rows: list[ValidationRow] = []

    for kv, published_delta in PUBLISHED_TABLE3.items():
        dawn = first_visible_delta(
            kv,
            hour_angle_step_deg=hour_angle_step_deg,
        )
        residual = dawn.delta_days - published_delta
        rows.append(
            ValidationRow(
                kv=kv,
                published_delta_days=published_delta,
                calculated_delta_days=dawn.delta_days,
                residual_days=residual,
                absolute_residual_days=abs(residual),
                passed=(abs(residual) <= tolerance_days),
                dawn=dawn,
            )
        )

    return rows


def write_csv(rows: Iterable[ValidationRow], path: Path) -> None:
    """Write reviewer-auditable validation diagnostics."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "kv",
                "published_delta_days",
                "calculated_delta_days",
                "residual_days",
                "absolute_residual_days",
                "passed",
                "sun_hour_angle_deg",
                "sun_altitude_deg",
                "sun_azimuth_deg",
                "star_hour_angle_deg",
                "star_altitude_deg",
                "star_azimuth_deg",
                "sun_star_azimuth_difference_deg",
                "limiting_magnitude_v",
                "visibility_margin_mag",
                "airmass",
                "background_nl",
                "within_schaefer1987_calibration_domain",
            ],
        )
        writer.writeheader()
        for row in rows:
            d = row.dawn
            writer.writerow(
                {
                    "kv": f"{row.kv:.2f}",
                    "published_delta_days": f"{row.published_delta_days:.3f}",
                    "calculated_delta_days": f"{row.calculated_delta_days:.6f}",
                    "residual_days": f"{row.residual_days:+.6f}",
                    "absolute_residual_days": f"{row.absolute_residual_days:.6f}",
                    "passed": row.passed,
                    "sun_hour_angle_deg": f"{d.sun_hour_angle_deg:.6f}",
                    "sun_altitude_deg": f"{d.sun_altitude_deg:.6f}",
                    "sun_azimuth_deg": f"{d.sun_azimuth_deg:.6f}",
                    "star_hour_angle_deg": f"{d.star_hour_angle_deg:.6f}",
                    "star_altitude_deg": f"{d.star_altitude_deg:.6f}",
                    "star_azimuth_deg": f"{d.star_azimuth_deg:.6f}",
                    "sun_star_azimuth_difference_deg": (
                        f"{d.sun_star_azimuth_difference_deg:.6f}"
                    ),
                    "limiting_magnitude_v": f"{d.limiting_magnitude_v:.6f}",
                    "visibility_margin_mag": f"{d.visibility_margin_mag:+.8f}",
                    "airmass": f"{d.airmass:.6f}",
                    "background_nl": f"{d.background_nl:.6f}",
                    "within_schaefer1987_calibration_domain": d.calibration_domain,
                }
            )


def print_report(rows: list[ValidationRow], tolerance_days: float) -> None:
    print("Schaefer (1987) Table 3 validation: Sirius heliacal rising")
    print("Parameters: m_V=-1.46, RA=101.09 deg, Dec=-16.69 deg, latitude=+30 deg")
    print(f"Acceptance tolerance: +/- {tolerance_days:.2f} day")
    print()
    print(
        " k_v   published Δ   calculated Δ   residual    star alt   Sun alt   status"
    )
    print(
        "-----  -----------   ------------   --------    --------   -------   ------"
    )

    for row in rows:
        d = row.dawn
        status = "PASS" if row.passed else "FAIL"
        print(
            f" {row.kv:0.2f}      {row.published_delta_days:7.1f}"
            f"        {row.calculated_delta_days:8.3f}"
            f"     {row.residual_days:+7.3f} d"
            f"     {d.star_altitude_deg:6.2f} deg"
            f"   {d.sun_altitude_deg:6.2f} deg"
            f"   {status}"
        )

    print()
    if all(row.passed for row in rows):
        print("RESULT: PASS -- all published k_v cases are reproduced within tolerance.")
    else:
        print("RESULT: FAIL -- one or more cases exceed the stated tolerance.")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate schaefer_visibility.py against Schaefer (1987) Table 3."
    )
    parser.add_argument(
        "--tolerance-days",
        type=float,
        default=DEFAULT_TOLERANCE_DAYS,
        help=f"PASS tolerance in days (default {DEFAULT_TOLERANCE_DAYS}).",
    )
    parser.add_argument(
        "--hour-angle-step",
        type=float,
        default=0.05,
        help="Dawn scan step in degrees of solar hour angle (default 0.05).",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=Path("validation_schaefer1987.csv"),
        help="Output CSV path (default validation_schaefer1987.csv).",
    )
    args = parser.parse_args()

    if args.tolerance_days <= 0.0:
        parser.error("--tolerance-days must be positive")
    if args.hour_angle_step <= 0.0:
        parser.error("--hour-angle-step must be positive")

    rows = run_validation(
        tolerance_days=args.tolerance_days,
        hour_angle_step_deg=args.hour_angle_step,
    )
    print_report(rows, args.tolerance_days)
    write_csv(rows, args.csv)
    print(f"Validation diagnostics written to: {args.csv}")

    return 0 if all(row.passed for row in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
