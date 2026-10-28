#!/usr/bin/env python3
"""End-to-end validation of heliacal_search.py against Schaefer (1987).

This validator reproduces the *standard Sirius heliacal-rising case* printed
in Schaefer (1987), Table 3, using the same visibility engine and day-by-day
dawn search used by the Chiang Saen analysis.

Published validation target
---------------------------
Schaefer (1987), pp. S24-S25, defines Delta as days elapsed since the vernal
equinox and gives for Sirius:

    epoch               A.D. 1950
    m_V                 -1.46
    RA (epoch 1950)     101.09 deg
    Dec (epoch 1950)    -16.69 deg
    observer latitude   +30 deg
    k_v                  0.20 mag/airmass
    m_z                  6.0 mag
    published Delta     133.8 days

The public heliacal_search.py interface expects a J2000 StarJ2000 object.
For this validation only, we numerically invert the same precession routine so
that the validation star has Schaefer's quoted RA/Dec at the 1950 vernal
equinox.  This avoids substituting a modern catalogue coordinate for the
coordinate explicitly printed in the validation paper.

The search is run with geometric altitudes (apply_refraction=False), because
Schaefer (1987) states that solar altitude h0 and stellar zenith distance Z*
are derived by spherical trigonometry from the celestial coordinates and
local sidereal time.

Why the tolerance is not zero
-----------------------------
This is an independent reconstruction, not execution of Schaefer's original
1985/1987 program.  Our solar ephemeris, precession implementation, numerical
sampling, and equinox calculation are independent.  A tolerance of +/-0.5 d
is therefore adopted for the end-to-end cross-check.  This is stricter than
one calendar day and appropriate to a first-visible-morning calculation.

No third-party packages are required.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import astronomy_geometry as ag
import heliacal_search as hs


# ---------------------------------------------------------------------------
# Schaefer (1987) published standard Sirius case
# ---------------------------------------------------------------------------

PUBLISHED_DELTA_DAYS = 133.8
SIRIUS_MV = -1.46
SIRIUS_RA_1950_DEG = 101.09
SIRIUS_DEC_1950_DEG = -16.69
VALIDATION_YEAR = 1950
VALIDATION_KV = 0.20
VALIDATION_MZ = 6.0
VALIDATION_LATITUDE_DEG = 30.0
DEFAULT_TOLERANCE_DAYS = 0.50

VALIDATION_SITE = ag.Site(
    name="Schaefer (1987) Sirius validation site",
    latitude_deg=VALIDATION_LATITUDE_DEG,
    longitude_deg=0.0,  # longitude is arbitrary when reporting LMST
    elevation_m=None,
)


@dataclass(frozen=True)
class ValidationResult:
    published_delta_days: float
    calculated_delta_days: float
    residual_days: float
    absolute_residual_days: float
    tolerance_days: float
    passed: bool
    equinox_jd: float
    first_visible_date: str
    first_visible_time_lmst: str
    maximum_visibility_time_lmst: str
    sun_altitude_at_max_deg: float
    star_altitude_at_max_deg: float
    star_zenith_distance_at_max_deg: float
    star_azimuth_at_max_deg: float
    visibility_margin_at_max_mag: float


# ---------------------------------------------------------------------------
# Vernal equinox and validation-coordinate adapter
# ---------------------------------------------------------------------------


def vernal_equinox_jd(year: int) -> float:
    """Find the March equinox as the Sun's declination crosses northward.

    The root uses the same compact solar ephemeris as astronomy_geometry.py,
    so the validation remains internally self-consistent.
    """
    lo = ag.julian_day(year, 3, 18, 0.0, calendar="gregorian")
    hi = ag.julian_day(year, 3, 23, 0.0, calendar="gregorian")

    dec_lo = ag.sun_equatorial_mean_deg(lo)[1]
    dec_hi = ag.sun_equatorial_mean_deg(hi)[1]
    if not (dec_lo < 0.0 < dec_hi):
        raise RuntimeError("Failed to bracket the March equinox.")

    for _ in range(70):
        mid = 0.5 * (lo + hi)
        dec_mid = ag.sun_equatorial_mean_deg(mid)[1]
        if dec_mid >= 0.0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def _wrap180(angle_deg: float) -> float:
    return (angle_deg + 180.0) % 360.0 - 180.0


def make_schaefer_sirius_j2000(reference_jd: float) -> ag.StarJ2000:
    """Construct a J2000 adapter reproducing Schaefer's epoch-1950 RA/Dec.

    heliacal_search.py accepts J2000 coordinates and precesses them to date.
    Schaefer instead prints the validation coordinate directly for epoch 1950.
    We therefore invert our own precession numerically. Proper motion is set
    to zero because the target itself is already the coordinate specified for
    the validation epoch.
    """
    ra0 = SIRIUS_RA_1950_DEG
    dec0 = SIRIUS_DEC_1950_DEG

    for _ in range(30):
        ra_date, dec_date = ag.precess_j2000_to_date_deg(ra0, dec0, reference_jd)
        dra = _wrap180(SIRIUS_RA_1950_DEG - ra_date)
        ddec = SIRIUS_DEC_1950_DEG - dec_date
        ra0 += dra
        dec0 += ddec
        if abs(dra) < 1e-12 and abs(ddec) < 1e-12:
            break

    # Verify the adapter itself before using it in the heliacal search.
    check_ra, check_dec = ag.precess_j2000_to_date_deg(ra0, dec0, reference_jd)
    if abs(_wrap180(check_ra - SIRIUS_RA_1950_DEG)) > 1e-8:
        raise RuntimeError("Failed to reproduce Schaefer's Sirius RA at epoch 1950.")
    if abs(check_dec - SIRIUS_DEC_1950_DEG) > 1e-8:
        raise RuntimeError("Failed to reproduce Schaefer's Sirius Dec at epoch 1950.")

    return ag.StarJ2000(
        name="Sirius — Schaefer (1987) standard validation case",
        ra_deg=ra0 % 360.0,
        dec_deg=dec0,
        pm_ra_cosdec_mas_per_yr=0.0,
        pm_dec_mas_per_yr=0.0,
    )


# ---------------------------------------------------------------------------
# End-to-end validation
# ---------------------------------------------------------------------------


def run_validation(
    *,
    tolerance_days: float = DEFAULT_TOLERANCE_DAYS,
    step_seconds: float = 2.0,
) -> ValidationResult:
    if tolerance_days <= 0.0:
        raise ValueError("tolerance_days must be positive")
    if step_seconds <= 0.0:
        raise ValueError("step_seconds must be positive")

    equinox = vernal_equinox_jd(VALIDATION_YEAR)
    sirius = make_schaefer_sirius_j2000(equinox)

    result = hs.search_heliacal_rising(
        sirius,
        SIRIUS_MV,
        VALIDATION_YEAR,
        VALIDATION_KV,
        mz=VALIDATION_MZ,
        site=VALIDATION_SITE,
        calendar="gregorian",
        step_seconds=step_seconds,
        max_days_after_conjunction=90,
        # Sirius is bright enough that Schaefer's own standard solution lies
        # outside the -18..-6 deg empirical-fit range quoted for Eq. (1).
        # The 1987 paper nevertheless evaluates this bright-star branch.
        require_calibration_domain=False,
        apply_refraction=False,
    )

    dawn = result.first_visible_morning
    if dawn.maximum is None or dawn.first_visible is None:
        raise RuntimeError("End-to-end search returned no visible point.")

    # Delta is defined by Schaefer as days since the vernal equinox.  Use the
    # maximum-visibility point on the first visible morning to attach a
    # continuous time coordinate to that first-visible date.
    pmax = dawn.maximum
    calculated_delta = pmax.jd_ut1 - equinox
    residual = calculated_delta - PUBLISHED_DELTA_DAYS

    return ValidationResult(
        published_delta_days=PUBLISHED_DELTA_DAYS,
        calculated_delta_days=calculated_delta,
        residual_days=residual,
        absolute_residual_days=abs(residual),
        tolerance_days=tolerance_days,
        passed=abs(residual) <= tolerance_days,
        equinox_jd=equinox,
        first_visible_date=f"{dawn.year:04d}-{dawn.month:02d}-{dawn.day:02d}",
        first_visible_time_lmst=dawn.first_visible.local_datetime.hhmmss,
        maximum_visibility_time_lmst=pmax.local_datetime.hhmmss,
        sun_altitude_at_max_deg=pmax.sun_altitude_deg,
        star_altitude_at_max_deg=pmax.star_altitude_deg,
        star_zenith_distance_at_max_deg=90.0 - pmax.star_altitude_deg,
        star_azimuth_at_max_deg=pmax.star_azimuth_deg,
        visibility_margin_at_max_mag=pmax.visibility_margin_mag,
    )


def write_csv(result: ValidationResult, path: Path) -> None:
    row = {
        "case": "Schaefer1987_Table3_Sirius_standard",
        "epoch": VALIDATION_YEAR,
        "latitude_deg": VALIDATION_LATITUDE_DEG,
        "m_v": SIRIUS_MV,
        "ra_epoch1950_deg": SIRIUS_RA_1950_DEG,
        "dec_epoch1950_deg": SIRIUS_DEC_1950_DEG,
        "k_v": VALIDATION_KV,
        "m_z": VALIDATION_MZ,
        "published_delta_days": result.published_delta_days,
        "calculated_delta_days": result.calculated_delta_days,
        "residual_days": result.residual_days,
        "absolute_residual_days": result.absolute_residual_days,
        "tolerance_days": result.tolerance_days,
        "passed": result.passed,
        "first_visible_date": result.first_visible_date,
        "first_visible_time_lmst": result.first_visible_time_lmst,
        "maximum_visibility_time_lmst": result.maximum_visibility_time_lmst,
        "sun_altitude_at_max_deg": result.sun_altitude_at_max_deg,
        "star_altitude_at_max_deg": result.star_altitude_at_max_deg,
        "star_zenith_distance_at_max_deg": result.star_zenith_distance_at_max_deg,
        "star_azimuth_at_max_deg": result.star_azimuth_at_max_deg,
        "visibility_margin_at_max_mag": result.visibility_margin_at_max_mag,
    }

    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(row.keys()))
        writer.writeheader()
        writer.writerow(row)


if __name__ == "__main__":
    out = Path(__file__).with_name("validation_heliacal_search_schaefer1987.csv")
    r = run_validation()
    write_csv(r, out)

    print("End-to-end heliacal_search.py validation")
    print("Case: Schaefer (1987), Table 3, standard Sirius heliacal rising")
    print(f"Published Delta:   {r.published_delta_days:.3f} d after vernal equinox")
    print(f"Calculated Delta:  {r.calculated_delta_days:.3f} d after vernal equinox")
    print(f"Residual:          {r.residual_days:+.3f} d")
    print(f"Tolerance:         +/-{r.tolerance_days:.3f} d")
    print(f"Status:            {'PASS' if r.passed else 'FAIL'}")
    print(f"First visible:     {r.first_visible_date} {r.first_visible_time_lmst} LMST")
    print(f"Maximum margin at: {r.maximum_visibility_time_lmst} LMST")
    print(f"Sun altitude:      {r.sun_altitude_at_max_deg:+.3f} deg")
    print(f"Sirius altitude:   {r.star_altitude_at_max_deg:.3f} deg")
    print(f"Sirius Z*:         {r.star_zenith_distance_at_max_deg:.3f} deg")
    print(f"Sirius azimuth:    {r.star_azimuth_at_max_deg:.3f} deg")
    print(f"CSV:               {out}")

    raise SystemExit(0 if r.passed else 1)
