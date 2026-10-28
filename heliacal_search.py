#!/usr/bin/env python3
"""Search for Schaefer-style heliacal rising at Chiang Saen.

This module connects the two independently testable layers used in the
Chiang Saen revision:

    astronomy_geometry.py  -> Sun/star position at each instant
    schaefer_visibility.py -> naked-eye visibility at that instant

The heliacal-rising search follows the operational definition used by
Schaefer (1987, 1993): the heliacal rising is the *first morning after solar
conjunction on which the star is visible at some time during dawn twilight*.

Unlike the manuscript's superseded fixed-Sun-altitude method, the code does
NOT set the Sun to -4 deg.  It scans dawn, evaluates atmospheric extinction,
twilight brightness, and the human-eye threshold at each step, and then
searches day by day for the first morning containing a visible interval.

Project conventions
-------------------
* Site: Pa Sak Gate, Chiang Saen, 20.2736 N, 100.0786 E (WGS 84).
* Historical calendar labels: Julian calendar for epochs before the Gregorian
  reform, matching the historical-date convention used in the Stellarium
  sanity check.  The calendar is always explicit in the public API.
* Time labels: local mean solar time (LMST), not modern UTC+7 civil time.
* Azimuth: clockwise from geographic (true) north.
* Chiang Saen extinction scenarios: k_v = 0.20 and 0.30 mag/airmass.
* Default dark-sky zenith limiting magnitude: m_z = 6.0.

Calibration-domain policy
-------------------------
Schaefer (1987) states that the empirical twilight relation used by
``schaefer_visibility.py`` was optimized approximately for

    -18 deg < Sun altitude < -6 deg,
    stellar zenith distance > 60 deg (star altitude < 30 deg),
    Sun-star azimuth difference < 90 deg.

The code always reports whether a point lies inside that quoted empirical-fit
domain, but it does not reject out-of-domain points by default.  This follows
Schaefer's own 1987 treatment: his bright-star examples (notably Sirius) can
reach first visibility with the Sun higher than -6 deg.  A strict-domain mode
remains available as an optional diagnostic sensitivity test.

No third-party packages are required.
"""

from __future__ import annotations

from dataclasses import dataclass
import argparse
import math
from typing import Iterable, Optional, Sequence, Tuple

import astronomy_geometry as ag
import schaefer_visibility as sv


# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class LocalMeanSolarDateTime:
    """Calendar date and local mean solar clock time."""

    year: int
    month: int
    day: int
    hours: float

    @property
    def hhmmss(self) -> str:
        total_seconds = int(round(self.hours * 3600.0)) % 86400
        hh = total_seconds // 3600
        mm = (total_seconds % 3600) // 60
        ss = total_seconds % 60
        return f"{hh:02d}:{mm:02d}:{ss:02d}"

    @property
    def ymd(self) -> str:
        return f"{self.year:04d}-{self.month:02d}-{self.day:02d}"


@dataclass(frozen=True)
class VisibilityPoint:
    """Geometry + visibility diagnostics for one instant."""

    jd_ut1: float
    local_datetime: LocalMeanSolarDateTime
    star_name: str
    kv: float

    sun_altitude_deg: float
    sun_azimuth_deg: float
    star_altitude_deg: float
    star_azimuth_deg: float
    sun_star_azimuth_difference_deg: float

    visible: bool
    visibility_margin_mag: float
    limiting_magnitude_v: float
    airmass: float
    extinction_mag: float
    background_nl: float
    within_calibration_domain: bool


@dataclass(frozen=True)
class DawnScan:
    """Visibility summary for one local morning."""

    year: int
    month: int
    day: int
    kv: float
    any_visible: bool
    maximum: Optional[VisibilityPoint]
    first_visible: Optional[VisibilityPoint]
    last_visible: Optional[VisibilityPoint]
    evaluated_points: int


@dataclass(frozen=True)
class ConjunctionResult:
    """Approximate annual solar conjunction / minimum angular separation."""

    jd_ut1: float
    local_datetime: LocalMeanSolarDateTime
    angular_separation_deg: float


@dataclass(frozen=True)
class HeliacalRisingResult:
    """First visible morning after conjunction for one k_v scenario."""

    star_name: str
    star_magnitude_v: float
    nominal_epoch_year: int
    calendar: str
    kv: float
    mz: float
    conjunction: ConjunctionResult
    first_visible_morning: DawnScan
    days_after_conjunction: float


# ---------------------------------------------------------------------------
# Calendar / local mean solar time helpers
# ---------------------------------------------------------------------------


def _jd_to_calendar(jd: float, *, calendar: str) -> Tuple[int, int, int, float]:
    """Convert JD to calendar date + decimal hours.

    The conversion is the standard inverse Julian-day algorithm.  ``calendar``
    is explicit so historical labels are never silently switched between the
    Julian and Gregorian calendars.
    """
    if calendar not in {"julian", "gregorian"}:
        raise ValueError("calendar must be 'julian' or 'gregorian'")

    x = jd + 0.5
    Z = math.floor(x)
    F = x - Z

    if calendar == "gregorian":
        alpha = math.floor((Z - 1867216.25) / 36524.25)
        A = Z + 1 + alpha - math.floor(alpha / 4)
    else:
        A = Z

    B = A + 1524
    C = math.floor((B - 122.1) / 365.25)
    D = math.floor(365.25 * C)
    E = math.floor((B - D) / 30.6001)

    day_decimal = B - D - math.floor(30.6001 * E) + F
    day = math.floor(day_decimal)
    frac_day = day_decimal - day

    if E < 14:
        month = E - 1
    else:
        month = E - 13

    if month > 2:
        year = C - 4716
    else:
        year = C - 4715

    hours = frac_day * 24.0
    return int(year), int(month), int(day), hours


def local_mean_solar_datetime(
    jd_ut1: float,
    site: ag.Site = ag.PA_SAK_GATE,
    *,
    calendar: str = "julian",
) -> LocalMeanSolarDateTime:
    """Convert JD(UT1-like) to a local-mean-solar date/time label."""
    local_jd = jd_ut1 + site.longitude_deg / 360.0
    y, m, d, h = _jd_to_calendar(local_jd, calendar=calendar)
    return LocalMeanSolarDateTime(y, m, d, h)


def _local_midnight_jd(
    year: int,
    month: int,
    day: int,
    site: ag.Site,
    *,
    calendar: str,
) -> float:
    return ag.jd_for_local_mean_solar_time(
        year, month, day, 0.0, site, calendar=calendar
    )


# ---------------------------------------------------------------------------
# Angular separation / conjunction search
# ---------------------------------------------------------------------------


def _great_circle_separation_deg(
    ra1_deg: float, dec1_deg: float, ra2_deg: float, dec2_deg: float
) -> float:
    a1 = math.radians(ra1_deg)
    d1 = math.radians(dec1_deg)
    a2 = math.radians(ra2_deg)
    d2 = math.radians(dec2_deg)

    c = (
        math.sin(d1) * math.sin(d2)
        + math.cos(d1) * math.cos(d2) * math.cos(a1 - a2)
    )
    c = max(-1.0, min(1.0, c))
    return math.degrees(math.acos(c))


def solar_star_separation_deg(star: ag.StarJ2000, jd_ut1: float) -> float:
    """Geocentric angular separation between Sun and star."""
    star_ra, star_dec = ag.star_equatorial_of_date_deg(star, jd_ut1)
    sun_ra, sun_dec = ag.sun_equatorial_mean_deg(jd_ut1)
    return _great_circle_separation_deg(star_ra, star_dec, sun_ra, sun_dec)


def find_annual_conjunction(
    star: ag.StarJ2000,
    year: int,
    site: ag.Site = ag.PA_SAK_GATE,
    *,
    calendar: str = "julian",
    coarse_step_hours: float = 6.0,
    refine_tolerance_seconds: float = 2.0,
) -> ConjunctionResult:
    """Find minimum Sun-star angular separation during a labelled year.

    A coarse full-year scan is followed by a golden-section minimization over
    the interval surrounding the best coarse point.  For Nakshatra stars this
    minimum is the annual solar conjunction used as the start of the heliacal
    rising search.
    """
    if coarse_step_hours <= 0.0:
        raise ValueError("coarse_step_hours must be positive")
    if refine_tolerance_seconds <= 0.0:
        raise ValueError("refine_tolerance_seconds must be positive")

    start = _local_midnight_jd(year, 1, 1, site, calendar=calendar)
    end = _local_midnight_jd(year + 1, 1, 1, site, calendar=calendar)
    step = coarse_step_hours / 24.0

    best_jd = start
    best_sep = solar_star_separation_deg(star, start)
    jd = start + step
    while jd < end:
        sep = solar_star_separation_deg(star, jd)
        if sep < best_sep:
            best_sep = sep
            best_jd = jd
        jd += step

    # Bracket the local minimum by one coarse step either side.
    lo = max(start, best_jd - step)
    hi = min(end, best_jd + step)

    # Golden-section minimization.
    gr = (math.sqrt(5.0) - 1.0) / 2.0
    c = hi - gr * (hi - lo)
    d = lo + gr * (hi - lo)
    fc = solar_star_separation_deg(star, c)
    fd = solar_star_separation_deg(star, d)
    tol_days = refine_tolerance_seconds / 86400.0

    while (hi - lo) > tol_days:
        if fc < fd:
            hi, d, fd = d, c, fc
            c = hi - gr * (hi - lo)
            fc = solar_star_separation_deg(star, c)
        else:
            lo, c, fc = c, d, fd
            d = lo + gr * (hi - lo)
            fd = solar_star_separation_deg(star, d)

    jd_min = 0.5 * (lo + hi)
    sep_min = solar_star_separation_deg(star, jd_min)
    return ConjunctionResult(
        jd_ut1=jd_min,
        local_datetime=local_mean_solar_datetime(
            jd_min, site, calendar=calendar
        ),
        angular_separation_deg=sep_min,
    )


# ---------------------------------------------------------------------------
# One-instant and one-dawn visibility evaluation
# ---------------------------------------------------------------------------


def evaluate_instant(
    star: ag.StarJ2000,
    star_magnitude_v: float,
    jd_ut1: float,
    kv: float,
    *,
    mz: float = sv.DEFAULT_MZ,
    site: ag.Site = ag.PA_SAK_GATE,
    calendar: str = "julian",
    pressure_hpa: float = 1010.0,
    temperature_c: float = 20.0,
    apply_refraction: bool = False,
) -> Optional[VisibilityPoint]:
    """Evaluate one instant, returning None when geometry is inadmissible.

    By default, geometric Sun/star altitudes are used, matching Schaefer
    (1987), who derives h0 and Z* from spherical trigonometry.  Optional
    refraction can be enabled explicitly for sensitivity tests.
    """
    g = ag.geometry_snapshot(
        star,
        jd_ut1,
        site,
        apply_refraction=apply_refraction,
        pressure_hpa=pressure_hpa,
        temperature_c=temperature_c,
    )

    if apply_refraction:
        sun_alt = g.sun_altitude_apparent_deg
        star_alt = g.star.altitude_apparent_deg
    else:
        sun_alt = g.sun_altitude_geometric_deg
        star_alt = g.star.altitude_geometric_deg

    if not (-90.0 <= sun_alt <= 0.0):
        return None
    if not (0.0 <= star_alt <= 90.0):
        return None

    try:
        v = sv.evaluate_visibility(
            star_magnitude_v=star_magnitude_v,
            sun_altitude_deg=sun_alt,
            star_altitude_deg=star_alt,
            sun_star_azimuth_difference_deg=g.sun_star_azimuth_difference_deg,
            kv=kv,
            mz=mz,
        )
    except ValueError:
        # The empirical background adjustment can become non-physical when
        # extrapolated far outside its useful twilight domain.
        return None

    return VisibilityPoint(
        jd_ut1=jd_ut1,
        local_datetime=local_mean_solar_datetime(
            jd_ut1, site, calendar=calendar
        ),
        star_name=star.name,
        kv=kv,
        sun_altitude_deg=sun_alt,
        sun_azimuth_deg=g.sun_azimuth_deg,
        star_altitude_deg=star_alt,
        star_azimuth_deg=g.star.azimuth_deg,
        sun_star_azimuth_difference_deg=g.sun_star_azimuth_difference_deg,
        visible=v.visible,
        visibility_margin_mag=v.visibility_margin_mag,
        limiting_magnitude_v=v.limiting_magnitude_v,
        airmass=v.airmass,
        extinction_mag=v.extinction_mag,
        background_nl=v.background_nl,
        within_calibration_domain=v.within_schaefer1987_calibration_domain,
    )


def scan_dawn(
    star: ag.StarJ2000,
    star_magnitude_v: float,
    year: int,
    month: int,
    day: int,
    kv: float,
    *,
    mz: float = sv.DEFAULT_MZ,
    site: ag.Site = ag.PA_SAK_GATE,
    calendar: str = "julian",
    lmst_start_hours: float = 2.0,
    lmst_end_hours: float = 8.0,
    step_seconds: float = 30.0,
    require_calibration_domain: bool = False,
    pressure_hpa: float = 1010.0,
    temperature_c: float = 20.0,
    apply_refraction: bool = False,
) -> DawnScan:
    """Scan one dawn and summarize any visible interval.

    ``maximum`` is the evaluated point with the greatest visibility margin.
    ``first_visible`` and ``last_visible`` are the first/last sampled points
    with margin >= 0.  The reported threshold times are therefore quantized by
    ``step_seconds``; the heliacal *date* is not affected once the visible
    interval is sampled.  A later reporting layer may refine threshold times.
    """
    if step_seconds <= 0.0:
        raise ValueError("step_seconds must be positive")
    if lmst_end_hours <= lmst_start_hours:
        raise ValueError("lmst_end_hours must exceed lmst_start_hours")

    first: Optional[VisibilityPoint] = None
    last: Optional[VisibilityPoint] = None
    maximum: Optional[VisibilityPoint] = None
    evaluated = 0

    n = int(math.floor((lmst_end_hours - lmst_start_hours) * 3600.0 / step_seconds))
    for i in range(n + 1):
        local_hours = lmst_start_hours + i * step_seconds / 3600.0
        jd = ag.jd_for_local_mean_solar_time(
            year, month, day, local_hours, site, calendar=calendar
        )
        p = evaluate_instant(
            star,
            star_magnitude_v,
            jd,
            kv,
            mz=mz,
            site=site,
            calendar=calendar,
            pressure_hpa=pressure_hpa,
            temperature_c=temperature_c,
            apply_refraction=apply_refraction,
        )
        if p is None:
            continue
        if require_calibration_domain and not p.within_calibration_domain:
            continue

        evaluated += 1
        if maximum is None or p.visibility_margin_mag > maximum.visibility_margin_mag:
            maximum = p

        if p.visible:
            if first is None:
                first = p
            last = p

    return DawnScan(
        year=year,
        month=month,
        day=day,
        kv=kv,
        any_visible=(first is not None),
        maximum=maximum,
        first_visible=first,
        last_visible=last,
        evaluated_points=evaluated,
    )


# ---------------------------------------------------------------------------
# Day-by-day heliacal-rising search
# ---------------------------------------------------------------------------


def _local_date_from_jd(
    jd_ut1: float,
    site: ag.Site,
    *,
    calendar: str,
) -> Tuple[int, int, int]:
    dt = local_mean_solar_datetime(jd_ut1, site, calendar=calendar)
    return dt.year, dt.month, dt.day


def search_heliacal_rising(
    star: ag.StarJ2000,
    star_magnitude_v: float,
    nominal_epoch_year: int,
    kv: float,
    *,
    mz: float = sv.DEFAULT_MZ,
    site: ag.Site = ag.PA_SAK_GATE,
    calendar: str = "julian",
    step_seconds: float = 30.0,
    max_days_after_conjunction: int = 90,
    require_calibration_domain: bool = False,
    pressure_hpa: float = 1010.0,
    temperature_c: float = 20.0,
    apply_refraction: bool = False,
) -> HeliacalRisingResult:
    """Find the first visible dawn after annual solar conjunction.

    The conjunction is found within ``nominal_epoch_year``.  The search then
    proceeds forward one local-mean-solar calendar day at a time.  The first
    morning containing at least one visible point is returned.
    """
    conjunction = find_annual_conjunction(
        star, nominal_epoch_year, site, calendar=calendar
    )

    # Start at the local midnight containing/just following conjunction.
    cy, cm, cd = _local_date_from_jd(
        conjunction.jd_ut1, site, calendar=calendar
    )
    day0 = _local_midnight_jd(cy, cm, cd, site, calendar=calendar)

    for offset in range(max_days_after_conjunction + 1):
        local_midnight = day0 + offset
        y, m, d = _local_date_from_jd(
            local_midnight + 1.0e-7, site, calendar=calendar
        )

        dawn = scan_dawn(
            star,
            star_magnitude_v,
            y,
            m,
            d,
            kv,
            mz=mz,
            site=site,
            calendar=calendar,
            step_seconds=step_seconds,
            require_calibration_domain=require_calibration_domain,
            pressure_hpa=pressure_hpa,
            temperature_c=temperature_c,
            apply_refraction=apply_refraction,
        )
        if dawn.any_visible:
            assert dawn.maximum is not None
            return HeliacalRisingResult(
                star_name=star.name,
                star_magnitude_v=star_magnitude_v,
                nominal_epoch_year=nominal_epoch_year,
                calendar=calendar,
                kv=kv,
                mz=mz,
                conjunction=conjunction,
                first_visible_morning=dawn,
                days_after_conjunction=dawn.maximum.jd_ut1 - conjunction.jd_ut1,
            )

    raise RuntimeError(
        f"No heliacal rising found for {star.name!r} within "
        f"{max_days_after_conjunction} days after conjunction at k_v={kv:.2f}."
    )


def search_project_kv_scenarios(
    star: ag.StarJ2000,
    star_magnitude_v: float,
    nominal_epoch_year: int,
    *,
    kv_values: Sequence[float] = sv.PROJECT_KV_SCENARIOS,
    **kwargs,
) -> Tuple[HeliacalRisingResult, ...]:
    """Run the Chiang Saen project's k_v scenarios for one star/epoch."""
    return tuple(
        search_heliacal_rising(
            star,
            star_magnitude_v,
            nominal_epoch_year,
            kv,
            **kwargs,
        )
        for kv in kv_values
    )


# ---------------------------------------------------------------------------
# Reporting helpers
# ---------------------------------------------------------------------------


def _format_point(label: str, p: Optional[VisibilityPoint]) -> str:
    if p is None:
        return f"  {label}: none"
    dt = p.local_datetime
    return (
        f"  {label}: {dt.ymd} {dt.hhmmss} LMST | "
        f"Sun alt={p.sun_altitude_deg:+.3f} deg | "
        f"star alt={p.star_altitude_deg:.3f} deg | "
        f"az={p.star_azimuth_deg:.3f} deg | "
        f"margin={p.visibility_margin_mag:+.3f} mag"
    )


def format_result(result: HeliacalRisingResult) -> str:
    """Human-readable audit summary."""
    c = result.conjunction
    d = result.first_visible_morning
    lines = [
        f"{result.star_name} | nominal epoch AD {result.nominal_epoch_year} | "
        f"k_v={result.kv:.2f}",
        (
            f"  conjunction/minimum separation: {c.local_datetime.ymd} "
            f"{c.local_datetime.hhmmss} LMST | separation={c.angular_separation_deg:.3f} deg"
        ),
        f"  first visible morning: {d.year:04d}-{d.month:02d}-{d.day:02d}",
        f"  days after conjunction (to dawn maximum): {result.days_after_conjunction:.3f}",
        _format_point("first sampled visible", d.first_visible),
        _format_point("maximum visibility", d.maximum),
        _format_point("last sampled visible", d.last_visible),
    ]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Demonstration / internal checks
# ---------------------------------------------------------------------------


# Alcyone (eta Tauri), ICRS/J2000 values used for the AD 500 geometry sanity
# check.  The final project catalogue will keep source metadata in stars.csv;
# this built-in object is only a reproducible demonstration target.
ALCYONE = ag.StarJ2000(
    name="Krittika / Alcyone (eta Tauri)",
    ra_deg=ag.ra_hms_to_deg(3, 47, 29.0766),
    dec_deg=ag.dec_dms_to_deg(24, 6, 18.49),
    pm_ra_cosdec_mas_per_yr=19.34,
    pm_dec_mas_per_yr=-43.67,
)
ALCYONE_V = 2.87


def _self_test() -> None:
    # Calendar round-trip at the exact historical convention used by the
    # previous Stellarium sanity check.
    jd = ag.jd_for_local_mean_solar_time(
        500, 5, 6, 5 + 10 / 60 + 24 / 3600,
        ag.PA_SAK_GATE,
        calendar="julian",
    )
    dt = local_mean_solar_datetime(jd, ag.PA_SAK_GATE, calendar="julian")
    assert (dt.year, dt.month, dt.day) == (500, 5, 6)
    assert abs(dt.hours - (5 + 10 / 60 + 24 / 3600)) < 2e-7

    # Great-circle identity.
    assert _great_circle_separation_deg(10.0, 20.0, 10.0, 20.0) < 1e-6


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Schaefer-style heliacal rising search for Chiang Saen"
    )
    parser.add_argument(
        "--demo-alcyone-ad500",
        action="store_true",
        help="run Alcyone AD 500 for k_v=0.20 and 0.30",
    )
    parser.add_argument(
        "--step-seconds",
        type=float,
        default=30.0,
        help="dawn sampling interval for the demo (default 30 s)",
    )
    args = parser.parse_args()

    _self_test()
    print("heliacal_search.py: internal checks passed")

    if args.demo_alcyone_ad500:
        for r in search_project_kv_scenarios(
            ALCYONE,
            ALCYONE_V,
            500,
            calendar="julian",
            step_seconds=args.step_seconds,
        ):
            print()
            print(format_result(r))

# ---------------------------------------------------------------------------
# Continuous maximum-margin dawn search
# ---------------------------------------------------------------------------


def _visibility_point_at_lmst_hours(
    star: ag.StarJ2000,
    star_magnitude_v: float,
    year: int,
    month: int,
    day: int,
    local_hours: float,
    kv: float,
    *,
    mz: float = sv.DEFAULT_MZ,
    site: ag.Site = ag.PA_SAK_GATE,
    calendar: str = "julian",
    pressure_hpa: float = 1010.0,
    temperature_c: float = 20.0,
    apply_refraction: bool = False,
) -> Optional[VisibilityPoint]:
    jd = ag.jd_for_local_mean_solar_time(
        year, month, day, local_hours, site, calendar=calendar
    )
    return evaluate_instant(
        star, star_magnitude_v, jd, kv,
        mz=mz, site=site, calendar=calendar,
        pressure_hpa=pressure_hpa, temperature_c=temperature_c,
        apply_refraction=apply_refraction,
    )


def _golden_maximize_visibility_margin(
    eval_func,
    a_hours: float,
    b_hours: float,
    *,
    tolerance_seconds: float = 0.25,
    max_iter: int = 80,
) -> Optional[VisibilityPoint]:
    """Maximize visibility margin on a local-time bracket.

    The bracket should surround a smooth local maximum.  Discontinuities at
    the Hecht day/night branch boundary are handled separately by
    ``maximize_dawn_visibility``.
    """
    if b_hours <= a_hours:
        return None
    phi = (math.sqrt(5.0) - 1.0) / 2.0
    c = b_hours - phi * (b_hours - a_hours)
    d = a_hours + phi * (b_hours - a_hours)
    pc = eval_func(c)
    pd = eval_func(d)
    def val(p):
        return -1.0e99 if p is None else p.visibility_margin_mag
    tol_h = tolerance_seconds / 3600.0
    it = 0
    while (b_hours - a_hours) > tol_h and it < max_iter:
        if val(pc) > val(pd):
            b_hours, d, pd = d, c, pc
            c = b_hours - phi * (b_hours - a_hours)
            pc = eval_func(c)
        else:
            a_hours, c, pc = c, d, pd
            d = a_hours + phi * (b_hours - a_hours)
            pd = eval_func(d)
        it += 1
    candidates = [p for p in (pc, pd, eval_func((a_hours+b_hours)/2.0)) if p is not None]
    return max(candidates, key=lambda p: p.visibility_margin_mag) if candidates else None


def maximize_dawn_visibility(
    star: ag.StarJ2000,
    star_magnitude_v: float,
    year: int,
    month: int,
    day: int,
    kv: float,
    *,
    mz: float = sv.DEFAULT_MZ,
    site: ag.Site = ag.PA_SAK_GATE,
    calendar: str = "julian",
    lmst_start_hours: float = 2.0,
    lmst_end_hours: float = 8.0,
    coarse_step_seconds: float = 120.0,
    optimize_tolerance_seconds: float = 0.25,
    require_calibration_domain: bool = False,
    pressure_hpa: float = 1010.0,
    temperature_c: float = 20.0,
    apply_refraction: bool = False,
) -> Optional[VisibilityPoint]:
    """Return the continuous-time maximum visibility point in one dawn.

    A coarse grid locates smooth local maxima.  Each is refined by a
    golden-section search.  In addition, every crossing of the Hecht
    point-source branch boundary (log10 B = 3.17) is explicitly located by
    bisection and sampled on both sides.  This matters because the empirical
    Hecht fit is piecewise and can create a sharp cusp near the photopic /
    scotopic transition; fixed time grids can otherwise shift the predicted
    first-visible date by one or more days.
    """
    if coarse_step_seconds <= 0:
        raise ValueError("coarse_step_seconds must be positive")
    if lmst_end_hours <= lmst_start_hours:
        raise ValueError("lmst_end_hours must exceed lmst_start_hours")

    def ev(h):
        p = _visibility_point_at_lmst_hours(
            star, star_magnitude_v, year, month, day, h, kv,
            mz=mz, site=site, calendar=calendar,
            pressure_hpa=pressure_hpa, temperature_c=temperature_c,
            apply_refraction=apply_refraction,
        )
        if p is not None and require_calibration_domain and not p.within_calibration_domain:
            return None
        return p

    step_h = coarse_step_seconds / 3600.0
    n = int(math.floor((lmst_end_hours-lmst_start_hours)/step_h))
    times = [lmst_start_hours + i*step_h for i in range(n+1)]
    if times[-1] < lmst_end_hours - 1e-12:
        times.append(lmst_end_hours)
    pts = [ev(t) for t in times]
    candidates = [p for p in pts if p is not None]
    if not candidates:
        return None

    def margin(i):
        p=pts[i]
        return -1.0e99 if p is None else p.visibility_margin_mag

    # Refine all coarse local maxima.  Also include the globally best coarse
    # point even if it sits adjacent to an invalid boundary.
    local_indices=[]
    for i in range(len(times)):
        if pts[i] is None:
            continue
        left = margin(i-1) if i>0 else -1.0e99
        right = margin(i+1) if i+1<len(times) else -1.0e99
        if margin(i) >= left and margin(i) >= right:
            local_indices.append(i)
    best_i=max((i for i,p in enumerate(pts) if p is not None), key=lambda i: margin(i))
    if best_i not in local_indices:
        local_indices.append(best_i)

    for i in local_indices:
        a=times[max(0,i-1)]; b=times[min(len(times)-1,i+1)]
        q=_golden_maximize_visibility_margin(ev,a,b,tolerance_seconds=optimize_tolerance_seconds)
        if q is not None:
            candidates.append(q)

    # Explicitly locate Hecht branch crossings.  Use log10(background)-3.17.
    break_log = sv.HECHT_LOG_B_BREAK
    for i in range(len(times)-1):
        p1,p2=pts[i],pts[i+1]
        if p1 is None or p2 is None or p1.background_nl<=0 or p2.background_nl<=0:
            continue
        f1=math.log10(p1.background_nl)-break_log
        f2=math.log10(p2.background_nl)-break_log
        if f1 == 0.0:
            root=times[i]
        elif f1*f2>0:
            continue
        else:
            lo,hi=times[i],times[i+1]
            flo,fhi=f1,f2
            for _ in range(50):
                if (hi-lo)*3600.0 <= 0.05:
                    break
                mid=(lo+hi)/2.0
                pm=ev(mid)
                if pm is None or pm.background_nl<=0:
                    break
                fm=math.log10(pm.background_nl)-break_log
                if flo*fm <= 0:
                    hi,fhi=mid,fm
                else:
                    lo,flo=mid,fm
            root=(lo+hi)/2.0
        # Evaluate immediately on both sides so a one-sided cusp is retained.
        eps_h=0.02/3600.0
        for h in (root-eps_h, root, root+eps_h):
            if lmst_start_hours <= h <= lmst_end_hours:
                q=ev(h)
                if q is not None:
                    candidates.append(q)

    return max(candidates,key=lambda p:p.visibility_margin_mag)


def search_heliacal_rising_continuous(
    star: ag.StarJ2000,
    star_magnitude_v: float,
    nominal_epoch_year: int,
    kv: float,
    *,
    mz: float = sv.DEFAULT_MZ,
    site: ag.Site = ag.PA_SAK_GATE,
    calendar: str = "julian",
    coarse_step_seconds: float = 120.0,
    max_days_after_conjunction: int = 120,
    require_calibration_domain: bool = False,
    pressure_hpa: float = 1010.0,
    temperature_c: float = 20.0,
    apply_refraction: bool = False,
    conjunction: Optional[ConjunctionResult] = None,
) -> HeliacalRisingResult:
    """Find first visible *date* using continuous maximum-margin optimization.

    This is the preferred production date search because it is insensitive to
    arbitrary sampling phase near the piecewise Hecht threshold cusp.
    ``DawnScan.maximum`` is the optimized maximum-margin point.  The
    first/last sampled threshold points are retained from a 5-s diagnostic
    scan of the selected date when available; they are not used to choose the
    heliacal date.
    """
    if conjunction is None:
        conjunction=find_annual_conjunction(star,nominal_epoch_year,site,calendar=calendar)
    cy,cm,cd=_local_date_from_jd(conjunction.jd_ut1,site,calendar=calendar)
    day0=_local_midnight_jd(cy,cm,cd,site,calendar=calendar)

    for offset in range(max_days_after_conjunction+1):
        midnight=day0+offset
        y,m,d=_local_date_from_jd(midnight+1e-7,site,calendar=calendar)
        pmax=maximize_dawn_visibility(
            star,star_magnitude_v,y,m,d,kv,mz=mz,site=site,calendar=calendar,
            coarse_step_seconds=coarse_step_seconds,
            require_calibration_domain=require_calibration_domain,
            pressure_hpa=pressure_hpa,temperature_c=temperature_c,
            apply_refraction=apply_refraction,
        )
        if pmax is None or pmax.visibility_margin_mag < 0.0:
            continue
        diag=scan_dawn(
            star,star_magnitude_v,y,m,d,kv,mz=mz,site=site,calendar=calendar,
            step_seconds=5.0,require_calibration_domain=require_calibration_domain,
            pressure_hpa=pressure_hpa,temperature_c=temperature_c,
            apply_refraction=apply_refraction,
        )
        # If a mathematically valid interval is narrower than the 5-s
        # diagnostic grid, preserve the optimized point as first/last rather
        # than falsely reporting no visibility.
        first=diag.first_visible if diag.first_visible is not None else pmax
        last=diag.last_visible if diag.last_visible is not None else pmax
        dawn=DawnScan(
            year=y,month=m,day=d,kv=kv,any_visible=True,
            maximum=pmax,first_visible=first,last_visible=last,
            evaluated_points=diag.evaluated_points,
        )
        return HeliacalRisingResult(
            star_name=star.name,star_magnitude_v=star_magnitude_v,
            nominal_epoch_year=nominal_epoch_year,calendar=calendar,kv=kv,mz=mz,
            conjunction=conjunction,first_visible_morning=dawn,
            days_after_conjunction=pmax.jd_ut1-conjunction.jd_ut1,
        )
    raise RuntimeError(
        f"No heliacal rising found for {star.name!r} within {max_days_after_conjunction} days "
        f"after conjunction at k_v={kv:.2f}."
    )
