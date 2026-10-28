#!/usr/bin/env python3
"""Astronomical geometry for the Chiang Saen heliacal-rising revision.

This module supplies the *geometry* layer used by the Chiang Saen
archaeoastronomy calculations.  It is designed to feed quantities such as
solar altitude, stellar altitude, stellar azimuth, and the Sun-star azimuth
separation into ``schaefer_visibility.py``.

Site
----
Pa Sak Gate (also called Chiang Saen Gate), WGS 84:

    latitude  = +20.2736 deg
    longitude = +100.0786 deg  (east-positive)

These are the coordinates previously used for the project simulations.

Model scope
-----------
The code is intentionally transparent and has no third-party dependencies.
It implements:

1. Julian-day conversion for Gregorian or Julian calendars.
2. Mean solar equatorial coordinates using the standard low-order solar
   formulae commonly tabulated by Meeus.
3. Linear proper-motion propagation in the J2000 inertial tangent plane.
4. Precession from J2000.0 to the epoch of date using the classical
   IAU-1976 / Meeus precession formulae.
5. Greenwich mean sidereal time and local hour angle.
6. Horizontal coordinates (altitude, azimuth measured clockwise from true
   north).
7. Optional Saemundsson atmospheric refraction for an apparent altitude.

The module does NOT implement heliacal visibility physics; that belongs in
``schaefer_visibility.py``.  It also does not attempt to reconstruct local
historical meteorology.

Accuracy note
-------------
The algorithms here are appropriate for the degree/sub-degree angular scales
of this archaeological alignment study, but they are not intended as a
replacement for a modern high-precision astrometric ephemeris.  The code
therefore keeps the geometry layer separate so that a future ephemeris
backend can be substituted without changing the visibility model.

Calendar note
-------------
For historical epochs, calendar labels can be ambiguous.  Functions below
accept either ``calendar='gregorian'`` or ``calendar='julian'`` explicitly.
Internally, all geometry is performed on Julian Day numbers.

Angle convention
----------------
Azimuth is measured clockwise from geographic (true) north:

    North =   0 deg
    East  =  90 deg
    South = 180 deg
    West  = 270 deg

All public angles are in degrees unless otherwise stated.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Optional, Tuple


# ---------------------------------------------------------------------------
# Project site
# ---------------------------------------------------------------------------

PA_SAK_LATITUDE_DEG = 20.2736
PA_SAK_LONGITUDE_DEG = 100.0786


@dataclass(frozen=True)
class Site:
    """Observing site.

    Longitude is east-positive.  Elevation is retained for metadata and
    possible later atmospheric refinements; it is not required for the
    spherical geometry implemented here.
    """

    name: str
    latitude_deg: float
    longitude_deg: float
    elevation_m: Optional[float] = None


PA_SAK_GATE = Site(
    name="Pa Sak Gate (Chiang Saen Gate)",
    latitude_deg=PA_SAK_LATITUDE_DEG,
    longitude_deg=PA_SAK_LONGITUDE_DEG,
    elevation_m=None,
)


@dataclass(frozen=True)
class StarJ2000:
    """J2000 astrometric input for one stellar candidate.

    Parameters
    ----------
    name : str
        Human-readable name.
    ra_deg : float
        J2000 right ascension in degrees.
    dec_deg : float
        J2000 declination in degrees.
    pm_ra_cosdec_mas_per_yr : float, default 0
        Proper motion in RA*cos(dec), milliarcseconds/year.
    pm_dec_mas_per_yr : float, default 0
        Proper motion in declination, milliarcseconds/year.

    Notes
    -----
    Radial velocity and perspective acceleration are ignored.  This is a
    deliberate simplification for the present naked-eye alignment problem.
    """

    name: str
    ra_deg: float
    dec_deg: float
    pm_ra_cosdec_mas_per_yr: float = 0.0
    pm_dec_mas_per_yr: float = 0.0


@dataclass(frozen=True)
class HorizontalPosition:
    """Horizontal coordinates for one object at one instant."""

    ra_of_date_deg: float
    dec_of_date_deg: float
    hour_angle_deg: float
    altitude_geometric_deg: float
    altitude_apparent_deg: float
    azimuth_deg: float
    zenith_distance_geometric_deg: float
    zenith_distance_apparent_deg: float


@dataclass(frozen=True)
class GeometrySnapshot:
    """Sun + star geometry at one instant for one observing site."""

    jd_ut1: float
    site: Site
    star_name: str

    sun_ra_deg: float
    sun_dec_deg: float
    sun_altitude_geometric_deg: float
    sun_altitude_apparent_deg: float
    sun_azimuth_deg: float

    star: HorizontalPosition
    sun_star_azimuth_difference_deg: float


# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------


def _wrap360(angle_deg: float) -> float:
    return angle_deg % 360.0


def _wrap180(angle_deg: float) -> float:
    """Wrap angle to [-180, +180)."""
    return (angle_deg + 180.0) % 360.0 - 180.0


def shortest_azimuth_separation_deg(az1_deg: float, az2_deg: float) -> float:
    """Absolute shortest azimuth separation in the range 0..180 deg."""
    return abs(_wrap180(az1_deg - az2_deg))


def ra_hms_to_deg(hours: float, minutes: float, seconds: float) -> float:
    """Convert right ascension h:m:s to degrees."""
    if hours < 0 or minutes < 0 or seconds < 0:
        raise ValueError("RA components must be non-negative")
    return 15.0 * (hours + minutes / 60.0 + seconds / 3600.0)


def dec_dms_to_deg(degrees: float, arcminutes: float, arcseconds: float) -> float:
    """Convert signed declination d:m:s to decimal degrees.

    The sign is taken from ``degrees``.  Example: (-23, 15, 30) -> -23.2583.
    """
    sign = -1.0 if degrees < 0 else 1.0
    return sign * (abs(degrees) + arcminutes / 60.0 + arcseconds / 3600.0)


# ---------------------------------------------------------------------------
# Julian day and calendar conversion
# ---------------------------------------------------------------------------


def julian_day(
    year: int,
    month: int,
    day: int,
    ut1_hours: float = 0.0,
    *,
    calendar: str = "gregorian",
) -> float:
    """Return Julian Day for a civil calendar date.

    Parameters
    ----------
    year : int
        Astronomical year numbering is supported (year 1 = AD 1, year 0 =
        1 BC) although the Chiang Saen project currently uses AD epochs.
    month, day : int
        Calendar month/day.
    ut1_hours : float
        UT1-like decimal hours.  For ancient epochs, exact civil clock time
        is uncertain because Delta-T is uncertain; a full-dawn scan is much
        less sensitive to that uncertainty than a single quoted clock time.
    calendar : {'gregorian', 'julian'}
        Calendar used only to label/convert the date.
    """
    if calendar not in {"gregorian", "julian"}:
        raise ValueError("calendar must be 'gregorian' or 'julian'")
    if not (1 <= month <= 12):
        raise ValueError("month must be 1..12")
    if not (1 <= day <= 31):
        raise ValueError("day must be 1..31")

    y = year
    m = month
    if m <= 2:
        y -= 1
        m += 12

    A = math.floor(y / 100.0)
    if calendar == "gregorian":
        B = 2 - A + math.floor(A / 4.0)
    else:
        B = 0

    jd0 = (
        math.floor(365.25 * (y + 4716))
        + math.floor(30.6001 * (m + 1))
        + day
        + B
        - 1524.5
    )
    return jd0 + ut1_hours / 24.0


# ---------------------------------------------------------------------------
# Solar coordinates
# ---------------------------------------------------------------------------


def mean_obliquity_deg(jd: float) -> float:
    """Mean obliquity of the ecliptic, degrees.

    Polynomial form widely used in Meeus-style low-order solar calculations.
    """
    T = (jd - 2451545.0) / 36525.0
    seconds = 21.448 - 46.8150 * T - 0.00059 * T * T + 0.001813 * T**3
    return 23.0 + 26.0 / 60.0 + seconds / 3600.0


def sun_equatorial_mean_deg(jd: float) -> Tuple[float, float]:
    """Return approximate geocentric solar RA, Dec of date [degrees].

    Uses a standard low-order solar longitude expansion adequate for the
    visibility geometry of this project.  Nutation and annual aberration are
    not required at the precision of the present archaeological comparison.
    """
    T = (jd - 2451545.0) / 36525.0

    L0 = _wrap360(280.46646 + 36000.76983 * T + 0.0003032 * T * T)
    M = _wrap360(
        357.52911
        + 35999.05029 * T
        - 0.0001537 * T * T
        + T**3 / 24490000.0
    )

    Mr = math.radians(M)
    C = (
        (1.914602 - 0.004817 * T - 0.000014 * T * T) * math.sin(Mr)
        + (0.019993 - 0.000101 * T) * math.sin(2.0 * Mr)
        + 0.000289 * math.sin(3.0 * Mr)
    )
    true_long = _wrap360(L0 + C)

    eps = math.radians(mean_obliquity_deg(jd))
    lam = math.radians(true_long)

    ra = math.degrees(math.atan2(math.cos(eps) * math.sin(lam), math.cos(lam)))
    dec = math.degrees(math.asin(math.sin(eps) * math.sin(lam)))
    return _wrap360(ra), dec


# ---------------------------------------------------------------------------
# Stellar proper motion and precession
# ---------------------------------------------------------------------------


def _normalize_vector(v: Tuple[float, float, float]) -> Tuple[float, float, float]:
    n = math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
    if n == 0.0:
        raise ValueError("Cannot normalize zero vector")
    return (v[0] / n, v[1] / n, v[2] / n)


def propagate_proper_motion_j2000(
    ra_deg: float,
    dec_deg: float,
    years_from_j2000: float,
    pm_ra_cosdec_mas_per_yr: float = 0.0,
    pm_dec_mas_per_yr: float = 0.0,
) -> Tuple[float, float]:
    """Propagate J2000 direction using tangent-plane proper motion.

    The method treats the two proper-motion components as a constant angular
    velocity vector in the J2000 tangent plane.  It ignores radial velocity
    and perspective acceleration.
    """
    ra = math.radians(ra_deg)
    dec = math.radians(dec_deg)

    r = (
        math.cos(dec) * math.cos(ra),
        math.cos(dec) * math.sin(ra),
        math.sin(dec),
    )
    e_ra = (-math.sin(ra), math.cos(ra), 0.0)
    e_dec = (
        -math.sin(dec) * math.cos(ra),
        -math.sin(dec) * math.sin(ra),
        math.cos(dec),
    )

    mas_to_rad = math.radians(1.0 / 3_600_000.0)
    mu_a = pm_ra_cosdec_mas_per_yr * mas_to_rad
    mu_d = pm_dec_mas_per_yr * mas_to_rad

    v = (
        mu_a * e_ra[0] + mu_d * e_dec[0],
        mu_a * e_ra[1] + mu_d * e_dec[1],
        mu_a * e_ra[2] + mu_d * e_dec[2],
    )

    rp = _normalize_vector(
        (
            r[0] + years_from_j2000 * v[0],
            r[1] + years_from_j2000 * v[1],
            r[2] + years_from_j2000 * v[2],
        )
    )

    ra2 = _wrap360(math.degrees(math.atan2(rp[1], rp[0])))
    dec2 = math.degrees(math.asin(max(-1.0, min(1.0, rp[2]))))
    return ra2, dec2


def precess_j2000_to_date_deg(
    ra_j2000_deg: float,
    dec_j2000_deg: float,
    jd: float,
) -> Tuple[float, float]:
    """Precess equatorial coordinates from J2000.0 to equinox of date.

    Classical IAU-1976 / Meeus precession formulae with the starting epoch
    fixed at J2000.0.  Adequate for the broad alignment scales used here.
    """
    t = (jd - 2451545.0) / 36525.0

    zeta_arcsec = 2306.2181 * t + 0.30188 * t * t + 0.017998 * t**3
    z_arcsec = 2306.2181 * t + 1.09468 * t * t + 0.018203 * t**3
    theta_arcsec = 2004.3109 * t - 0.42665 * t * t - 0.041833 * t**3

    zeta = math.radians(zeta_arcsec / 3600.0)
    z = math.radians(z_arcsec / 3600.0)
    theta = math.radians(theta_arcsec / 3600.0)

    ra = math.radians(ra_j2000_deg)
    dec = math.radians(dec_j2000_deg)

    A = math.cos(dec) * math.sin(ra + zeta)
    B = (
        math.cos(theta) * math.cos(dec) * math.cos(ra + zeta)
        - math.sin(theta) * math.sin(dec)
    )
    C = (
        math.sin(theta) * math.cos(dec) * math.cos(ra + zeta)
        + math.cos(theta) * math.sin(dec)
    )

    ra_date = math.degrees(math.atan2(A, B) + z)
    dec_date = math.degrees(math.asin(max(-1.0, min(1.0, C))))
    return _wrap360(ra_date), dec_date


def star_equatorial_of_date_deg(star: StarJ2000, jd: float) -> Tuple[float, float]:
    """Return stellar RA/Dec of date after proper motion + precession."""
    years = (jd - 2451545.0) / 365.25
    ra_pm, dec_pm = propagate_proper_motion_j2000(
        star.ra_deg,
        star.dec_deg,
        years,
        star.pm_ra_cosdec_mas_per_yr,
        star.pm_dec_mas_per_yr,
    )
    return precess_j2000_to_date_deg(ra_pm, dec_pm, jd)


# ---------------------------------------------------------------------------
# Sidereal time and horizontal coordinates
# ---------------------------------------------------------------------------


def greenwich_mean_sidereal_time_deg(jd_ut1: float) -> float:
    """Greenwich mean sidereal time, degrees."""
    T = (jd_ut1 - 2451545.0) / 36525.0
    theta = (
        280.46061837
        + 360.98564736629 * (jd_ut1 - 2451545.0)
        + 0.000387933 * T * T
        - T**3 / 38710000.0
    )
    return _wrap360(theta)


def local_mean_sidereal_time_deg(jd_ut1: float, longitude_deg: float) -> float:
    """Local mean sidereal time, east-positive longitude."""
    return _wrap360(greenwich_mean_sidereal_time_deg(jd_ut1) + longitude_deg)


def horizontal_from_equatorial_deg(
    ra_deg: float,
    dec_deg: float,
    jd_ut1: float,
    site: Site = PA_SAK_GATE,
    *,
    apply_refraction: bool = True,
    pressure_hpa: float = 1010.0,
    temperature_c: float = 20.0,
) -> HorizontalPosition:
    """Convert RA/Dec of date to altitude/azimuth at the site.

    Azimuth is clockwise from true north.  The geometric altitude is always
    returned.  ``altitude_apparent_deg`` includes a Saemundsson refraction
    correction when requested and when the source is close enough to the
    horizon for the formula to be applicable.
    """
    if not (-90.0 <= site.latitude_deg <= 90.0):
        raise ValueError("site latitude must lie in [-90, 90] deg")
    if pressure_hpa <= 0:
        raise ValueError("pressure_hpa must be positive")

    lst = local_mean_sidereal_time_deg(jd_ut1, site.longitude_deg)
    H_deg = _wrap180(lst - ra_deg)

    H = math.radians(H_deg)
    dec = math.radians(dec_deg)
    lat = math.radians(site.latitude_deg)

    sin_h = (
        math.sin(lat) * math.sin(dec)
        + math.cos(lat) * math.cos(dec) * math.cos(H)
    )
    sin_h = max(-1.0, min(1.0, sin_h))
    h_geo = math.degrees(math.asin(sin_h))

    # Azimuth, clockwise from north.
    # Adding 180 deg places H=-90, dec=0, lat=0 at due east (90 deg).
    az = math.degrees(
        math.atan2(
            math.sin(H),
            math.cos(H) * math.sin(lat) - math.tan(dec) * math.cos(lat),
        )
    )
    az = _wrap360(az + 180.0)

    if apply_refraction:
        h_app = h_geo + saemundsson_refraction_deg(
            h_geo,
            pressure_hpa=pressure_hpa,
            temperature_c=temperature_c,
        )
    else:
        h_app = h_geo

    return HorizontalPosition(
        ra_of_date_deg=_wrap360(ra_deg),
        dec_of_date_deg=dec_deg,
        hour_angle_deg=H_deg,
        altitude_geometric_deg=h_geo,
        altitude_apparent_deg=h_app,
        azimuth_deg=az,
        zenith_distance_geometric_deg=90.0 - h_geo,
        zenith_distance_apparent_deg=90.0 - h_app,
    )


# ---------------------------------------------------------------------------
# Atmospheric refraction
# ---------------------------------------------------------------------------


def saemundsson_refraction_deg(
    geometric_altitude_deg: float,
    *,
    pressure_hpa: float = 1010.0,
    temperature_c: float = 20.0,
) -> float:
    """Approximate atmospheric refraction in degrees.

    Saemundsson-style formula.  For altitudes below -1 deg the formula is not
    used and zero is returned, because it is not reliable far below the
    horizon.  The result is the amount added to geometric altitude.
    """
    h = geometric_altitude_deg
    if h < -1.0 or h >= 90.0:
        return 0.0

    # Formula is singular only far outside the intended range; h >= -1 keeps
    # us in the practical near-horizon domain.
    angle = math.radians(h + 10.3 / (h + 5.11))
    R_arcmin = 1.02 / math.tan(angle)
    R_arcmin *= (pressure_hpa / 1010.0) * (283.0 / (273.0 + temperature_c))
    return R_arcmin / 60.0


# ---------------------------------------------------------------------------
# Combined Sun-star snapshot
# ---------------------------------------------------------------------------


def geometry_snapshot(
    star: StarJ2000,
    jd_ut1: float,
    site: Site = PA_SAK_GATE,
    *,
    apply_refraction: bool = True,
    pressure_hpa: float = 1010.0,
    temperature_c: float = 20.0,
) -> GeometrySnapshot:
    """Return Sun and star geometry for one instant.

    This is the main interface intended for the heliacal-search layer.
    """
    star_ra, star_dec = star_equatorial_of_date_deg(star, jd_ut1)
    star_hor = horizontal_from_equatorial_deg(
        star_ra,
        star_dec,
        jd_ut1,
        site,
        apply_refraction=apply_refraction,
        pressure_hpa=pressure_hpa,
        temperature_c=temperature_c,
    )

    sun_ra, sun_dec = sun_equatorial_mean_deg(jd_ut1)
    sun_hor = horizontal_from_equatorial_deg(
        sun_ra,
        sun_dec,
        jd_ut1,
        site,
        apply_refraction=apply_refraction,
        pressure_hpa=pressure_hpa,
        temperature_c=temperature_c,
    )

    return GeometrySnapshot(
        jd_ut1=jd_ut1,
        site=site,
        star_name=star.name,
        sun_ra_deg=sun_ra,
        sun_dec_deg=sun_dec,
        sun_altitude_geometric_deg=sun_hor.altitude_geometric_deg,
        sun_altitude_apparent_deg=sun_hor.altitude_apparent_deg,
        sun_azimuth_deg=sun_hor.azimuth_deg,
        star=star_hor,
        sun_star_azimuth_difference_deg=shortest_azimuth_separation_deg(
            sun_hor.azimuth_deg,
            star_hor.azimuth_deg,
        ),
    )


# ---------------------------------------------------------------------------
# Convenience helpers for dawn scanning
# ---------------------------------------------------------------------------


def jd_for_local_mean_solar_time(
    year: int,
    month: int,
    day: int,
    local_mean_solar_hours: float,
    site: Site = PA_SAK_GATE,
    *,
    calendar: str = "gregorian",
) -> float:
    """Approximate JD corresponding to a local mean solar clock reading.

    Local mean solar time differs from UT1 by longitude/15 hours.  This
    helper is convenient for creating a broad dawn scan window.  It should
    not be interpreted as a reconstruction of historical civil time.
    """
    ut1_hours = local_mean_solar_hours - site.longitude_deg / 15.0
    # Let julian_day naturally handle a UT value outside 0..24 by adding the
    # fractional day; the JD remains continuous.
    return julian_day(year, month, day, ut1_hours, calendar=calendar)


# ---------------------------------------------------------------------------
# Internal checks
# ---------------------------------------------------------------------------


def _self_test() -> None:
    # J2000 reference instant.
    assert abs(julian_day(2000, 1, 1, 12.0, calendar="gregorian") - 2451545.0) < 1e-9

    # Site constants.
    assert abs(PA_SAK_GATE.latitude_deg - 20.2736) < 1e-12
    assert abs(PA_SAK_GATE.longitude_deg - 100.0786) < 1e-12

    # No proper motion and no elapsed time -> same direction.
    ra2, dec2 = propagate_proper_motion_j2000(56.0, 24.0, 0.0)
    assert abs(ra2 - 56.0) < 1e-12
    assert abs(dec2 - 24.0) < 1e-12

    # Precession at J2000 -> identity.
    ra3, dec3 = precess_j2000_to_date_deg(56.0, 24.0, 2451545.0)
    assert abs(ra3 - 56.0) < 1e-10
    assert abs(dec3 - 24.0) < 1e-10

    # Equatorial sanity check: at the equator, dec=0, H=-90 -> east horizon.
    test_site = Site("Equator", 0.0, 0.0)
    jd = 2451545.0
    lst = local_mean_sidereal_time_deg(jd, 0.0)
    ra = _wrap360(lst + 90.0)  # H = LST - RA = -90 deg
    p = horizontal_from_equatorial_deg(
        ra, 0.0, jd, test_site, apply_refraction=False
    )
    assert abs(p.altitude_geometric_deg) < 1e-9
    assert abs(p.azimuth_deg - 90.0) < 1e-9

    # Refraction near the horizon should be positive and of order 0.5 deg.
    r0 = saemundsson_refraction_deg(0.0)
    assert 0.4 < r0 < 0.7

    # Solar declination around the March 2000 equinox should be close to zero.
    jd_eq = julian_day(2000, 3, 20, 12.0, calendar="gregorian")
    _, sun_dec = sun_equatorial_mean_deg(jd_eq)
    assert abs(sun_dec) < 1.0


if __name__ == "__main__":
    _self_test()
    print("astronomy_geometry.py: internal checks passed")
    print(
        f"Pa Sak Gate: {PA_SAK_GATE.latitude_deg:.4f} deg N, "
        f"{PA_SAK_GATE.longitude_deg:.4f} deg E"
    )
