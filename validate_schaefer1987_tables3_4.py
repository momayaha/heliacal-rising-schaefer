#!/usr/bin/env python3
"""Validate the heliacal-search implementation against Schaefer (1987).

This script performs two independent published-result checks using the same
geometry, visibility, and continuous dawn-search code used for Chiang Saen.

1. Table 3: Sirius at latitude +30 deg, epoch AD 1950, m_z=6.0, with
   k_v = 0.20, 0.30, 0.40.  Published Delta is days after the vernal equinox.
2. Table 4: the same Sirius geometry and standard sky (k_v=0.20, m_z=6.0),
   while the adopted stellar magnitude is varied from -4 to +5.  Schaefer
   publishes Delta, stellar zenith distance Z*, and solar altitude h0.

The goal is not bit-for-bit reproduction of Schaefer's 1980s code.  Solar
ephemerides, precession, and numerical optimization are independently
implemented here.  We therefore use explicit validation tolerances:

    |Delta_calc - Delta_pub| <= 1.0 day
    |Z_calc - Z_pub|         <= 1.0 deg  (Table 4)
    |h0_calc - h0_pub|       <= 1.0 deg  (Table 4)

These tolerances are intentionally reported in the output rather than hidden.
"""
from __future__ import annotations

import csv
from pathlib import Path

import heliacal_search as hs
from validate_heliacal_search_schaefer1987 import (
    VALIDATION_SITE,
    VALIDATION_YEAR,
    VALIDATION_MZ,
    vernal_equinox_jd,
    make_schaefer_sirius_j2000,
)

DELTA_TOL_D = 1.0
ANGLE_TOL_DEG = 1.0

TABLE3 = [
    # k_v, published Delta [d]
    (0.20, 133.8),
    (0.30, 136.3),
    (0.40, 138.5),
]

TABLE4 = [
    # m_V, Delta [d], Z* [deg], h0 [deg]
    (-4.0, 130.6, 87.6, -2.6),
    (-3.0, 131.8, 87.5, -3.5),
    (-2.0, 133.1, 87.3, -4.5),
    (-1.0, 134.5, 87.0, -5.5),
    ( 0.0, 136.3, 86.4, -6.5),
    ( 1.0, 138.5, 85.6, -7.7),
    ( 2.0, 141.9, 83.8, -9.1),
    ( 3.0, 145.4, 86.7, -15.1),
    ( 4.0, 147.4, 85.2, -15.5),
    ( 5.0, 152.1, 80.6, -15.3),
]


def calculate(m_v: float, k_v: float, equinox_jd: float, sirius):
    r = hs.search_heliacal_rising_continuous(
        sirius,
        m_v,
        VALIDATION_YEAR,
        k_v,
        mz=VALIDATION_MZ,
        site=VALIDATION_SITE,
        calendar="gregorian",
        max_days_after_conjunction=100,
        require_calibration_domain=False,
        apply_refraction=False,
    )
    p = r.first_visible_morning.maximum
    if p is None:
        raise RuntimeError("No maximum-visibility point returned")
    return {
        "delta_calc_d": p.jd_ut1 - equinox_jd,
        "z_calc_deg": 90.0 - p.star_altitude_deg,
        "h0_calc_deg": p.sun_altitude_deg,
        "date": p.local_datetime.ymd,
        "time_lmst": p.local_datetime.hhmmss,
        "star_alt_deg": p.star_altitude_deg,
        "star_az_deg": p.star_azimuth_deg,
        "background_nl": p.background_nl,
        "margin_mag": p.visibility_margin_mag,
    }


def main():
    eq = vernal_equinox_jd(VALIDATION_YEAR)
    sirius = make_schaefer_sirius_j2000(eq)
    rows=[]

    for kv, dpub in TABLE3:
        c=calculate(-1.46,kv,eq,sirius)
        dr=c["delta_calc_d"]-dpub
        rows.append({
            "source_table":"Table 3", "m_v":-1.46, "k_v":kv,
            "delta_published_d":dpub, "delta_calculated_d":c["delta_calc_d"],
            "delta_residual_d":dr,
            "z_published_deg":"", "z_calculated_deg":c["z_calc_deg"], "z_residual_deg":"",
            "h0_published_deg":"", "h0_calculated_deg":c["h0_calc_deg"], "h0_residual_deg":"",
            "first_visible_date":c["date"], "max_visibility_time_lmst":c["time_lmst"],
            "star_alt_deg":c["star_alt_deg"], "star_az_deg":c["star_az_deg"],
            "background_nl":c["background_nl"], "visibility_margin_mag":c["margin_mag"],
            "delta_pass":abs(dr)<=DELTA_TOL_D, "z_pass":"", "h0_pass":"",
            "overall_pass":abs(dr)<=DELTA_TOL_D,
        })

    for mv,dpub,zpub,hpub in TABLE4:
        c=calculate(mv,0.20,eq,sirius)
        dr=c["delta_calc_d"]-dpub
        zr=c["z_calc_deg"]-zpub
        hr=c["h0_calc_deg"]-hpub
        ok=abs(dr)<=DELTA_TOL_D and abs(zr)<=ANGLE_TOL_DEG and abs(hr)<=ANGLE_TOL_DEG
        rows.append({
            "source_table":"Table 4", "m_v":mv, "k_v":0.20,
            "delta_published_d":dpub, "delta_calculated_d":c["delta_calc_d"],
            "delta_residual_d":dr,
            "z_published_deg":zpub, "z_calculated_deg":c["z_calc_deg"], "z_residual_deg":zr,
            "h0_published_deg":hpub, "h0_calculated_deg":c["h0_calc_deg"], "h0_residual_deg":hr,
            "first_visible_date":c["date"], "max_visibility_time_lmst":c["time_lmst"],
            "star_alt_deg":c["star_alt_deg"], "star_az_deg":c["star_az_deg"],
            "background_nl":c["background_nl"], "visibility_margin_mag":c["margin_mag"],
            "delta_pass":abs(dr)<=DELTA_TOL_D, "z_pass":abs(zr)<=ANGLE_TOL_DEG,
            "h0_pass":abs(hr)<=ANGLE_TOL_DEG, "overall_pass":ok,
        })

    out=Path(__file__).with_name("validation_schaefer1987_tables3_4.csv")
    with out.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys()));w.writeheader();w.writerows(rows)

    print("Schaefer (1987) Tables 3 & 4 validation")
    print(f"Tolerances: Delta ±{DELTA_TOL_D:.1f} d; Z*, h0 ±{ANGLE_TOL_DEG:.1f} deg")
    for r in rows:
        if r["source_table"]=="Table 3":
            print(f"Table 3 k={r['k_v']:.2f}: Delta pub={r['delta_published_d']:.1f}, calc={r['delta_calculated_d']:.3f}, residual={r['delta_residual_d']:+.3f} d -> {'PASS' if r['overall_pass'] else 'FAIL'}")
        else:
            print(f"Table 4 m={r['m_v']:+.0f}: dD={r['delta_residual_d']:+.3f} d, dZ={r['z_residual_deg']:+.3f} deg, dh0={r['h0_residual_deg']:+.3f} deg -> {'PASS' if r['overall_pass'] else 'FAIL'}")
    passed=sum(bool(r['overall_pass']) for r in rows)
    print(f"Overall: {passed}/{len(rows)} cases passed")
    print(f"CSV: {out}")
    if passed != len(rows):
        raise SystemExit(1)

if __name__=="__main__":
    main()
