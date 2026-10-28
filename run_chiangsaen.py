#!/usr/bin/env python3
"""Run the Chiang Saen Schaefer-style heliacal-rising reconstruction.

Primary analysis: 27 Nakshatra candidates x 5 epochs x k_v={0.20,0.30}.
The annual conjunction is computed once per star/epoch and shared by the two
extinction scenarios.  Each dawn is scanned at the requested cadence.

Historical date convention: Julian calendar before the Gregorian reform;
Gregorian calendar for AD 2000.  Times are local mean solar time (LMST).
Geometric altitudes are used in the Schaefer-1987 core calculation.

The CSV deliberately retains first, maximum-margin, and last visible points
on the first-visible morning.  The formal heliacal event is the first visible
*date*; retaining all three within-morning diagnostics avoids silently
imposing a unique time definition not required by Schaefer's date criterion.
"""
from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import astronomy_geometry as ag
import heliacal_search as hs
import schaefer_visibility as sv

EPOCHS=(1,500,1000,1500,2000)
KV_VALUES=(0.20,0.30)


def calendar_for_epoch(year:int)->str:
    return 'julian' if year < 1582 else 'gregorian'


def load_primary_stars(path:Path):
    rows=[]
    with path.open(newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f):
            if int(r['primary_analysis']) != 1:
                continue
            star=ag.StarJ2000(
                name=f"{r['nakshatra']} | {r['modern_identification']}",
                ra_deg=float(r['ra_j2000_deg']),
                dec_deg=float(r['dec_j2000_deg']),
                pm_ra_cosdec_mas_per_yr=float(r['pm_ra_cosdec_mas_yr']),
                pm_dec_mas_per_yr=float(r['pm_dec_mas_yr']),
            )
            rows.append((r,star,float(r['v_mag_adopted'])))
    if len(rows)!=27:
        raise RuntimeError(f'Expected 27 primary candidates, found {len(rows)}')
    return rows


def _point_fields(prefix,p):
    if p is None:
        return {f'{prefix}_{k}':'' for k in (
            'date','time_lmst','sun_alt_deg','sun_az_deg','star_alt_deg','star_az_deg',
            'sun_star_azdiff_deg','visibility_margin_mag','limiting_mag_v','airmass',
            'extinction_mag','background_nl','inside_1987_fit_domain')}
    return {
        f'{prefix}_date':p.local_datetime.ymd,
        f'{prefix}_time_lmst':p.local_datetime.hhmmss,
        f'{prefix}_sun_alt_deg':f'{p.sun_altitude_deg:.6f}',
        f'{prefix}_sun_az_deg':f'{p.sun_azimuth_deg:.6f}',
        f'{prefix}_star_alt_deg':f'{p.star_altitude_deg:.6f}',
        f'{prefix}_star_az_deg':f'{p.star_azimuth_deg:.6f}',
        f'{prefix}_sun_star_azdiff_deg':f'{p.sun_star_azimuth_difference_deg:.6f}',
        f'{prefix}_visibility_margin_mag':f'{p.visibility_margin_mag:.6f}',
        f'{prefix}_limiting_mag_v':f'{p.limiting_magnitude_v:.6f}',
        f'{prefix}_airmass':f'{p.airmass:.6f}',
        f'{prefix}_extinction_mag':f'{p.extinction_mag:.6f}',
        f'{prefix}_background_nl':f'{p.background_nl:.6f}',
        f'{prefix}_inside_1987_fit_domain':int(p.within_calibration_domain),
    }


def _run_star_epoch(job):
    meta,star_dict,vmag,epoch,step_seconds,strict_domain=job
    star=ag.StarJ2000(**star_dict)
    cal=calendar_for_epoch(epoch)
    conjunction=hs.find_annual_conjunction(star,epoch,ag.PA_SAK_GATE,calendar=cal)
    outs=[]
    for kv in KV_VALUES:
        result=hs.search_heliacal_rising_continuous(
            star, vmag, epoch, kv, mz=sv.DEFAULT_MZ, site=ag.PA_SAK_GATE,
            calendar=cal, coarse_step_seconds=120.0, max_days_after_conjunction=120,
            require_calibration_domain=strict_domain, apply_refraction=False,
            conjunction=conjunction,
        )
        dawn=result.first_visible_morning
        days_after=result.days_after_conjunction
        out={
            'no':meta['no'],'nakshatra':meta['nakshatra'],
            'lanna_transliteration':meta['lanna_transliteration'],
            'modern_identification':meta['modern_identification'],
            'v_mag':f'{vmag:.3f}','epoch_ad':epoch,'calendar':cal,
            'kv_mag_per_airmass':f'{kv:.2f}','mz':f'{sv.DEFAULT_MZ:.2f}',
            'diagnostic_threshold_step_seconds':'5.0',
            'continuous_maximum_optimization':'yes',
            'conjunction_date':conjunction.local_datetime.ymd,
            'conjunction_time_lmst':conjunction.local_datetime.hhmmss,
            'conjunction_separation_deg':f'{conjunction.angular_separation_deg:.6f}',
            'first_visible_morning':f'{dawn.year:04d}-{dawn.month:02d}-{dawn.day:02d}',
            'days_after_conjunction_at_max':f'{days_after:.6f}',
            'evaluated_points_on_first_visible_morning':dawn.evaluated_points,
        }
        out.update(_point_fields('first',dawn.first_visible))
        out.update(_point_fields('max',dawn.maximum))
        out.update(_point_fields('last',dawn.last_visible))
        outs.append(out)
    return outs


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--stars',default='stars.csv')
    ap.add_argument('--output',default='chiangsaen_heliacal_results.csv')
    ap.add_argument('--step-seconds',type=float,default=30.0)
    ap.add_argument('--workers',type=int,default=1)
    ap.add_argument('--strict-domain',action='store_true')
    ap.add_argument('--job-start',type=int,default=0,help='0-based start index in the 135 star-epoch job list')
    ap.add_argument('--job-stop',type=int,default=None,help='exclusive stop index in the 135 star-epoch job list')
    args=ap.parse_args()

    stars=load_primary_stars(Path(args.stars))
    jobs=[]
    for meta,star,vmag in stars:
        sd={'name':star.name,'ra_deg':star.ra_deg,'dec_deg':star.dec_deg,
            'pm_ra_cosdec_mas_per_yr':star.pm_ra_cosdec_mas_per_yr,
            'pm_dec_mas_per_yr':star.pm_dec_mas_per_yr}
        ms={k:meta[k] for k in ('no','nakshatra','lanna_transliteration','modern_identification')}
        for epoch in EPOCHS:
            jobs.append((ms,sd,vmag,epoch,args.step_seconds,args.strict_domain))

    all_jobs=jobs
    stop=len(all_jobs) if args.job_stop is None else min(args.job_stop,len(all_jobs))
    jobs=all_jobs[max(0,args.job_start):stop]
    print(f'Running job slice {max(0,args.job_start)}:{stop} = {len(jobs)} star-epoch jobs -> {len(jobs)*2} result rows')
    results=[]
    if args.workers<=1:
        for i,j in enumerate(jobs,1):
            results.extend(_run_star_epoch(j))
            if i%10==0 or i==len(jobs): print(f'  completed {i}/{len(jobs)} star-epoch jobs',flush=True)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futs=[ex.submit(_run_star_epoch,j) for j in jobs]
            for i,fut in enumerate(as_completed(futs),1):
                results.extend(fut.result())
                if i%10==0 or i==len(jobs): print(f'  completed {i}/{len(jobs)} star-epoch jobs',flush=True)

    results.sort(key=lambda r:(int(r['no']),int(r['epoch_ad']),float(r['kv_mag_per_airmass'])))
    out=Path(args.output)
    with out.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(results[0].keys())); w.writeheader(); w.writerows(results)
    print(f'Wrote {out} ({len(results)} rows)')

if __name__=='__main__': main()
