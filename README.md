# Heliacal-rising simulations: Schaefer (1987) implementation


> **Public research software repository.**
> This repository provides seven modular Python scripts,
> input data, validation procedures, and numerical results
> supporting the reproducible heliacal-rising simulations
> presented in the Chiang Saen archaeoastronomy study.


Python implementation of historical heliacal-rising computations based on the visibility model of Schaefer (1987), demonstrated with 27 Nakshatra candidates and archaeological sight lines at Chiang Saen, Thailand.

## Research context

Software associated with the manuscript:

*A Workflow for Evaluating Stellar-Alignment Candidates Using Multi-Epoch Heliacal-Rising Simulations: A Case Study of Chiang Saen, Thailand*, submitted to *Astronomy and Computing*, manuscript ASCOM-D-26-00288.

This archive contains seven Python scripts, input astrometry, validation outputs, and example derived results. This is a reproducibility-oriented research implementation, not a universal historical sky or full-fidelity ephemeris package.

## How the software works

The model consists of: (1) astronomical geometry, including precession and proper-motion propagation; (2) stellar detectability during dawn using a Schaefer-style visibility model; (3) search for the first visible morning after conjunction; (4) repetition over 27 primary stars, five epochs (AD 1, 500, 1000, 1500, 2000), and extinction coefficients k_V = 0.20 and 0.30 mag/airmass; and (5) validation against published Schaefer (1987) cases.

The code uses local mean solar time (LMST) labels and explicitly chosen historical Julian/Gregorian calendars. Heliacal visibility is **not** determined by imposing a universal fixed Sun altitude. The calculation uses a modeled visibility threshold, and its empirical calibration-domain flags are included in output fields.

The adopted observing location is Pa Sak / Chiang Saen Gate (20.2736 degrees N, 100.0786 degrees E). These simulations do not incorporate a measured local landscape-horizon profile or historically reconstructed obstacles.

## Installation

Tested in this software-package audit with Python 3.13.5. All seven scripts use the Python standard library and do not require pip-installed third-party packages.


1. Download this repository using **Code > Download ZIP**,
   or clone it using:

   git clone https://github.com/momayaha/heliacal-rising-schaefer.git

2. Open Terminal / Command Prompt in the repository folder.
3. Run the validation scripts or the Chiang Saen simulation
   commands described below.

No special installation is needed for the code as supplied here. See `requirements.txt`.

## Quick self-tests

```bash
python astronomy_geometry.py
python schaefer_visibility.py
python heliacal_search.py
```

## Validation against Schaefer (1987)

```bash
python validate_schaefer1987.py
python validate_heliacal_search_schaefer1987.py
python validate_schaefer1987_tables3_4.py
```

Three scripts compare independent numerical calculations with the published Sirius example / Tables 3 and 4. The existing CSV outputs are archived in `validation/`. Each validation script also writes a CSV into the script directory when executed. This independent package audit on Python 3.13.5 obtained PASS in all published cases tested (the 13-case Tables 3/4 comparison yielded 13/13 PASS), but additional independent scientific verification remains valuable.

## Reproduce a short sample

```bash
python run_chiangsaen.py --stars stars.csv --output trial_one_job.csv --job-stop 1
```

This runs one star-epoch job and writes two rows (one each for k_V = 0.20 and 0.30). The full calculation evaluates 27 x 5 x 2 = 270 results:

```bash
python run_chiangsaen.py --stars stars.csv --output trial_full_results.csv --workers 1
```

This full run can be computationally intensive. Use `--job-start` and `--job-stop` to run subsets. The archival output is provided as `results/reproduced_results.csv`. A sample job was executed during this packaging audit, but a full 270-row rerun was **not** completed as part of the audit.

**Important:** The command-line `--step-seconds` option appears in `run_chiangsaen.py`, but its parsed value is not passed through to the actual search call in the supplied code. Do not represent it as an active configurable scientific setting without fixing and retesting it. The science run currently passes `coarse_step_seconds=120.0`, and the output records additional 5-second diagnostics; see source for the precise internal search/optimization procedure.

## Files

| Item | Contents |
|---|---|
| `astronomy_geometry.py` | Historical dates, solar/stellar geometry, proper motion, precession |
| `schaefer_visibility.py` | Visibility physics, atmospheric effects and detection threshold |
| `heliacal_search.py` | Numerical search for heliacal rising |
| `run_chiangsaen.py` | Batch driver for the primary 27 x 5 x 2 calculation |
| `validate_schaefer1987.py` | Reference Sirius validation |
| `validate_heliacal_search_schaefer1987.py` | End-to-end heliacal-search validation |
| `validate_schaefer1987_tables3_4.py` | Published Table 3 / Table 4 validation |
| `stars.csv` | 28 candidate data records: 27 primary rows and one alternative Krittika (Maia) row |
| `results/reproduced_results.csv` | 270-row main simulation table, 56 columns |
| `results/kv_sensitivity_reproduced.csv` | 135-row paired k_V sensitivity table, 112 columns |
| `results/archaeological_matches_reproduced.csv` | Archived architectural alignment matches |
| `results/table3_replacement_k020.csv` | Derived manuscript Table 3 at k_V = 0.20 |
| `validation/*.csv` | Archived outputs of the three verification scripts |

The primary analysis includes Alcyone as the Krittika representative; Maia is retained as a secondary catalogue alternative but is not run by `run_chiangsaen.py` when `primary_analysis=0`.

The original archive also contains a second 270-row file, `chiangsaen_heliacal_results_FINAL.csv`. Its entries were compared field-by-field against `reproduced_results.csv`: values agree, and differences are number formatting only. To avoid ambiguity, this draft repository retains only `reproduced_results.csv`.

## What is reproducible directly, and what still needs documenting

The primary 270-row run and validation CSVs have associated executable scripts. The provided source bundle does **not** include an independent driver that regenerates all three *derived* exports (`kv_sensitivity_reproduced.csv`, `archaeological_matches_reproduced.csv`, and `table3_replacement_k020.csv`) from the main CSV. Those outputs are included for inspection and audit; post-processing procedures or scripts should be documented before claiming that every secondary table can be recreated with one command.

## Limitations

- Model visibility depends on extinction, assumed dark-sky limiting magnitude, and the empirical domain of Schaefer's twilight model.
- The astrometric / solar approximations prioritize interpretability, not precision-ephemeris parity.
- A local Eastern-horizon terrain profile and any historical vegetation/building obstructions are not encoded in this version.
- The program evaluates astronomical plausibility; it does not establish cultural intention or confirm any archaeological alignment.

## References

- Schaefer, B. E. (1986), *Atmospheric Extinction Effects on Stellar Alignments*.
- Schaefer, B. E. (1987), *Heliacal Rise Phenomena*.
- Schaefer, B. E. (1993), *Astronomy and the Limits of Vision*.
- Associated Chiang Saen manuscript, *Astronomy and Computing*, ASCOM-D-26-00288 (under revision as of October 2026).


## Software authorship, citation and license

Software author: Farung Surina
Institution: Chiang Rai Rajabhat University
License: MIT

GitHub:
https://github.com/momayaha/heliacal-rising-schaefer

Zenodo:
https://doi.org/10.5281/zenodo.23724651

Please cite the archived software release when using
this implementation in published research.
Citation metadata are provided in `CITATION.cff`.

## Software development and AI-use statement

Farung Surina developed this research software with assistance from OpenAI's ChatGPT for code generation, debugging, and refinement. Google Colaboratory was used as an environment for running and testing the scripts. Farung Surina retains responsibility for scientific decisions and the validation and interpretation of results. This acknowledgement does not list ChatGPT or Google Colaboratory as software authors.
