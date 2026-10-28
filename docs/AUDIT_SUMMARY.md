# Technical audit summary — 2026-10-26

Input: `data+code+results.zip` supplied by manuscript authors.

- Seven distinct Python scripts; no third-party package imports found.
- `stars.csv` has 28 records; 27 marked `primary_analysis=1`, Maia retained as extra alternative `primary_analysis=0`.
- Full archival run is 270 rows, two k_V values x five epochs x 27 stars; sensitivity pairs make 135 rows.
- Internal test entrypoints (geometry, visibility, heliacal search): PASS.
- `validate_schaefer1987.py`: PASS 3/3 Table 3 Sirius extinction cases (±0.50 days).
- `validate_heliacal_search_schaefer1987.py`: PASS 1/1 end-to-end case (±0.50 days).
- `validate_schaefer1987_tables3_4.py`: PASS 13/13 published-case checks (±1.0 day and ±1.0 deg as applicable).
- `run_chiangsaen.py --job-stop 1`: PASS, produced 2 rows; agrees with first two rows of archival `reproduced_results.csv` after numeric-format normalization.
- `chiangsaen_heliacal_results_FINAL.csv` and `reproduced_results.csv` share the same scientific values in all 270 rows; formatting only differs.
- `--step-seconds` argument appears unused in the supplied batch driver; do not assume it is active.
- Full batch of all 270 results NOT re-executed during audit.
- No statement is made here that every source-model equation was independently checked.
- Python scripts were copied byte-for-byte without modification in this draft package.
