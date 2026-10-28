Colab-friendly flat bundle for the Chiang Saen heliacal-visibility R2 code.
All .py files and stars.csv are in the same directory, so imports work directly.
Suggested checks:
  python schaefer_visibility.py
  python astronomy_geometry.py
  python heliacal_search.py
  python heliacal_search.py --demo-alcyone-ad500
  python validate_schaefer1987.py
  python validate_heliacal_search_schaefer1987.py
  python validate_schaefer1987_tables3_4.py
  python run_chiangsaen.py --stars stars.csv --output reproduced_results.csv --job-stop 1
Full run:
  python run_chiangsaen.py --stars stars.csv --output reproduced_results.csv --workers 1
