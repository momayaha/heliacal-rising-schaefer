# Checklist before publishing publicly

1. **Author confirmed.** `CITATION.cff` and `LICENSE` list Farung Surina as the sole software creator and copyright holder. Before publishing, verify ownership/permission for input datasets and any third-party material.
2. **Resolve reproducibility gap.** Decide whether to provide scripts/notebook(s) for post-processing the secondary CSV tables or add a precise, verified methodology for their derivation. See README.
3. **CLI check.** `--step-seconds` is accepted by `run_chiangsaen.py` but does not control the algorithm in the currently supplied code. Either wire it correctly and rerun tests or remove the misleading option; do not silently change scientific outputs.
4. **Scientific fidelity.** Confirm the published results and assumptions match the latest manuscript, especially how the local horizon is handled, what definition of heliacal event and azimuth is used, and any line corrections requested in minor revision.
5. **Re-run tests** in fresh environment. A short sample and three validation scripts passed Python 3.13.5 in this audit; a full 270-row rerun was not performed.
6. **Credit datasets.** Verify rights/attribution for astrometry, magnitudes, and archaeological alignment inputs, and add citations in README.
7. **AI transparency.** README identifies ChatGPT coding assistance and Google Colab computing environment; verify final phrasing reflects actual use and add corresponding manuscript disclosure.
8. **Start with GitHub Private repository.** Upload the *contents* of this extracted folder (not the ZIP itself), check all files and secrets; after items 1–7, make it Public.
9. **Connect Zenodo GitHub integration before first stable release.** After validation, create release v1.0.0 and verify the Zenodo version DOI/metadata.
10. **Submit to ASCL**, using code URL and reference; registration is an editorial process, not instant.
11. **Update manuscript** Code Availability, software citations, reviewer response with actual persistent identifiers.

This package is safe to stage in a **Private** repository; do not make it public or create a stable release until pre-publication checks are completed. The source calculation scripts are unchanged.
