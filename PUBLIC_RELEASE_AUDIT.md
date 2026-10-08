# Public release audit

1. **Public scripts:** 24 (`.R`, `.py`, `.Rmd` and `.sh` under `code/`).
2. **Derived data files:** 46, including 13 Supplementary Table CSV files.
3. **Figures 1–7:** present in both PDF and SVG (14 final figure files).
4. **Supplementary Tables S1–S11:** present; split tables S4a/S4b and S5a/S5b give 13 CSV files.
5. **Reproducibility matrix:** 20 READY; 0 PARTIAL.
6. **Portable-script verification:** PASSED. Organoid-versus-2D and respiratory-submodule outputs were exact; CPTAC outputs were exact or numerically equivalent within absolute `1e-10` and relative `1e-8` tolerance.
7. **Remaining absolute private paths:** none detected in public code or documentation.
8. **External inputs users must download:** large organoid/Figshare source matrices, DepMap 24Q2 dependency matrices, MSigDB 2025.1.Hs source gene sets, HGNC source table, TCGA COAD/READ Xena files, GSE132465 processed data and CPTAC-2 Colon source matrices. The exact DepMap 24Q2 `Model.csv` is retained as small version-specific metadata.
9. **Licence decision required:** YES. The author must select a code licence; no third-party data licence is assigned.
10. **Repository size:** 38,043,846 bytes (36.28 MiB).

Security scan result: no user-specific Windows/macOS paths, local session paths, API keys, tokens, passwords, private email addresses or direct patient identifiers were detected. Pseudonymous public sample/model identifiers required for reproducibility remain in derived tables.

No Git repository was initialized and nothing was uploaded.

