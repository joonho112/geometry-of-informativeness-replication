# Saved inputs for Figures 1–8

These files contain the results used to draw the figures. Figure
creation does not refit models or access participant-level OCRS records.

| File | Contents |
|---|---|
| `tbl_ratios_shape.csv` | Paired simulation loss ratios |
| `tbl_fair.csv` | The 27 held-out-cell metamodel comparisons, using true generating reliability profiles |
| `design-illustration-data.csv` | Hypothetical lognormal and local design curves |
| `ocrs-corrected-summaries.csv` | Empirical summaries calculated by adaptive quadrature |
| `ocrs_fits_corrected.rds` | Site estimates, variances, posterior summaries, and draws used in Figures 3 and 8 |

The CSVs match those used for the paper. The RDS retains the values used for
plotting and omits unused posterior matrices and historical sensitivity
results. Input checksums are in `manifest/figure-source-lock.json`.
Recomputed analyses write to `analysis/outputs/` without replacing these files.
