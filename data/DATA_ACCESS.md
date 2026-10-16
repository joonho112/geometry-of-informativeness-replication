# OCRS data and analysis

The repository includes aggregate OCRS inputs but no participant records or
identifiers. The default `./run_all.sh --reproduce` command uses these inputs
to recompute the posterior and sensitivity summaries.

## Run without participant data

From the repository directory:

```sh
./run_all.sh --posterior
./run_all.sh --sensitivity
```

The first command uses site estimates and sampling variances for six analysis
specifications. It calculates posterior summaries by adaptive quadrature and
writes them to `analysis/outputs/numerical-stability/`.

The second command uses saved site-level inputs from 500 bootstrap resamples
per outcome and all 256 one-site-per-paired-school selections. It recomputes
their posterior summaries and the variance-floor comparisons in
`analysis/outputs/ocrs-audit/`. It does not draw a new participant bootstrap.

These aggregate inputs cannot reconstruct participant eligibility or the
first-stage regressions. Use the participant file for those steps.

## Obtain and use the participant file

Download `OCRS_publicuse_impactanalysisdata.csv` from
[openICPSR project 202181](https://www.openicpsr.org/openicpsr/project/202181),
subject to the source's access and use terms. The scripts do not download it.
Set `OCRS_DATA` to the file location:

```sh
export OCRS_DATA=/path/to/OCRS_publicuse_impactanalysisdata.csv
./run_all.sh --empirical
```

Alternatively, place the file in this ignored `data/` directory. The script
checks that its SHA-256 checksum matches the file used for the paper:

```text
d20cbbaaf07ee8e7aad7e5d0e516b74f0a1114617fb3c7db9a4c7c2afcb681c7
```

The command reconstructs the first-stage estimates and HC2 standard errors,
repeats 500 participant bootstrap resamples per outcome and all 256 subset
selections, and recalculates the quadrature and variance-floor summaries.
Resampling occurs within site-by-arm-by-randomization-block cells. The original
sampler is also run to reproduce the random-number sequence; posterior
summaries use quadrature. Only aggregate results are written.

The analysis sample has 1,737 student-by-subject records from 1,683 students,
24 schools, and 32 school-by-subject sites (13 Algebra 1 and 19 English 9), with
144 within-site randomization blocks. The script checks these counts and
compares the first-stage estimates with the saved site inputs.

## Input files and interpretation

- `analysis/inputs/ocrs_fits.rds` contains site estimates, sampling variances,
  labels, and results from the original grid calculation.
- `analysis/inputs/ocrs-aggregate/` contains aggregate first-stage results and
  saved bootstrap and subset inputs.
- `analysis/inputs/tables/tbl_ocrs_sites.csv` contains school-by-subject estimates,
  variances, and counts. Its grid-based posterior columns are not used for
  the paper's posterior summaries.
- `analysis/inputs/figure-reference/` contains the plotting inputs;
  drawing the figures does not rerun the empirical analysis.

Posterior summaries use the normal working model with a reference
measure flat in the mean and standard deviation, conditional on supplied
sampling variances. Variation across bootstrap refits is separate from these
conditional posterior intervals. The shift in the credits bootstrap and
dependence due to shared students remain limitations; reproducing the numbers
does not resolve them. Further numerical details are in the
[running guide](../docs/running-the-analyses.md#interpretation).

Study reference: Rickles, J., Clements, M., Brodziak de los Reyes, I., Lachowicz,
M., Lin, S., & Heppen, J. (2024). A multisite randomized study of an online
learning approach to high school credit recovery. *Journal of Research on
Educational Effectiveness, 17*(3), 467–490.
