# Explicit dependency installation; never invoked by the analysis driver.
# Rscript setup.R --check
# Rscript setup.R --install
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !args %in% c('--check', '--install')) stop('Use --check or --install')
needed <- c('jsonlite', 'digest', 'ggplot2', 'patchwork', 'tikzDevice')
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (args == '--install' && length(missing)) install.packages(missing, repos = 'https://cloud.r-project.org')
versions <- vapply(needed, function(p) if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else 'MISSING', character(1))
print(versions)
cat('R:', as.character(getRversion()), '\n')
cat('Also required: Python3; figures: TeX Live/MacTeX (latexmk,xelatex,TikZ,fontspec); figures also need Poppler; optional PDF QA: pypdf,pdfplumber.\n')
if (any(versions == 'MISSING')) quit(status = 1)
