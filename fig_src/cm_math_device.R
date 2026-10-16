# Keep ordinary figure text in Source Sans 3 and mathematical labels in native CM.
# fontspec's no-math option is required to retain CM roman math/operators.
suppressPackageStartupMessages(library(tikzDevice))
CM_FONT_PATH <- paste0(normalizePath(file.path(ROOT, 'fig_src/fonts')), '/')
CM_FONT_OPTIONS <- paste0(
  '[Path={', CM_FONT_PATH,
  '},BoldFont=SourceSans3-Bold.otf]'
)
CM_PACKAGES <- c(
  '\\usepackage{tikz}',
  '\\usepackage[active,tightpage,xetex]{preview}',
  '\\usepackage[no-math]{fontspec}',
  '\\usepackage{xunicode}',
  '\\usepackage{amsmath,amssymb}',
  paste0('\\setmainfont{SourceSans3-Regular.otf}', CM_FONT_OPTIONS),
  paste0('\\setsansfont{SourceSans3-Regular.otf}', CM_FONT_OPTIONS),
  '\\PreviewEnvironment{pgfpicture}',
  '\\setlength\\PreviewBorder{0pt}'
)
FIGURE_PYTHON <- Sys.getenv('PYTHON', unset = 'python3')
if (!nzchar(Sys.which(FIGURE_PYTHON))) stop('Python 3 not found: set PYTHON or add python3 to PATH')
FIGURE_XELATEX <- Sys.which('xelatex')
if (!nzchar(FIGURE_XELATEX)) stop('XeLaTeX not found: add your TeX distribution binaries to PATH')
options(
  tikzDocumentDeclaration = '\\documentclass[12pt]{article}',
  tikzXelatexPackages = CM_PACKAGES,
  tikzXelatex = unname(FIGURE_XELATEX)
)
# Quartz scales R dash runs by the stroke width; tikzDevice writes them
# as fixed TeX-point runs. Restore the original dash/width ratio.
scale_dash_runs <- function(line) {
  if (!grepl('dash pattern=', line, fixed = TRUE)) return(line)
  width_match <- regmatches(line, regexec('line width= *([0-9.]+)pt', line))[[1]]
  stopifnot(length(width_match) == 2L)
  width <- as.numeric(width_match[2])
  positions <- gregexpr('(on|off) [0-9.]+pt', line, perl = TRUE)
  runs <- regmatches(line, positions)[[1]]
  replacements <- vapply(runs, function(run) {
    parts <- strsplit(run, ' ', fixed = TRUE)[[1]]
    sprintf('%s %.5fpt', parts[1], as.numeric(sub('pt$', '', parts[2])) * width)
  }, '')
  regmatches(line, positions) <- list(replacements)
  line
}
save_cm_figure <- function(p, name, h) {
  generated <- file.path(ROOT, 'fig_src/generated-cm')
  dir.create(generated, recursive = TRUE, showWarnings = FALSE)
  tex <- file.path(generated, paste0(name, '.tex'))
  tikzDevice::tikz(
    file = tex, width = 6.5, height = h, pointsize = PT,
    standAlone = TRUE, sanitize = FALSE, engine = 'xetex', packages = CM_PACKAGES,
    documentDeclaration = getOption('tikzDocumentDeclaration'),
    lwdUnit = 72.27 / 96, timestamp = FALSE, verbose = FALSE
  )
  tryCatch(if (is.function(p)) p() else print(p), finally = dev.off())
  lines <- readLines(tex, warn = FALSE)
  # Metrics run in a temporary directory; generated sources compile beside fonts.
  # Keep their font paths portable while using the same files for both stages.
  lines <- gsub(CM_FONT_PATH, '../fonts/', lines, fixed = TRUE)
  writeLines(vapply(lines, scale_dash_runs, ''), tex, useBytes = TRUE)
  status <- system2(FIGURE_PYTHON, c(
    shQuote(file.path(ROOT, 'fig_src/compile_cm_figure.py')),
    shQuote(tex)
  ))
  assert(status == 0L, paste('Native math figure build failed:', name))
}
