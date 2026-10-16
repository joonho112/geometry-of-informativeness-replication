# Check the supplied inputs using the Python checksum utility.
arg <- grep('^--file=', commandArgs(FALSE), value = TRUE)
root <- dirname(dirname(normalizePath(sub('^--file=', '', arg[1]), mustWork = TRUE)))
python <- Sys.getenv('PYTHON', 'python3')
status <- system2(python, shQuote(file.path(root, 'manifest', 'verify.py')))
quit(status = status)
