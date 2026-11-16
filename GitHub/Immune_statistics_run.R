# Run using Rscript; input paths are relative to this script.
arg <- grep('^--file=', commandArgs(FALSE), value=TRUE)
if (!length(arg)) stop('Run with Rscript Immune_statistics_run.R')
base <- dirname(normalizePath(sub('^--file=', '', arg[1]), winslash='/'))
Sys.setenv(IMMUNE_INPUT_DIR=file.path(base,'Immune_inputs'),
           IMMUNE_OUTPUT_DIR=file.path(base,'Immune_outputs'))
source(file.path(base,'Immune_panel_statistics.R'), encoding='UTF-8')
