# Analysis uses the original Xena-derived survival file supplied by the author.
revision_dir <- "."
input_dir <- file.path(revision_dir,"P3H2_prognosis_inputs")
output_dir <- file.path(revision_dir,"P3H2_prognosis_outputs")
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
Sys.setenv(P3H2_PROG_INPUTS=input_dir,P3H2_PROG_OUTPUTS=output_dir)
for(p in c("survival","survminer"))if(!requireNamespace(p,quietly=TRUE))stop("Missing R package: ",p)
source(file.path(revision_dir,"P3H2_survival_immune_analysis.R"),encoding="UTF-8")
message("Completed using original Xena endpoints: 85 patients, 27 deaths. Results: ",output_dir)
