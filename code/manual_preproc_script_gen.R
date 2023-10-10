#Script for phenotypes that keep hitting the stack limit on R and need the additional --max-ppsize=500000 command
library(tidyverse)

#Pull log files ----
preproc_logs_paths <- list.files("log/",
                                 "preprocess",
                                 full.names = TRUE,
                                 recursive = TRUE)

preproc_logs <- lapply(preproc_logs_paths,
                       read_lines)
names(preproc_logs) <- preproc_logs_paths


preproc_stack_fail_vec <- unlist(lapply(preproc_logs,
                                        function(x){any(x == "Error: protect(): protection stack overflow")}))

preproc_stack_fail <- preproc_logs[preproc_stack_fail_vec]

pheno <- str_split(names(preproc_stack_fail),
                   "\\/",
                   simplify = TRUE)[, 3, drop = TRUE]

group <- str_split(str_split(names(preproc_stack_fail),
                             "\\/",
                             simplify = TRUE)[, 4, drop = TRUE],
                   "\\.",
                   simplify = TRUE)[, 1, drop = TRUE]

geno <- str_split(str_split(names(preproc_stack_fail),
                            "\\/",
                            simplify = TRUE)[, 4, drop = TRUE],
                  "\\.",
                  simplify = TRUE)[, 2, drop = TRUE]

# rule preprocess_data:
#   input: code/preproc.R, data/mikropml/Pragmatic.full/complete.pan.csv
# output: data/mikropml/Pragmatic.full/complete.pan.dat_proc.Rds
# log: log/Pragmatic.full/complete.pan.preprocess_data.txt
# jobid: 118
# benchmark: benchmarks/Pragmatic.full/complete.pan.preprocess_data.txt
# reason: Missing output files: data/mikropml/Pragmatic.full/complete.pan.dat_proc.Rds
# wildcards: phenotype=Pragmatic.full, group=complete, genome=pan
# resources: mem_mb=8GB, disk_mb=1000, disk_mib=954, tmpdir=<TBD>, ncores=20


for(i in 1:length(preproc_stack_fail)){
#generate specific pre-processing script

script <- c("library(mikropml)",
            paste0("data_raw <- readr::read_csv(\"data/mikropml/", pheno[i], "/", group[i], ".", geno[i], ".csv\")"),
            "dim(data_raw)",
            paste0("data_processed <- preprocess_data(data_raw, outcome_colname = \"", pheno[i], "\", group_neg_corr = TRUE, remove_var = 'nzv')"),
            "summary(data_processed)",
            paste0("saveRDS(data_processed, file = \"data/mikropml/", pheno[i], "/", group[i], ".", geno[i], ".dat_proc.Rds\")"))

write_lines(script,
            file = paste0("code/expand_stack_preproc.", pheno[i], ".", group[i], ".", geno[i], ".R"),
            sep = "\n")

#generate specific pre-processing sbat

sh <- c("#!/bin/sh",
        "#SBATCH --job-name=expand_stack_preproc",
        "#SBATCH --mail-user=emilycma@umich.edu",
        "#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE",
        paste0("#SBATCH --output=log/", pheno[i], "/", group[i], ".", geno[i], ".expand_stack_preproc.txt"),
        "#SBATCH --export=ALL",
        "#SBATCH --partition=standard",
        "#SBATCH --account=esnitkin1",
        "#SBATCH --nodes=1 --ntasks=1 --mem=10g --time=24:00:00",
        "",
        paste0("Rscript --max-ppsize=500000 code/expand_stack_preproc.", pheno[i], ".", group[i], ".", geno[i], ".R"))

write_lines(sh,
            file = paste0("code/expand_stack_preproc.", pheno[i], ".", group[i], ".", geno[i], ".sbat"),
            sep = "\n")

}
