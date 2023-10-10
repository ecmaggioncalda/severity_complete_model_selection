#LIBRARIES ----
library(tidyverse)
library(mikropml)

#METADATA ----
metadata <- read_csv("../../data/patient_metadata.csv")

metadata_sub <- metadata %>%
  select(genome_id,
         clade) %>%
  mutate(across(everything(), ~as.character(.x)))

metadata_sub <- as.data.frame(metadata_sub)
rownames(metadata_sub) <- metadata_sub$genome_id
metadata_sub <- metadata_sub[, -1, drop = FALSE]

#GENOS ----
core_df <- read_delim("../../data/core_mat_sift.tsv")

core_df_sub <- core_df %>%
  mutate(variant = gsub("_$", "", variant)) %>%
  column_to_rownames("variant")

core_df_merge <- t(core_df_sub)

if(sum(rownames(metadata_sub) %in% rownames(core_df_merge)) != length(rownames(metadata_sub))){
  stop("mismatch between pheno and geno contents")
}

index <- match(rownames(metadata_sub), rownames(core_df_merge))

core_df_ordered <- core_df_merge[index, , drop = FALSE]

if(sum(rownames(metadata_sub) == rownames(core_df_ordered)) != length(rownames(metadata_sub))){
  stop("mismatch between pheno and geno contents")
}

core_df_frame <- cbind(metadata_sub,
                       core_df_ordered)

core_preproc <- preprocess_data(dataset = core_df_frame,
                                outcome_colname = "clade",
                                method = NULL,
                                collapse_corr_feats = TRUE,
                                group_neg_corr = TRUE,
                                to_numeric = FALSE,
                                remove_var = "zv")

core_preproc_df <- core_preproc$dat_transformed
colnames(core_preproc_df)[grep("grp", colnames(core_preproc_df))] <- gsub("$", "_core", colnames(core_preproc_df)[grep("grp", colnames(core_preproc_df))])

gene_df <- read_delim("../../data/gene_mat_sift.tsv")

gene_df_sub <- gene_df %>%
  mutate(variant = gsub("_$", "", variant)) %>%
  column_to_rownames("variant")

gene_df_merge <- t(gene_df_sub)

if(sum(rownames(metadata_sub) %in% rownames(gene_df_merge)) != length(rownames(metadata_sub))){
  stop("mismatch between pheno and geno contents")
}

index <- match(rownames(metadata_sub), rownames(gene_df_merge))

gene_df_ordered <- gene_df_merge[index, , drop = FALSE]

if(sum(rownames(metadata_sub) == rownames(gene_df_ordered)) != length(rownames(metadata_sub))){
  stop("mismatch between pheno and geno contents")
}

gene_df_frame <- cbind(metadata_sub,
                       gene_df_ordered)

gene_preproc <- preprocess_data(dataset = gene_df_frame,
                                outcome_colname = "clade",
                                method = NULL,
                                collapse_corr_feats = TRUE,
                                group_neg_corr = TRUE,
                                to_numeric = FALSE,
                                remove_var = "zv")

gene_preproc_df <- gene_preproc$dat_transformed
colnames(gene_preproc_df)[grep("grp", colnames(gene_preproc_df))] <- gsub("$", "_gene", colnames(gene_preproc_df)[grep("grp", colnames(gene_preproc_df))])

pan_df <- read_delim("../../data/pan_mat.tsv")

pan_df_sub <- pan_df %>%
  mutate(variant = gsub("_$", "", variant)) %>%
  column_to_rownames("variant")

pan_df_merge <- t(pan_df_sub)

if(sum(rownames(metadata_sub) %in% rownames(pan_df_merge)) != length(rownames(metadata_sub))){
  stop("mismatch between pheno and geno contents")
}

index <- match(rownames(metadata_sub), rownames(pan_df_merge))

pan_df_ordered <- pan_df_merge[index, , drop = FALSE]

if(sum(rownames(metadata_sub) == rownames(pan_df_ordered)) != length(rownames(metadata_sub))){
  stop("mismatch between pheno and geno contents")
}

pan_df_frame <- cbind(metadata_sub,
                      pan_df_ordered)

pan_preproc <- preprocess_data(dataset = pan_df_frame,
                               outcome_colname = "clade",
                               method = NULL,
                               collapse_corr_feats = TRUE,
                               group_neg_corr = TRUE,
                               to_numeric = FALSE,
                               remove_var = "zv")

pan_preproc_df <- pan_preproc$dat_transformed
colnames(pan_preproc_df)[grep("grp", colnames(pan_preproc_df))] <- gsub("$", "_pan", colnames(pan_preproc_df)[grep("grp", colnames(pan_preproc_df))])

struct_df <- read_delim("../../data/pan_struct_mat.tsv")

struct_df_sub <- struct_df %>%
  mutate(variant = gsub("_$", "", variant)) %>%
  column_to_rownames("variant")

struct_df_merge <- t(struct_df_sub)

if(sum(rownames(metadata_sub) %in% rownames(struct_df_merge)) != length(rownames(metadata_sub))){
  stop("mismatch between pheno and geno contents")
}

index <- match(rownames(metadata_sub), rownames(struct_df_merge))

struct_df_ordered <- struct_df_merge[index, , drop = FALSE]

if(sum(rownames(metadata_sub) == rownames(struct_df_ordered)) != length(rownames(metadata_sub))){
  stop("mismatch between pheno and geno contents")
}

struct_df_frame <- cbind(metadata_sub,
                         struct_df_ordered)

struct_preproc <- preprocess_data(dataset = struct_df_frame,
                                  outcome_colname = "clade",
                                  method = NULL,
                                  collapse_corr_feats = TRUE,
                                  group_neg_corr = TRUE,
                                  to_numeric = FALSE,
                                  remove_var = "zv")

struct_preproc_df <- struct_preproc$dat_transformed
colnames(struct_preproc_df)[grep("grp", colnames(struct_preproc_df))] <- gsub("$", "_struct", colnames(struct_preproc_df)[grep("grp", colnames(struct_preproc_df))])

geno_df <- bind_cols(metadata_sub,
                     core_preproc_df[,-1],
                     gene_preproc_df[,-1],
                     pan_preproc_df[,-1],
                     struct_preproc_df[,-1]) %>%
  mutate(clade = NULL)

geno_pivot <- as.data.frame(t(geno_df)) %>%
  rownames_to_column("variant")

write_delim(geno_pivot,
            "data/combined_mat.tsv")

geno_preprocessed <- list("core" = core_preproc,
                          "gene" = gene_preproc,
                          "pan" = pan_preproc,
                          "struct" = struct_preproc)

save(geno_preprocessed,
     file = "data/geno_frames_preprocessed.RData")

#PHENOTYPE ----
pheno_dirs <- list.files("../../2023_01_05_snakemake_sift_core_analysis/severity_core_sift/data/pheno",
                         full.names = TRUE)
pheno_path <- unlist(lapply(pheno_dirs, function(x){list.files(x,
                                                               pattern = "*.tsv",
                                                               full.names = TRUE)}))
pheno <- read_delim(pheno_path[1])
colnames(pheno)[2] <- gsub("\\.tsv",
                           "",
                           paste0(str_split(pheno_path[1],
                                            "/",
                                            simplify = TRUE)[,7],
                                  ".",
                                  str_split(pheno_path[1],
                                            "/",
                                            simplify = TRUE)[,8]))

for(i in 2:length(pheno_path)){
  
  a <- read_delim(pheno_path[i])
  colnames(a)[2]<- gsub("\\.tsv",
                        "",
                        paste0(str_split(pheno_path[i],
                                         "/",
                                         simplify = TRUE)[,7],
                               ".",
                               str_split(pheno_path[i],
                                         "/",
                                         simplify = TRUE)[,8]))
  
  pheno <- full_join(pheno,
                     a,
                     by = "genome_id")
  
}

colnames(pheno)[grep("q2_3", colnames(pheno))] <- gsub("q2_3", "q23", colnames(pheno)[grep("q2_3", colnames(pheno))])

pheno_chr <- pheno  %>%
  mutate(across(!genome_id,  ~replace(.x, .x == 0, "not_severe")),
         across(!genome_id,  ~replace(.x, .x == 1, "severe")))

write_csv(pheno_chr,
          "data/pheno_full.csv")