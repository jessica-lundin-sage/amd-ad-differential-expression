## ============================================================================
## UNIFIED DE PIPELINE — MSBB / ROSMAP / MAYO  (SEX AS COVARIATE)
## Consolidates DE_MSBB_FINAL_diagnosis_temp.R, DE_ROSMAP_FINAL_models.R and
## DE_MAYO_FINAL_models.R into one engine driven by per-cohort config lists.
##
## Each original script did the same thing end to end:
##   read a cqn'd dge_cqn/metadata object from Synapse -> cohort-specific
##   metadata prep (diagnosis derivation, transforms) -> set factors -> build
##   a cell-means (~0 + group + covariates + random effects) formula ->
##   voomWithDreamWeights -> contrasts (AD vs CT, by tissue) -> dream ->
##   eBayes -> topTable per contrast -> save an .rds and synStore it.
##
## SEX HANDLING (changed from the sex-stratified version):
##   - `group` is now diagnosis x tissue only (e.g. "AD2.STG"), NOT
##     diagnosis x sex x tissue.
##   - sex enters the model as a fixed-effect covariate via `+ sex` in each
##     cfg$formula_str. Males and females are fit together.
##   - Contrasts are one per tissue (AD vs CT, adjusted for sex), e.g.
##     "AD2_vs_CT2.STG", instead of one per sex x tissue.
##   - cfg$sex_levels is kept only to set the factor level order (the first
##     level is the reference, so the sex coefficient is "sexfemale").
##
## What differed between cohorts (and is now just config, not code):
##   - input Synapse ID / output parent Synapse ID / output filename
##   - tissue levels, diagnosis variable + levels, case/control labels
##   - the model formula (covariates + which variable is the random-effect batch)
##   - cohort-specific metadata prep (MSBB: Braak/CERAD/CDR diagnosis derivation,
##     age + PMI transforms, PC-tissue exclusion; MAYO: diagnosis text -> diag2;
##     ROSMAP: none needed beyond factoring final_batch)
##
## To run a cohort: edit the Synapse IDs in the matching cfg block below (or
## override them right before calling run_DE_pipeline(), see bottom of file),
## then call run_DE_pipeline(msbb_cfg) / run_DE_pipeline(rosmap_cfg) /
## run_DE_pipeline(mayo_cfg).
##
## NOTE ON RESULTS NAMING: each run returns DE_save$results, a named list keyed
## by contrast name (e.g. "AD2_vs_CT2.STG"). DE_save$results_sex holds the
## topTable for the sex coefficient itself (female vs male, all tissues).
## ============================================================================

pacman::p_load(tidyverse, limma, edgeR, biomaRt, DESeq2, vsn, sva, pamr)
pacman::p_load(synapser, dplyr, purrr, readr, lubridate, stringr, tibble, ggplot2)
pacman::p_load(variancePartition, BiocParallel, png, grid, knitr)
pacman::p_load(dplyr, ggplot2, viridis, patchwork, matrixStats, stringr, forcats, plotly, knitr, kableExtra)
pacman::p_load(cqn, ggrepel, RColorBrewer)

library(synapser)
synLogin()

## ----------------------------------------------------------------------------
## Shared helpers
## ----------------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

# Rank-based inverse normal transformation (used for ageDeath, which is
# heavily truncated at "90+").
rank_inverse_normal <- function(x, k = 3 / 8) {
  r <- rank(x, na.last = "keep", ties.method = "average")
  n <- sum(!is.na(x))
  qnorm((r - k) / (n - 2 * k + 1))
}

# Builds the AD_vs_CT-style contrasts (case - control) for every tissue.
# `group` is diagnosis x tissue, so each contrast is sex-adjusted (sex is a
# covariate in the formula) rather than sex-specific. `case_label`/
# `control_label` let the contrast NAME differ from the factor LEVEL.
build_group_contrasts <- function(tissue_levels,
                                  case_level, control_level,
                                  case_label = case_level,
                                  control_label = control_level) {
  exprs <- sprintf(
    "(group%s.%s - group%s.%s)",
    case_level, tissue_levels,
    control_level, tissue_levels
  )
  nms <- sprintf(
    "%s_vs_%s.%s",
    case_label, control_label, tissue_levels
  )
  setNames(exprs, nms)
}

# Minimal defensive check that a cfg has everything the engine needs.
validate_cfg <- function(cfg) {
  required <- c(
    "name", "input_synid", "output_parent_synid", "output_filename",
    "tissue_levels", "diagnosis_var", "diagnosis_levels",
    "case_level", "control_level", "formula_str", "prep_metadata", "n_workers"
  )
  missing <- setdiff(required, names(cfg))
  if (length(missing) > 0) {
    stop("cfg for '", cfg$name %||% "?", "' is missing: ", paste(missing, collapse = ", "))
  }
  if (!grepl("\\bsex\\b", cfg$formula_str)) {
    warning("cfg for '", cfg$name, "': formula_str does not include `sex` as a covariate")
  }
}

## ----------------------------------------------------------------------------
## Shared engine — identical steps for every cohort
## ----------------------------------------------------------------------------

run_DE_sexadj_pipeline <- function(cfg, qc_plots = TRUE) {
  validate_cfg(cfg)
  message("== Running DE pipeline (sex as covariate): ", cfg$name, " ==")
  
  # 1. Pull the cqn'd counts + metadata object from Synapse
  file_entity <- synGet(cfg$input_synid)
  cqn_obj <- readRDS(file_entity$path)
  md_sv   <- cqn_obj$metadata
  dge_cqn <- cqn_obj$dge_cqn
  
  # 2. Cohort-specific metadata prep (may add columns and/or drop rows,
  #    e.g. MSBB excluding the PC tissue and deriving diagnosis/age/PMI terms)
  md_sv <- cfg$prep_metadata(md_sv, qc_plots = qc_plots)
  
  # 3. Keep only samples present in the (possibly filtered) metadata, then
  #    align dge_cqn <-> md_sv exactly as all three originals did
  keep_samples <- colnames(dge_cqn) %in% md_sv$specimenID
  dge_cqn <- dge_cqn[, keep_samples]
  dge_cqn$E <- dge_cqn$E[, keep_samples, drop = FALSE]
  
  rownames(md_sv) <- md_sv$specimenID
  specID <- colnames(dge_cqn$counts)
  geneID <- rownames(dge_cqn$counts)
  colnames(dge_cqn$offset) <- specID; rownames(dge_cqn$offset) <- geneID
  colnames(dge_cqn$E)      <- specID; rownames(dge_cqn$E)      <- geneID
  md_sv <- md_sv[specID, , drop = FALSE]
  
  stopifnot(identical(colnames(dge_cqn$counts), colnames(dge_cqn$offset)))
  stopifnot(identical(colnames(dge_cqn$counts), rownames(md_sv)))
  
  # 4. Factors + cell-means group variable (diagnosis x tissue; sex is a
  #    covariate in the formula, not part of the group)
  md_sv$tissue2 <- factor(md_sv$tissue2, levels = cfg$tissue_levels)
  md_sv$sex     <- factor(md_sv$sex,     levels = cfg$sex_levels %||% c("male", "female"))
  md_sv[[cfg$diagnosis_var]] <- factor(md_sv[[cfg$diagnosis_var]], levels = cfg$diagnosis_levels)
  
  md_sv$group <- interaction(md_sv[[cfg$diagnosis_var]], md_sv$tissue2, drop = TRUE)
  group <- md_sv$group
  # some formula terms are resolved via the global env as a fallback, matching
  # the original scripts' `environment(form_final0) <- globalenv()` pattern
  assign("group", group, envir = globalenv())
  
  if (qc_plots) {
    print(table(md_sv$group))
    print(table(md_sv$group, md_sv$sex, useNA = "ifany"))
  }
  
  form_final0 <- as.formula(cfg$formula_str)
  environment(form_final0) <- globalenv()
  md_sv <- as.data.frame(md_sv)
  
  # 5. voom weights + dream fit
  param_n <- MulticoreParam(workers = cfg$n_workers, progressbar = TRUE)
  
  voom_gene_expression <- voomWithDreamWeights(
    counts  = dge_cqn,
    formula = form_final0,
    data    = md_sv,
    plot    = qc_plots,
    BPPARAM = param_n
  )
  #voom_gene_expression$E <- as.matrix(dge_cqn$E)
  voom_gene_expression$E <- as.matrix(dge_cqn$E)[names(voom_gene_expression$voom.xy$x) %||% rownames(dge_cqn$E),
                                                 colnames(voom_gene_expression), drop = FALSE]
  contrasts_vec <- build_group_contrasts(
    tissue_levels  = cfg$tissue_levels,
    case_level     = cfg$case_level,
    control_level  = cfg$control_level,
    case_label     = cfg$case_label %||% cfg$case_level,
    control_label  = cfg$control_label %||% cfg$control_level
  )
  
  L_group <- makeContrastsDream(form_final0, md_sv, contrasts = contrasts_vec)
  if (qc_plots) plotContrasts(L_group)
  
  fit_contrasts <- variancePartition::dream(
    exprObj = voom_gene_expression,
    formula = form_final0,
    data    = md_sv,
    L       = L_group,
    BPPARAM = param_n
  )
  ebayes_fit <- eBayes(fit_contrasts)
  
  # 6. One topTable per contrast, collected into a named list
  #    (selected by name, since sex/covariate coefs also live in the fit)
  de_results <- setNames(
    lapply(colnames(L_group), function(nm) topTable(ebayes_fit, coef = nm, number = Inf, confint = TRUE)),
    colnames(L_group)
  )
  
  # Sex main effect (female vs male, adjusted for group + covariates)
  sex_coef <- grep("^sex", colnames(coef(ebayes_fit)), value = TRUE)
  results_sex <- setNames(
    lapply(sex_coef, function(nm) topTable(ebayes_fit, coef = nm, number = Inf, confint = TRUE)),
    sex_coef
  )
  
  if (qc_plots) {
    for (nm in names(de_results)) {
      tt <- de_results[[nm]]
      plot(tt$logFC, -log10(tt$adj.P.Val), main = nm, xlab = "logFC", ylab = "-log10(adj.P.Val)")
    }
  }
  
  DE_save <- list(
    metadata      = md_sv,
    vobj_expr     = voom_gene_expression,
    dge_cqn       = dge_cqn,
    ebayes        = ebayes_fit,
    fit_contrasts = fit_contrasts,
    results       = de_results,
    results_sex   = results_sex
  )
  
  # 7. Save + upload to the cohort's Synapse parent folder
  saveRDS(DE_save, file = cfg$output_filename)
  invisible(synLogin(silent = TRUE))
  synStore(File(path = cfg$output_filename, parent = cfg$output_parent_synid))
  message("Saved '", cfg$output_filename, "' -> ", cfg$output_parent_synid)
  
  # optional secondary save under a different name/parent (off by default) —
  # this replaces MSBB's old ad-hoc "_dx_v9" comparison save
  if (isTRUE(cfg$save_diagnostic_copy) && !is.null(cfg$diagnostic_filename)) {
    saveRDS(DE_save, file = cfg$diagnostic_filename)
    synStore(File(path = cfg$diagnostic_filename,
                  parent = cfg$diagnostic_parent_synid %||% cfg$output_parent_synid), 
                  used = cfg$used, 
                  executed = cfg$executed, 
                  activityName = cfg$activityName, 
                  activityDescription = cfg$activityDescription
             )
    message("Also saved diagnostic copy '", cfg$diagnostic_filename, "'")
  }
  
  invisible(DE_save)
}

## ----------------------------------------------------------------------------
## Cohort-specific metadata prep
## ----------------------------------------------------------------------------

msbb_prep_metadata <- function(md_sv, qc_plots = TRUE) {
  md_sv <- md_sv[md_sv$tissue2 != "PC", ]
  
  # Rank-based inverse normal transform for ageDeath, to account for the
  # large number of subjects truncated at "90+"
  md_sv$ageDeath_invnorm <- rank_inverse_normal(md_sv$ageDeath)
  
  md_sv$PMI_log <- log(md_sv$PMI)
  
  if (qc_plots) {
    hist(md_sv$ageDeath_invnorm, main = "ageDeath_invnorm")
    hist(md_sv$RIN, main = "RIN")
    hist(md_sv$PMI, main = "PMI")
    hist(md_sv$PMI_log, main = "log(PMI)")
  }
  
  if (qc_plots) {
    print(table(md_sv$diag2, useNA = "always"))
    print(table(md_sv$diag2, md_sv$sex, useNA = "always"))
  }
  
  md_sv
}

rosmap_prep_metadata <- function(md_sv, qc_plots = TRUE) {
  md_sv$final_batch <- factor(md_sv$final_batch)
  md_sv$ageDeath_invnorm <- rank_inverse_normal(md_sv$ageDeath)
  md_sv$PMI_log <- log(md_sv$PMI)
  md_sv
}

mayo_prep_metadata <- function(md_sv, qc_plots = TRUE) {
  md_sv$ageDeath_invnorm <- rank_inverse_normal(md_sv$ageDeath)
  md_sv$diag2 <- NA_character_
  md_sv$diag2[md_sv$diagnosis == "Alzheimer Disease"]              <- "AD2"
  md_sv$diag2[md_sv$diagnosis == "control"]                        <- "CT2"
  md_sv$diag2[md_sv$diagnosis == "pathological aging"]             <- "PathAg2"
  md_sv$diag2[md_sv$diagnosis == "progressive supranuclear palsy"] <- "PSP2"
  if (qc_plots) print(table(md_sv$diag2))
  md_sv
}

## ----------------------------------------------------------------------------
## Cohort configs — this is the part you edit per run
## ----------------------------------------------------------------------------

msbb_cfg_sexadj <- list(
  name                 = "MSBB",
  input_synid          = "syn77540009",  # MSBB_md_counts_cqn_FINAL.rds - new syn folder
  output_parent_synid  = "syn77539791",
  output_filename      = "amp-ad-de_MSBB_neuro.path.diag_sex_adj.rds",
  tissue_levels        = c("FP", "STG", "IFG", "PG"),
  sex_levels           = c("male", "female"),  # factor order only; male = reference
  diagnosis_var        = "diag2",
  diagnosis_levels     = c("CT2", "AD2", "OTHER2"),
  case_level           = "AD2",
  control_level        = "CT2",
  case_label           = "AD2",
  control_label        = "CT2",
  formula_str          = "~ 0 + group + sex + PMI_log + RIN + apoe4Status + ageDeath_invnorm + PC1_metrics + PC2_metrics + PC3_metrics + PC4_metrics + (1|sequencingBatch) + (1|individualID)",
  prep_metadata        = msbb_prep_metadata,
  n_workers            = 4,
  used                 = c("syn77540009"), 
  executed             = "https://github.com/jessica-lundin-sage/amd-ad-differential-expression", 
  activityName         = "Reprocessed RNAseq differential expression files", 
  activityDescription  = "Metadata, counts, differential expression results"
)

rosmap_cfg_sexadj <- list(
  name                = "ROSMAP",
  input_synid         = "syn77559937",  # ROSMAP_md_counts_cqn_DLPFC_CN_PCC_FINAL.rds - new syn folder
  output_parent_synid = "syn77539791",
  output_filename     = "amp-ad-de_ROSMAP_neuro.path.diag_sex_adj.rds",
  tissue_levels       = c("DLPFC", "PCC", "CN"),
  sex_levels          = c("male", "female"),  # factor order only; male = reference
  diagnosis_var       = "diag2",
  diagnosis_levels    = c("CT2", "AD2", "OTHER2"),
  case_level          = "AD2",
  control_level       = "CT2",
  formula_str         = "~ 0 + group + sex + apoe4Status + ageDeath_invnorm + PMI_log + RIN + PC1_metrics + PC2_metrics + PC3_metrics + PC4_metrics + (1|final_batch) + (1|individualID)",
  prep_metadata       = rosmap_prep_metadata,
  n_workers           = 4,
  used                = c("syn77559937"), 
  executed            = "https://github.com/jessica-lundin-sage/amd-ad-differential-expression", 
  activityName        = "Reprocessed RNAseq differential expression files", 
  activityDescription = "Metadata, counts, differential expression results"
)

mayo_cfg_sexadj <- list(
  name                = "MAYO",
  input_synid         = "syn77548692",  # MAYO_md_counts_cqn_FINAL.rds - new syn folder
  output_parent_synid = "syn77539791",
  output_filename     = "amp-ad-de_MAYO_neuro.path.diag_sex_adj.rds",
  tissue_levels       = c("CER", "TCX"),
  sex_levels          = c("male", "female"),  # factor order only; male = reference
  diagnosis_var       = "diag2",
  diagnosis_levels    = c("CT2", "AD2", "PathAg2", "PSP2"),
  case_level          = "AD2",
  control_level       = "CT2",
  formula_str         = "~ 0 + group + sex + apoe4Status + ageDeath_invnorm + RIN + PC1_metrics + PC2_metrics + PC3_metrics + (1|flowcell) + (1|individualID)",
  prep_metadata       = mayo_prep_metadata,
  n_workers           = 4,
  used                = c("syn77548692"), 
  executed            = "https://github.com/jessica-lundin-sage/amd-ad-differential-expression", 
  activityName        = "Reprocessed RNAseq differential expression files", 
  activityDescription = "Metadata, counts, differential expression results"
)


rosmap_cfg_sexadj_plus_clinical <- list(
  name                = "ROSMAP",
  input_synid         = "syn77559937",  # ROSMAP_md_counts_cqn_DLPFC_CN_PCC_FINAL.rds - new syn folder
  output_parent_synid = "syn77539791",
  output_filename     = "amp-ad-de_ROSMAP_neuro.path.diag.plus.clinical_sex_adj.rds",
  tissue_levels       = c("DLPFC", "PCC", "CN"),
  sex_levels          = c("male", "female"),  # factor order only; male = reference
  diagnosis_var       = "diagnosis",
  diagnosis_levels    = c("CT", "AD", "OTHER"),
  case_level          = "AD",
  control_level       = "CT",
  formula_str         = "~ 0 + group + sex + apoe4Status + ageDeath_invnorm + PMI_log + RIN + PC1_metrics + PC2_metrics + PC3_metrics + PC4_metrics + (1|final_batch) + (1|individualID)",
  prep_metadata       = rosmap_prep_metadata,
  n_workers           = 4,
  used                = c("syn77559937"), 
  executed            = "https://github.com/jessica-lundin-sage/amd-ad-differential-expression", 
  activityName        = "Reprocessed RNAseq differential expression files", 
  activityDescription = "Metadata, counts, differential expression results"
)

msbb_cfg_sexadj_plus_clinical  <- list(
  name                 = "MSBB",
  input_synid          = "syn77540009",  # MSBB_md_counts_cqn_FINAL.rds - new syn folder
  output_parent_synid  = "syn77539791",
  output_filename      = "amp-ad-de_MSBB_neuro.path.diag.plus.clinical_sex_adj.rds",
  tissue_levels        = c("FP", "STG", "IFG", "PG"),
  sex_levels           = c("male", "female"),  # factor order only; male = reference
  diagnosis_var        = "diagnosis",
  diagnosis_levels     = c("CT", "AD", "OTHER"),
  case_level           = "AD",
  control_level        = "CT",
  case_label           = "AD",
  control_label        = "CT",
  formula_str          = "~ 0 + group + sex + PMI_log + RIN + apoe4Status + ageDeath_invnorm + PC1_metrics + PC2_metrics + PC3_metrics + PC4_metrics + (1|sequencingBatch) + (1|individualID)",
  prep_metadata        = msbb_prep_metadata,
  n_workers            = 4,
  used                 = c("syn77540009"), 
  executed             = "https://github.com/jessica-lundin-sage/amd-ad-differential-expression", 
  activityName         = "Reprocessed RNAseq differential expression files", 
  activityDescription  = "Metadata, counts, differential expression results"
)
## ----------------------------------------------------------------------------
## RUN
## Edit Synapse IDs above (or override a field right before calling, e.g.
##   msbb_cfg$input_synid <- "synNEWID"
##   msbb_cfg$output_parent_synid <- "synNEWPARENT"
## ) then run one of:
## ----------------------------------------------------------------------------

result_msbb   <- run_DE_sexadj_pipeline(msbb_cfg_sexadj)
result_rosmap <- run_DE_sexadj_pipeline(rosmap_cfg_sexadj)
result_mayo   <- run_DE_sexadj_pipeline(mayo_cfg_sexadj)

result_rosmap_plus_clinical <- run_DE_sexadj_pipeline(rosmap_cfg_sexadj_plus_clinical)
result_msbb2_plus_clinical <- run_DE_sexadj_pipeline(msbb_cfg_sexadj_plus_clinical)
