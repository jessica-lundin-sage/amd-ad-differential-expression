## ============================================================================
## UNIFIED DE PIPELINE — MSBB / ROSMAP / MAYO
## Sex-STRATIFIED  and  sex-ADJUSTED  models from one engine.
##
## Consolidates DE_MSBB_FINAL_diagnosis_temp.R, DE_ROSMAP_FINAL_models.R and
## DE_MAYO_FINAL_models.R into one engine driven by per-cohort config lists.
##
## Each run does:
##   read a cqn'd dge_cqn/metadata object from Synapse -> cohort-specific
##   metadata prep (diagnosis derivation, transforms) -> set factors -> build
##   a cell-means (~0 + group + covariates + random effects) formula ->
##   voomWithDreamWeights -> contrasts (AD vs CT) -> dream -> eBayes ->
##   topTable per contrast -> merge gene info -> write CSV + synStore it
##   (plus the normalized log2 counts CSV, if configured).
##
## SEX HANDLING — set by cfg$sex_mode:
##   "stratified" (default)
##       group = diagnosis x sex x tissue (e.g. "AD2.male.STG")
##       one contrast per sex x tissue, e.g. "AD_vs_CT_male.STG"
##       -> Comparison column: "AD_male_STG - CT_male_STG"
##       formula_str must NOT contain `sex` (it is inside `group`)
##   "adjusted"
##       group = diagnosis x tissue (e.g. "AD2.STG")
##       sex enters as a fixed-effect covariate (`+ sex` in formula_str)
##       one contrast per tissue, e.g. "AD_vs_CT.STG"
##       -> Comparison column: "AD_STG - CT_STG"
##       optionally also writes the sex main effect (female - male, adjusted
##       for group + covariates) as Comparison "female - male"
##       (cfg$include_sex_effect, default TRUE)
##
## The easiest way to get a sex-adjusted run is make_sexadj_cfg(<stratified
## cfg>): it copies the cfg, sets sex_mode = "adjusted", inserts `+ sex` after
## `group` in the formula, and renames output files "sex_stratified" ->
## "sex_adj". See the RUN section at the bottom.
##
## Cohort differences live in config only:
##   - input / output Synapse IDs and filenames
##   - tissue levels, diagnosis variable + levels, case/control labels
##   - model formula (covariates + random-effect batch variable)
##   - cohort-specific metadata prep (MSBB: PC exclusion, PG -> PHG, age + PMI
##     transforms; ROSMAP: CN -> ACC, final_batch factor; MAYO: CER -> CBE,
##     diagnosis text -> diag2)
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

# Builds AD_vs_CT-style contrasts (case - control).
#   sex_mode = "stratified": one per sex x tissue, "tissue outer loop, sex
#                            inner loop (male, female)" order.
#   sex_mode = "adjusted"  : one per tissue (sex is a covariate).
# `case_label`/`control_label` let the contrast NAME differ from the factor
# LEVEL (levels "AD2"/"CT2", contrast names "AD_vs_CT_...").
build_group_contrasts <- function(tissue_levels, sex_levels,
                                  case_level, control_level,
                                  case_label = case_level,
                                  control_label = control_level,
                                  sex_mode = c("stratified", "adjusted")) {
  sex_mode <- match.arg(sex_mode)
  
  if (sex_mode == "stratified") {
    combos <- expand.grid(sex = sex_levels, tissue = tissue_levels, stringsAsFactors = FALSE)
    combos <- combos[order(match(combos$tissue, tissue_levels), match(combos$sex, sex_levels)), ]
    exprs <- sprintf(
      "(group%s.%s.%s - group%s.%s.%s)",
      case_level, combos$sex, combos$tissue,
      control_level, combos$sex, combos$tissue
    )
    nms <- sprintf("%s_vs_%s_%s.%s", case_label, control_label, combos$sex, combos$tissue)
  } else {
    exprs <- sprintf(
      "(group%s.%s - group%s.%s)",
      case_level, tissue_levels,
      control_level, tissue_levels
    )
    nms <- sprintf("%s_vs_%s.%s", case_label, control_label, tissue_levels)
  }
  setNames(exprs, nms)
}

# Minimal defensive check that a cfg has everything the engine needs.
validate_cfg <- function(cfg) {
  required <- c(
    "name", "input_synid", "output_parent_synid", "output_filename",
    "tissue_levels", "sex_levels", "diagnosis_var", "diagnosis_levels",
    "case_level", "control_level", "formula_str", "prep_metadata", "n_workers"
  )
  missing <- setdiff(required, names(cfg))
  if (length(missing) > 0) {
    stop("cfg for '", cfg$name %||% "?", "' is missing: ", paste(missing, collapse = ", "))
  }
  
  sex_mode <- cfg$sex_mode %||% "stratified"
  if (!sex_mode %in% c("stratified", "adjusted")) {
    stop("cfg for '", cfg$name, "': sex_mode must be 'stratified' or 'adjusted'")
  }
  has_sex <- grepl("\\bsex\\b", cfg$formula_str)
  if (sex_mode == "adjusted" && !has_sex) {
    warning("cfg for '", cfg$name, "': sex_mode = 'adjusted' but formula_str has no `sex` covariate")
  }
  if (sex_mode == "stratified" && has_sex) {
    stop("cfg for '", cfg$name, "': sex_mode = 'stratified' but formula_str includes `sex` ",
         "(sex is already in `group`; the design would be rank-deficient)")
  }
}

# Turn a sex-stratified cfg into the matching sex-adjusted cfg.
#   - sex_mode            -> "adjusted"
#   - formula_str         -> inserts "+ sex" right after "group"
#   - output_filename     -> "sex_stratified" replaced by "sex_adj"
#                            (or "_sex_adj" appended before the extension)
#   - out_counts_filename -> NULL by default: the normalized counts are the
#                            cqn E matrix and don't depend on the model, so
#                            they'd duplicate the stratified run's file.
#                            Pass keep_counts = TRUE to write them anyway.
# Any extra named args override fields, e.g.
#   make_sexadj_cfg(msbb_cfg, output_filename = "my_name.csv")
make_sexadj_cfg <- function(cfg, keep_counts = FALSE, ...) {
  new <- cfg
  new$sex_mode <- "adjusted"
  
  if (!grepl("\\bsex\\b", new$formula_str)) {
    new$formula_str <- sub("\\bgroup\\b", "group + sex", new$formula_str)
  }
  
  rename_out <- function(f) {
    if (is.null(f)) return(f)
    if (grepl("sex_stratified", f)) return(sub("sex_stratified", "sex_adj", f))
    sub("(\\.[^.]+)$", "_sex_adj\\1", f)
  }
  new$output_filename <- rename_out(new$output_filename)
  
  if (keep_counts) {
    new$out_counts_filename <- rename_out(new$out_counts_filename)
  } else {
    new$out_counts_filename <- NULL
  }
  
  overrides <- list(...)
  for (nm in names(overrides)) new[[nm]] <- overrides[[nm]]
  new
}

## ----------------------------------------------------------------------------
## Shared engine — identical steps for every cohort / sex_mode
## ----------------------------------------------------------------------------

run_DE_pipeline <- function(cfg, qc_plots = TRUE) {
  validate_cfg(cfg)
  sex_mode <- cfg$sex_mode %||% "stratified"
  message("== Running DE pipeline: ", cfg$name, " (sex ", sex_mode, ") ==")
  
  # 1. Pull the cqn'd counts + metadata object from Synapse
  file_entity <- synGet(cfg$input_synid)
  cqn_obj <- readRDS(file_entity$path)
  md_sv   <- cqn_obj$metadata
  dge_cqn <- cqn_obj$dge_cqn
  
  # 2. Cohort-specific metadata prep (may add columns and/or drop rows)
  md_sv <- cfg$prep_metadata(md_sv)
  
  # centering vars so scales are not so different to help with convergence (small number excluded due to this)
  num_covs <- c("RIN", "ageDeath_invnorm",
                grep("^PC\\d+_metrics$", colnames(md_sv), value = TRUE))
  num_covs <- intersect(num_covs, all.vars(form_final0))
  md_sv[num_covs] <- lapply(md_sv[num_covs], function(x) as.numeric(scale(x)))
  
  # 3. Keep only samples present in the (possibly filtered) metadata, then
  #    align dge_cqn <-> md_sv
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
  
  # 4. Factors + cell-means group variable
  md_sv$tissue2 <- factor(md_sv$tissue2, levels = cfg$tissue_levels)
  md_sv$sex     <- factor(md_sv$sex,     levels = cfg$sex_levels)  # 1st level = reference in adjusted mode
  md_sv[[cfg$diagnosis_var]] <- factor(md_sv[[cfg$diagnosis_var]], levels = cfg$diagnosis_levels)
  
  if (sex_mode == "stratified") {
    md_sv$group <- interaction(md_sv[[cfg$diagnosis_var]], md_sv$sex, md_sv$tissue2, drop = TRUE)
  } else {
    md_sv$group <- interaction(md_sv[[cfg$diagnosis_var]], md_sv$tissue2, drop = TRUE)
  }
  group <- md_sv$group
  # some formula terms are resolved via the global env as a fallback, matching
  # the original scripts' `environment(form_final0) <- globalenv()` pattern
  assign("group", group, envir = globalenv())
  
  if (qc_plots) {
    print(table(md_sv$group))
    if (sex_mode == "adjusted") print(table(md_sv$group, md_sv$sex, useNA = "ifany"))
  }
  
  form_final0 <- as.formula(cfg$formula_str)
  environment(form_final0) <- globalenv()
  md_sv <- as.data.frame(md_sv)
  
  # Set standard global seed for the main session
  set.seed(123)
  
  # 5. voom weights + dream fit
  param_n <- MulticoreParam(workers = cfg$n_workers, RNGseed = 123, progressbar = TRUE)
  
  voom_gene_expression <- voomWithDreamWeights(
    counts  = dge_cqn,
    formula = form_final0,
    data    = md_sv,
    plot    = qc_plots,
    BPPARAM = param_n
  )
  voom_gene_expression$E <- as.matrix(dge_cqn$E)[names(voom_gene_expression$voom.xy$x) %||% rownames(dge_cqn$E),
                                                 colnames(voom_gene_expression), drop = FALSE]
  
  contrasts_vec <- build_group_contrasts(
    tissue_levels  = cfg$tissue_levels,
    sex_levels     = cfg$sex_levels,
    case_level     = cfg$case_level,
    control_level  = cfg$control_level,
    case_label     = cfg$case_label %||% cfg$case_level,
    control_label  = cfg$control_label %||% cfg$control_level,
    sex_mode       = sex_mode
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
  
  # 6. One topTable per contrast, selected by NAME (covariate coefs also live
  #    in the fit, so positional indexing isn't safe in adjusted mode)
  de_results <- setNames(
    lapply(colnames(L_group), function(nm) topTable(ebayes_fit, coef = nm, number = Inf, confint = TRUE)),
    colnames(L_group)
  )
  
  # Sex main effect (adjusted mode only): female vs male, adjusted for
  # group + covariates
  results_sex <- list()
  if (sex_mode == "adjusted") {
    sex_coef <- grep("^sex", colnames(coef(ebayes_fit)), value = TRUE)
    results_sex <- setNames(
      lapply(sex_coef, function(nm) topTable(ebayes_fit, coef = nm, number = Inf, confint = TRUE)),
      sex_coef
    )
  }
  
  if (qc_plots) {
    for (nm in names(de_results)) {
      tt <- de_results[[nm]]
      plot(tt$logFC, -log10(tt$adj.P.Val), main = nm, xlab = "logFC", ylab = "-log10(adj.P.Val)")
    }
  }
  
  # 7. Rename contrasts to "case - control" Comparison labels.
  #    Each pattern only matches its own format, so the order doesn't matter.
  res_x <- de_results
  nm <- names(res_x)
  nm <- sub("^(AD2?)_vs_(CT2?)_(\\w+)\\.(\\w+)$", "\\1_\\3_\\4 - \\2_\\3_\\4", nm)  # sex-stratified
  nm <- sub("^(AD2?)_vs_(CT2?)\\.(\\w+)$",        "\\1_\\3 - \\2_\\3",         nm)  # sex-adjusted (tissue only)
  names(res_x) <- nm
  
  if (sex_mode == "adjusted" && isTRUE(cfg$include_sex_effect %||% TRUE) && length(results_sex) > 0) {
    ref_sex <- cfg$sex_levels[1]
    sex_nm  <- sprintf("%s - %s", sub("^sex", "", names(results_sex)), ref_sex)  # e.g. "female - male"
    res_x   <- c(res_x, setNames(results_sex, sex_nm))
  }
  
  # Stack into one long table + add gene info
  combined <- res_x |>
    lapply(rownames_to_column, var = "ensembl_gene_id") |>
    bind_rows(.id = "Comparison")
  
  syn_gene_file <- synapser::synGet("syn74907403") ## created using code from here: https://github.com/jaclynbeck-sage/sageRNAUtils/blob/main/R/gc_content_functions.R
  gene_info <- read.csv(syn_gene_file$path)
  combined <- merge(combined, gene_info, by = "ensembl_gene_id", all.x=T)
  
  # 8. Save DE results CSV and upload
  write.csv(combined, cfg$output_filename, row.names = FALSE)
  invisible(synLogin(silent = TRUE))
  synStore(File(cfg$output_filename, parent = cfg$output_parent_synid))
  message("Saved '", cfg$output_filename, "' -> ", cfg$output_parent_synid)
  
  # 9. Normalized log2 counts (skipped if out_counts_filename is NULL)
  if (!is.null(cfg$out_counts_filename)) {
    Normalized_log2_counts <- as.data.frame(voom_gene_expression$E)
    Normalized_log2_counts <- tibble::rownames_to_column(Normalized_log2_counts, var = "ensembl_gene_id")
    write.csv(Normalized_log2_counts, cfg$out_counts_filename, row.names = FALSE)
    invisible(synLogin(silent = TRUE))
    synStore(File(cfg$out_counts_filename,
                  parent = cfg$out_counts_parent_synid %||% cfg$output_parent_synid))
    message("Saved '", cfg$out_counts_filename, "' -> ",
            cfg$out_counts_parent_synid %||% cfg$output_parent_synid)
  }
  
  # Returned (not saved) so results can be inspected interactively
  invisible(list(
    cfg           = cfg,
    metadata      = md_sv,
    vobj_expr     = voom_gene_expression,
    ebayes        = ebayes_fit,
    fit_contrasts = fit_contrasts,
    results       = de_results,
    results_sex   = results_sex,
    combined      = combined
  ))
}

## ----------------------------------------------------------------------------
## Cohort-specific metadata prep
## ----------------------------------------------------------------------------

msbb_prep_metadata <- function(md_sv, qc_plots = TRUE) {
  md_sv <- md_sv[md_sv$tissue2 != "PC", ]
  md_sv[md_sv$tissue2 == "PG", "tissue2"] <- "PHG"
  
  # Rank-based inverse normal transform for ageDeath, to account for the
  # large number of subjects truncated at "90+"
  md_sv$ageDeath_invnorm <- rank_inverse_normal(md_sv$ageDeath)
  md_sv$PMI_log <- log(md_sv$PMI)
  
  if (qc_plots) {
    hist(md_sv$ageDeath_invnorm, main = "ageDeath_invnorm")
    hist(md_sv$RIN, main = "RIN")
    hist(md_sv$PMI, main = "PMI")
    hist(md_sv$PMI_log, main = "log(PMI)")
    print(table(md_sv$diagnosis, useNA = "always"))
    print(table(md_sv$diagnosis, md_sv$sex, useNA = "always"))
  }
  
  md_sv
}

rosmap_prep_metadata <- function(md_sv, qc_plots = TRUE) {
  md_sv[md_sv$tissue2 == "CN", "tissue2"] <- "ACC"
  md_sv$final_batch <- factor(md_sv$final_batch)
  md_sv$ageDeath_invnorm <- rank_inverse_normal(md_sv$ageDeath)
  md_sv$PMI_log <- log(md_sv$PMI)
  md_sv
}

mayo_prep_metadata <- function(md_sv, qc_plots = TRUE) {
  md_sv[md_sv$tissue2 == "CER", "tissue2"] <- "CBE"
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
## Cohort configs (sex-stratified) — this is the part you edit per run.
## Sex-adjusted versions are derived from these with make_sexadj_cfg() below.
## ----------------------------------------------------------------------------

msbb_cfg <- list(
  name                    = "MSBB",
  sex_mode                = "stratified",
  input_synid             = "syn77540009",  # MSBB_md_counts_cqn_FINAL.rds - new syn folder
  output_parent_synid     = "syn77965053",
  output_filename         = "de_MSBB_neuro.path.diag_sex_stratified.csv",
  out_counts_parent_synid = "syn77870712",
  out_counts_filename     = "norm-counts_msbb.csv",
  tissue_levels           = c("FP", "STG", "IFG", "PHG"),
  sex_levels              = c("male", "female"),
  diagnosis_var           = "diag2",
  diagnosis_levels        = c("CT2", "AD2", "OTHER2"),
  case_level              = "AD2",
  control_level           = "CT2",
  case_label              = "AD",     # contrast names use AD/CT, not AD2/CT2
  control_label           = "CT",
  formula_str             = "~ 0 + group + PMI_log + RIN + apoe4Status + ageDeath_invnorm + PC1_metrics + PC2_metrics + PC3_metrics + PC4_metrics + (1|sequencingBatch) + (1|individualID)",
  prep_metadata           = msbb_prep_metadata,
  n_workers               = 4,
  used                    = c("syn77540009"),
  executed                = "https://github.com/jessica-lundin-sage/amd-ad-differential-expression",
  activityName            = "Reprocessed RNAseq differential expression files",
  activityDescription     = "Metadata, counts, differential expression results"
)

rosmap_cfg <- list(
  name                    = "ROSMAP",
  sex_mode                = "stratified",
  input_synid             = "syn77964695",  # ROSMAP_md_counts_cqn_DLPFC_CN_PCC_FINAL.rds - new syn folder
  output_parent_synid     = "syn77965053",
  output_filename         = "de_ROSMAP_neuro.path.diag_sex_stratified.csv",
  out_counts_parent_synid = "syn77870712",
  out_counts_filename     = "norm-counts_rosmap.csv",
  tissue_levels           = c("DLPFC", "PCC", "ACC"),
  sex_levels              = c("male", "female"),
  diagnosis_var           = "diag2",
  diagnosis_levels        = c("CT2", "AD2", "OTHER2"),
  case_level              = "AD2",
  control_level           = "CT2",
  case_label              = "AD",
  control_label           = "CT",
  formula_str             = "~ 0 + group + apoe4Status + ageDeath_invnorm + PMI_log + RIN + PC1_metrics + PC2_metrics + PC3_metrics + PC4_metrics + (1|final_batch) + (1|individualID)",
  prep_metadata           = rosmap_prep_metadata,
  n_workers               = 4,
  used                    = c("syn77559937"),
  executed                = "https://github.com/jessica-lundin-sage/amd-ad-differential-expression",
  activityName            = "Reprocessed RNAseq differential expression files",
  activityDescription     = "Metadata, counts, differential expression results"
)

mayo_cfg <- list(
  name                    = "MAYO",
  sex_mode                = "stratified",
  input_synid             = "syn77548692",  # MAYO_md_counts_cqn_FINAL.rds - new syn folder
  output_parent_synid     = "syn77965053",
  output_filename         = "de_MAYO_neuro.path.diag_sex_stratified.csv",
  out_counts_parent_synid = "syn77870712",
  out_counts_filename     = "norm-counts_mayo_test.csv",
  tissue_levels           = c("CBE", "TCX"),
  sex_levels              = c("male", "female"),
  diagnosis_var           = "diag2",
  diagnosis_levels        = c("CT2", "AD2", "PathAg2", "PSP2"),
  case_level              = "AD2",
  control_level           = "CT2",
  case_label              = "AD",
  control_label           = "CT",
  formula_str             = "~ 0 + group + apoe4Status + ageDeath_invnorm + RIN + PC1_metrics + PC2_metrics + PC3_metrics + (1|flowcell) + (1|individualID)",
  prep_metadata           = mayo_prep_metadata,
  n_workers               = 4,
  used                    = c("syn77548692"),
  executed                = "https://github.com/jessica-lundin-sage/amd-ad-differential-expression",
  activityName            = "Reprocessed RNAseq differential expression files",
  activityDescription     = "Metadata, counts, differential expression results"
)

msbb_cfg_plus_clinical <- list(
  name                    = "MSBB",
  sex_mode                = "stratified",
  input_synid             = "syn77540009",  # MSBB_md_counts_cqn_FINAL.rds - new syn folder
  output_parent_synid     = "syn77965053",
  output_filename         = "de_MSBB_neuro.path.diag.plus.clinical_sex_stratified.csv",
  out_counts_parent_synid = "syn77870712",
  out_counts_filename     = "norm-counts_msbb.csv",
  tissue_levels           = c("FP", "STG", "IFG", "PHG"),
  sex_levels              = c("male", "female"),
  diagnosis_var           = "diagnosis",
  diagnosis_levels        = c("CT", "AD", "OTHER"),
  case_level              = "AD",
  control_level           = "CT",
  case_label              = "AD",
  control_label           = "CT",
  formula_str             = "~ 0 + group + PMI_log + RIN + apoe4Status + ageDeath_invnorm + PC1_metrics + PC2_metrics + PC3_metrics + PC4_metrics + (1|sequencingBatch) + (1|individualID)",
  prep_metadata           = msbb_prep_metadata,
  n_workers               = 4,
  used                    = c("syn77540009"),
  executed                = "https://github.com/jessica-lundin-sage/amd-ad-differential-expression",
  activityName            = "Reprocessed RNAseq differential expression files",
  activityDescription     = "Metadata, counts, differential expression results"
)

rosmap_cfg_plus_clinical <- list(
  name                    = "ROSMAP",
  sex_mode                = "stratified",
  input_synid             = "syn77964695",  # ROSMAP_md_counts_cqn_DLPFC_CN_PCC_FINAL.rds - new syn folder
  output_parent_synid     = "syn77965053",
  output_filename         = "de_ROSMAP_neuro.path.diag.plus.clinical_sex_stratified.csv",
  out_counts_parent_synid = "syn77870712",
  out_counts_filename     = "norm-counts_rosmap.csv",
  tissue_levels           = c("DLPFC", "PCC", "ACC"),
  sex_levels              = c("male", "female"),
  diagnosis_var           = "diagnosis",
  diagnosis_levels        = c("CT", "AD", "OTHER"),
  case_level              = "AD",
  control_level           = "CT",
  case_label              = "AD",
  control_label           = "CT",
  formula_str             = "~ 0 + group + apoe4Status + ageDeath_invnorm + PMI_log + RIN + PC1_metrics + PC2_metrics + PC3_metrics + PC4_metrics + (1|final_batch) + (1|individualID)",
  prep_metadata           = rosmap_prep_metadata,
  n_workers               = 4,
  used                    = c("syn77559937"),
  executed                = "https://github.com/jessica-lundin-sage/amd-ad-differential-expression",
  activityName            = "Reprocessed RNAseq differential expression files",
  activityDescription     = "Metadata, counts, differential expression results"
)

## ----------------------------------------------------------------------------
## Sex-adjusted configs — derived from the stratified ones above
## (formula gets "+ sex", outputs renamed "..._sex_adj.csv", counts file
## skipped since it would be identical to the stratified run's)
## ----------------------------------------------------------------------------

msbb_cfg_sexadj                 <- make_sexadj_cfg(msbb_cfg)
rosmap_cfg_sexadj               <- make_sexadj_cfg(rosmap_cfg)
mayo_cfg_sexadj                 <- make_sexadj_cfg(mayo_cfg)
msbb_cfg_plus_clinical_sexadj   <- make_sexadj_cfg(msbb_cfg_plus_clinical)
rosmap_cfg_plus_clinical_sexadj <- make_sexadj_cfg(rosmap_cfg_plus_clinical)

## ----------------------------------------------------------------------------
## RUN
## Edit Synapse IDs above (or override a field right before calling, e.g.
##   msbb_cfg$input_synid <- "synNEWID"
##   msbb_cfg_sexadj <- make_sexadj_cfg(msbb_cfg, output_parent_synid = "synNEWPARENT")
## ) then run any of:
## ----------------------------------------------------------------------------

## Sex-stratified
result_msbb   <- run_DE_pipeline(msbb_cfg)
result_rosmap <- run_DE_pipeline(rosmap_cfg)
result_mayo   <- run_DE_pipeline(mayo_cfg)

result_msbb_plus_clinical   <- run_DE_pipeline(msbb_cfg_plus_clinical)
result_rosmap_plus_clinical <- run_DE_pipeline(rosmap_cfg_plus_clinical)

## Sex-adjusted (sex as covariate)
result_msbb_sexadj   <- run_DE_pipeline(msbb_cfg_sexadj)
result_rosmap_sexadj <- run_DE_pipeline(rosmap_cfg_sexadj)
result_mayo_sexadj   <- run_DE_pipeline(mayo_cfg_sexadj)

result_msbb_plus_clinical_sexadj   <- run_DE_pipeline(msbb_cfg_plus_clinical_sexadj)
result_rosmap_plus_clinical_sexadj <- run_DE_pipeline(rosmap_cfg_plus_clinical_sexadj)


