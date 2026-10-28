## ============================================================================
## UNIFIED DE-ON-RESIDUALS PIPELINE — MSBB / ROSMAP / MAYO
## Consolidates DE_MSBB_Residuals_for_sharing.R, DE_ROSMAP_Residuals_for_sharing.R
## and DE_MAYO_Residuals_for_sharing.R into one engine driven by per-cohort
## config lists.
## J Lundin
##
## Each original script did the same thing end to end, in two stages:
##   STAGE 1 (residuals): read a cqn'd dge_cqn/metadata object from Synapse ->
##     cohort-specific metadata prep -> set factors -> voomWithDreamWeights on a
##     technical-covariate formula -> swap in the cqn expression -> dream with
##     computeResiduals = TRUE -> residuals -> save <COHORT>_DE_res.rds + synStore.
##   STAGE 2 (DE on residuals): cell-means (~0 + group + covariates + (1|individualID))
##     formula on the residualized matrix -> contrasts (AD vs CT, by sex x tissue)
##     -> dream -> eBayes -> topTable per contrast -> save <COHORT>_DE_res2.rds
##     + synStore.
##
## What differed between cohorts (and is now just config, not code):
##   - input Synapse ID / output parent Synapse ID / output filenames
##   - tissue levels, diagnosis variable + levels, case/control labels
##   - the residual (technical) formula, incl. which variable is the batch
##   - how residuals were extracted: stats::residuals(fit) (ROSMAP, MSBB) vs
##     variancePartition::residuals(fit, vobj) (MAYO) — kept per cohort so
##     results reproduce; worth checking once and standardizing
##   - cohort-specific metadata prep (MAYO: diagnosis text -> diag2;
##     ROSMAP: factor final_batch + diagnosis; MSBB: factor sequencingBatch)
##
## To run a cohort: edit the Synapse IDs in the matching cfg block below (or
## override them right before calling run_DE_residuals_pipeline(), see bottom
## of file), then call run_DE_residuals_pipeline(msbb_cfg) /
## run_DE_residuals_pipeline(rosmap_cfg) / run_DE_residuals_pipeline(mayo_cfg).
##
## NOTE ON RESULTS NAMING: instead of separate variables per contrast
## (males_DLPFC6_res3, females_FP_res3, ...), each run returns DE_save$results,
## a named list keyed by contrast name (e.g. "AD2_vs_CT2_male.DLPFC"). Downstream
## code that expects the old list elements will need `results[["<name>"]]`
## instead.
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

rank_inverse_normal <- function(x, k = 3 / 8) {
  r <- rank(x, na.last = "keep", ties.method = "average")
  n <- sum(!is.na(x))
  qnorm((r - k) / (n - 2 * k + 1))
}


# Builds the AD_vs_CT-style contrasts (case - control) for every sex x tissue
# combination, in the same "tissue outer loop, sex inner loop (male, female)"
# order used in all three original scripts. `case_label`/`control_label` let
# the contrast NAME differ from the factor LEVEL if ever needed.
build_group_contrasts <- function(tissue_levels, sex_levels,
                                  case_level, control_level,
                                  case_label = case_level,
                                  control_label = control_level,
                                  sex_mode = c("stratified", "adjusted"))  {
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
    "name", "input_synid", "output_parent_synid_res", "output_parent_synid_de",
    "output_filename_res", "output_filename_de",
    "tissue_levels", "sex_levels", "diagnosis_var", "diagnosis_levels",
    "case_level", "control_level", "residual_formula_str", "de_formula_str",
    "residuals_method", "prep_metadata", "n_workers"
  )
  missing <- setdiff(required, names(cfg))
  if (length(missing) > 0) {
    stop("cfg for '", cfg$name %||% "?", "' is missing: ", paste(missing, collapse = ", "))
  }
  if (!cfg$residuals_method %in% c("stats", "variancePartition")) {
    stop("cfg$residuals_method must be 'stats' or 'variancePartition'")
  }
  sex_mode <- cfg$sex_mode %||% "stratified"
  if (!sex_mode %in% c("stratified", "adjusted")) {
    stop("cfg for '", cfg$name, "': sex_mode must be 'stratified' or 'adjusted'")
  }
  has_sex <- grepl("\\bsex\\b", cfg$de_formula_str)
  if (sex_mode == "adjusted" && !has_sex) {
    warning("cfg for '", cfg$name, "': sex_mode = 'adjusted' but formula_str has no `sex` covariate")
  }
  if (sex_mode == "stratified" && has_sex) {
    stop("cfg for '", cfg$name, "': sex_mode = 'stratified' but formula_str includes `sex` ",
         "(sex is already in `group`; the design would be rank-deficient)")
  }
}
`%||%` <- function(a, b) if (is.null(a)) b else a


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
  
  if (!grepl("\\bsex\\b", new$de_formula_str)) {
    new$de_formula_str <- sub("\\bgroup\\b", "group + sex", new$de_formula_str)
  }
  
  
rename_out <- function(f) {
    if (is.null(f)) return(f)
    if (grepl("sex_stratified", f)) return(sub("sex_stratified", "sex_adj", f))
    sub("(\\.[^.]+)$", "_sex_adj\\1", f)
  }
  new$output_filename_de <- rename_out(new$output_filename_de)
  
  if (keep_counts) {
    new$out_counts_filename <- rename_out(new$out_counts_filename)
  } else {
    new$out_counts_filename <- NULL
  }
  
  overrides <- list(...)
  for (nm in names(overrides)) new[[nm]] <- overrides[[nm]]
  new
}


# saveRDS + synStore, recording the input file as provenance
#save_and_store <- function(obj, filename, cfg) {
#  saveRDS(obj, file = filename)
#  invisible(synLogin(silent = TRUE))
#  synStore(File(path = filename, parent = cfg$output_parent_synid), used = cfg$input_synid)
#  message("Saved '", filename, "' -> ", cfg$output_parent_synid)
#}

## ----------------------------------------------------------------------------
## Shared engine — identical steps for every cohort
## ----------------------------------------------------------------------------

run_DE_residuals_pipeline <- function(cfg, qc_plots = TRUE) {
  validate_cfg(cfg)
  sex_mode <- cfg$sex_mode %||% "stratified"
  message("== Running DE pipeline: ", cfg$name, " (sex ", sex_mode, ") ==")
  
  # 1. Pull the cqn'd counts + metadata object from Synapse
  file_entity <- synGet(cfg$input_synid)
  cqn_obj <- readRDS(file_entity$path)
  md_sv   <- cqn_obj$metadata
  dge_cqn <- cqn_obj$dge_cqn
  
  # 2. Cohort-specific metadata prep (diagnosis derivation, batch factors)
  md_sv <- cfg$prep_metadata(md_sv, qc_plots = qc_plots)
 
  # formula
  form_final0 <- as.formula(cfg$formula_str)
  
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
  
  # 3. Align dge_cqn <-> md_sv exactly as all three originals did
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
  md_sv$sex     <- factor(md_sv$sex,     levels = cfg$sex_levels)
  md_sv[[cfg$diagnosis_var]] <- factor(md_sv[[cfg$diagnosis_var]], levels = cfg$diagnosis_levels)
  md_sv <- as.data.frame(md_sv)
  
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
  
  param_n <- MulticoreParam(workers = cfg$n_workers, RNGseed = 123, progressbar = TRUE)
  
  ## --------------------------------------------------------------------------
  ## STAGE 1: residualize technical covariates
  ## --------------------------------------------------------------------------
  
  res_final <- as.formula(cfg$residual_formula_str)
  environment(res_final) <- globalenv()
  
  # 5. voom weights, then swap in the cqn expression (keeps voom weights)
  voom_res <- voomWithDreamWeights(
    counts  = dge_cqn,
    formula = res_final,
    data    = md_sv,
    plot    = qc_plots,
    BPPARAM = param_n
  )
  #voom_res$E <- as.matrix(dge_cqn$E)
  voom_res$E <- as.matrix(dge_cqn$E)[names(voom_res$voom.xy$x) %||% rownames(dge_cqn$E),
                                                 colnames(voom_res), drop = FALSE]  
  # 6. dream on technical covariates + residuals
  fit_res <- variancePartition::dream(
    exprObj          = voom_res,
    formula          = res_final,
    data             = md_sv,
    computeResiduals = TRUE,
    BPPARAM          = param_n
  )
  
  res_counts <- switch(cfg$residuals_method,
                       stats             = stats::residuals(fit_res),
                       variancePartition = variancePartition::residuals(fit_res, voom_res)
  )
  if (qc_plots) print(res_counts[1:5, 1:5])
  
  # Stack into one long table + add gene info
  res_counts2 <- res_counts |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "ensembl_gene_id")
  
  #res_counts2 <- res_counts |>
   # lapply(rownames_to_column, var = "ensembl_gene_id") |>
  #  bind_rows(.id = "Comparison")
  
  syn_gene_file <- synapser::synGet("syn74907403") ## created using code from here: https://github.com/jaclynbeck-sage/sageRNAUtils/blob/main/R/gc_content_functions.R
  gene_info <- read.csv(syn_gene_file$path)
  combined <- merge(res_counts2, gene_info, by = "ensembl_gene_id", all.x=T)
  
  # 8. Save DE results CSV and upload
  write.csv(combined, cfg$output_filename_res, row.names = FALSE)
  invisible(synLogin(silent = TRUE))
  synStore(File(cfg$output_filename_res, parent = cfg$output_parent_synid_res))
  message("Saved '", cfg$output_filename_res, "' -> ", cfg$output_parent_synid_res)
  
  
 # DE_res <- list(
 #   metadata   = md_sv,
 #   vobj_res   = voom_res,
  #  dge_cqn    = dge_cqn,
 #   fit_res    = fit_res,
 #   res_counts = res_counts
 # )
  
  # 7. Save + upload residuals
 # save_and_store(DE_res, cfg$output_filename_res, cfg)
  
  ## --------------------------------------------------------------------------
  ## STAGE 2: DE on the residualized matrix
  ## --------------------------------------------------------------------------
  
  #md_sv$group <- interaction(md_sv[[cfg$diagnosis_var]], md_sv$sex, md_sv$tissue2, drop = TRUE)
  group <- md_sv$group
  # some formula terms are resolved via the global env as a fallback, matching
  # the original scripts' `environment(form) <- globalenv()` pattern
  assign("group", group, envir = globalenv())
  
  if (qc_plots) print(table(md_sv$group))
  
  form_ck <- as.formula(cfg$de_formula_str)
  environment(form_ck) <- globalenv()
  
  contrasts_vec <- build_group_contrasts(
    tissue_levels  = cfg$tissue_levels,
    sex_levels     = cfg$sex_levels,
    case_level     = cfg$case_level,
    control_level  = cfg$control_level,
    case_label     = cfg$case_label %||% cfg$case_level,
    control_label  = cfg$control_label %||% cfg$control_level,
    sex_mode       = sex_mode
  )
  
  L_group <- makeContrastsDream(form_ck, md_sv, contrasts = contrasts_vec)
  if (qc_plots) print(plotContrasts(L_group))
  
  # 8. dream on residuals (res = residualized matrix from stage 1)
  fit_res_contrasts <- variancePartition::dream(
    exprObj = res_counts,
    formula = form_ck,
    data    = md_sv,
    L       = L_group,
    BPPARAM = param_n
  )
  ebayes_res <- variancePartition::eBayes(fit_res_contrasts)
  
  # 9. One topTable per contrast, collected into a named list
  #    (selected by contrast name, so coefficient order can't matter)
  de_results <- setNames(
    lapply(colnames(L_group), function(nm) topTable(ebayes_res, coef = nm, number = Inf, confint = TRUE)),
    colnames(L_group)
  )
  
  # Sex main effect (adjusted mode only): female vs male, adjusted for
  # group + covariates
  results_sex <- list()
  if (sex_mode == "adjusted") {
    sex_coef <- grep("^sex", colnames(coef(ebayes_res)), value = TRUE)
    results_sex <- setNames(
      lapply(sex_coef, function(nm) topTable(ebayes_res, coef = nm, number = Inf, confint = TRUE)),
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
  write.csv(combined, cfg$output_filename_de, row.names = FALSE)
  invisible(synLogin(silent = TRUE))
  synStore(File(cfg$output_filename_de, parent = cfg$output_parent_synid_de))
  message("Saved '", cfg$output_filename_de, "' -> ", cfg$output_parent_synid_de)
  
  
 # DE_save <- list(
#    metadata_res      = md_sv,
 #   fit_res_contrasts = fit_res_contrasts,
 #   ebayes_res        = ebayes_res,
 #   results           = de_results
#  )
  
  # 10. Save + upload DE results
 # save_and_store(DE_save, cfg$output_filename_de, cfg)
  
 # invisible(list(DE_res = DE_res, DE_save = DE_save))
}

## ----------------------------------------------------------------------------
## Cohort-specific metadata prep
## ----------------------------------------------------------------------------

msbb_prep_metadata <- function(md_sv, qc_plots = TRUE) {
  md_sv$ageDeath_invnorm <- rank_inverse_normal(md_sv$ageDeath)
  md_sv <- md_sv[md_sv$tissue2 != "PC", ]
  md_sv[md_sv$tissue2 == "PG", "tissue2"] <- "PHG"
  md_sv$PMI_log <- log(md_sv$PMI)
  md_sv$sequencingBatch <- factor(md_sv$sequencingBatch)
  if (qc_plots) {
    print(table(md_sv$diagnosis, useNA = "always"))
    print(table(md_sv$diagnosis, md_sv$sex, useNA = "always"))
  }
  md_sv
}

rosmap_prep_metadata <- function(md_sv, qc_plots = TRUE) {
  md_sv$ageDeath_invnorm <- rank_inverse_normal(md_sv$ageDeath)
  md_sv[md_sv$tissue2 == "CN", "tissue2"] <- "ACC"
  md_sv$PMI_log <- log(md_sv$PMI)
  md_sv$final_batch <- factor(md_sv$final_batch)
  md_sv$diagnosis   <- factor(md_sv$diagnosis, levels = c("CT", "AD", "OTHER"))
  if (qc_plots) print(table(md_sv$diag2, useNA = "always"))
  md_sv
}

mayo_prep_metadata <- function(md_sv, qc_plots = TRUE) {
  md_sv$ageDeath_invnorm <- rank_inverse_normal(md_sv$ageDeath)
  md_sv[md_sv$tissue2 == "CER", "tissue2"] <- "CBE"
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

de_formula_shared <- "~ 0 + group + apoe4Status + ageDeath_invnorm + (1|individualID)"

msbb_cfg <- list(
  name                 = "MSBB",
  sex_mode                = "stratified",
  input_synid          = "syn77540009",  # MSBB_md_counts_cqn_FINAL.rds - new syn folder
  #output_parent_synid  = "syn77539791", # Staging
  output_parent_synid_res  = "syn77802600", # res-counts
  output_parent_synid_de  = "syn77956686", # de-res-counts
  output_filename_res  = "res-counts_msbb.csv",
  output_filename_de   = "de-res_MSBB_neuro.path.diag_sex_stratified.csv",
  tissue_levels        = c("FP", "IFG", "PHG", "STG"),
  sex_levels           = c("male", "female"),
  diagnosis_var        = "diag2",
  diagnosis_levels     = c("CT2", "AD2", "OTHER2"),
  case_level           = "AD2",
  control_level        = "CT2",
  case_label              = "AD",     # contrast names use AD/CT, not AD2/CT2
  control_label           = "CT",
  residual_formula_str = "~ PMI_log + RIN + PC1_metrics + PC2_metrics + PC3_metrics + PC4_metrics + (1|sequencingBatch)",
  residuals_method     = "variancePartition",
  de_formula_str       = de_formula_shared,
  prep_metadata        = msbb_prep_metadata,
  n_workers            = 4
)

rosmap_cfg <- list(
  name                 = "ROSMAP",
  sex_mode                = "stratified",
  input_synid         = "syn78387829",  # ROSMAP_md_counts_cqn_DLPFC_CN_PCC_FINAL_updated_filtered.rds - single CQN run
  #output_parent_synid = "syn77539791", # Staging
  output_parent_synid_res  = "syn77802600", # res-counts
  output_parent_synid_de  = "syn77956686", # de-res-counts
  output_filename_res  = "res-counts_rosmap.csv",  
  output_filename_de   = "de-res_ROSMAP_neuro.path.diag_sex_stratified.csv",
  tissue_levels        = c("DLPFC", "PCC", "ACC"),
  sex_levels           = c("male", "female"),
  diagnosis_var        = "diag2",
  diagnosis_levels     = c("CT2", "AD2", "OTHER2"),
  case_level           = "AD2",
  control_level        = "CT2",
  case_label              = "AD",     # contrast names use AD/CT, not AD2/CT2
  control_label           = "CT",
  residual_formula_str = "~ PMI_log + RIN + libraryPrep + PC1_metrics + PC2_metrics + PC3_metrics + PC4_metrics + (1|final_batch)",
  residuals_method     = "variancePartition",
  de_formula_str       = de_formula_shared,
  prep_metadata        = rosmap_prep_metadata,
  n_workers            = 4
)

mayo_cfg <- list(
  name                 = "MAYO",
  sex_mode                = "stratified",
  input_synid          = "syn77548692",  # MAYO_md_counts_cqn_FINAL.rds - new syn folder
  output_parent_synid_res  = "syn77802600", # res-counts
  output_parent_synid_de  = "syn77956686", # de-res-counts
  output_filename_res  = "res-counts_mayo.csv",
  output_filename_de   = "de-res_MAYO_neuro.path.diag_sex_stratified.csv",
  tissue_levels        = c("CBE", "TCX"),
  sex_levels           = c("male", "female"),
  diagnosis_var        = "diag2",
  diagnosis_levels     = c("CT2", "AD2", "PathAg2", "PSP2"),
  case_level           = "AD2",
  control_level        = "CT2",
  case_label              = "AD",     # contrast names use AD/CT, not AD2/CT2
  control_label           = "CT",
  residual_formula_str = "~ RIN + PC1_metrics + PC2_metrics + PC3_metrics + (1|flowcell)",
  residuals_method     = "variancePartition",
  de_formula_str       = de_formula_shared,
  prep_metadata        = mayo_prep_metadata,
  n_workers            = 4
)

## ----------------------------------------------------------------------------
## Sex-adjusted configs — derived from the stratified ones above
## (formula gets "+ sex", outputs renamed "..._sex_adj.csv", counts file
## skipped since it would be identical to the stratified run's)
## ----------------------------------------------------------------------------

msbb_cfg_sexadj                 <- make_sexadj_cfg(msbb_cfg)
rosmap_cfg_sexadj               <- make_sexadj_cfg(rosmap_cfg)
mayo_cfg_sexadj                 <- make_sexadj_cfg(mayo_cfg)
#msbb_cfg_plus_clinical_sexadj   <- make_sexadj_cfg(msbb_cfg_plus_clinical)
#rosmap_cfg_plus_clinical_sexadj <- make_sexadj_cfg(rosmap_cfg_plus_clinical)

## ----------------------------------------------------------------------------
## RUN
## Edit Synapse IDs above (or override a field right before calling, e.g.
##   msbb_cfg$input_synid <- "synNEWID"
##   msbb_cfg$output_parent_synid <- "synNEWPARENT"
## ) then run one of:
## ----------------------------------------------------------------------------

## Sex-stratified
result_msbb   <- run_DE_residuals_pipeline(msbb_cfg)
result_rosmap <- run_DE_residuals_pipeline(rosmap_cfg)
result_mayo   <- run_DE_residuals_pipeline(mayo_cfg)

## Sex-adjusted (sex as covariate)
result_msbb_sexadj   <- run_DE_residuals_pipeline(msbb_cfg_sexadj)
result_rosmap_sexadj <- run_DE_residuals_pipeline(rosmap_cfg_sexadj)
result_mayo_sexadj   <- run_DE_residuals_pipeline(mayo_cfg_sexadj)



