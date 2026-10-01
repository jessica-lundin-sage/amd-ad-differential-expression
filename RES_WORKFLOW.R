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
                                  control_label = control_level) {
  combos <- expand.grid(sex = sex_levels, tissue = tissue_levels, stringsAsFactors = FALSE)
  combos <- combos[order(match(combos$tissue, tissue_levels), match(combos$sex, sex_levels)), ]
  
  exprs <- sprintf(
    "(group%s.%s.%s - group%s.%s.%s)",
    case_level, combos$sex, combos$tissue,
    control_level, combos$sex, combos$tissue
  )
  nms <- sprintf(
    "%s_vs_%s_%s.%s",
    case_label, control_label, combos$sex, combos$tissue
  )
  setNames(exprs, nms)
}

# Minimal defensive check that a cfg has everything the engine needs.
validate_cfg <- function(cfg) {
  required <- c(
    "name", "input_synid", "output_parent_synid",
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
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# saveRDS + synStore, recording the input file as provenance
save_and_store <- function(obj, filename, cfg) {
  saveRDS(obj, file = filename)
  invisible(synLogin(silent = TRUE))
  synStore(File(path = filename, parent = cfg$output_parent_synid), used = cfg$input_synid)
  message("Saved '", filename, "' -> ", cfg$output_parent_synid)
}

## ----------------------------------------------------------------------------
## Shared engine — identical steps for every cohort
## ----------------------------------------------------------------------------

run_DE_residuals_pipeline <- function(cfg, qc_plots = TRUE) {
  validate_cfg(cfg)
  message("== Running DE-on-residuals pipeline: ", cfg$name, " ==")
  
  # 1. Pull the cqn'd counts + metadata object from Synapse
  file_entity <- synGet(cfg$input_synid)
  cqn_obj <- readRDS(file_entity$path)
  md_sv   <- cqn_obj$metadata
  dge_cqn <- cqn_obj$dge_cqn
  
  # 2. Cohort-specific metadata prep (diagnosis derivation, batch factors)
  md_sv <- cfg$prep_metadata(md_sv, qc_plots = qc_plots)
  
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
  
  param_n <- MulticoreParam(workers = cfg$n_workers, progressbar = TRUE)
  
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
  
  DE_res <- list(
    metadata   = md_sv,
    vobj_res   = voom_res,
    dge_cqn    = dge_cqn,
    fit_res    = fit_res,
    res_counts = res_counts
  )
  
  # 7. Save + upload residuals
  save_and_store(DE_res, cfg$output_filename_res, cfg)
  
  ## --------------------------------------------------------------------------
  ## STAGE 2: DE on the residualized matrix
  ## --------------------------------------------------------------------------
  
  md_sv$group <- interaction(md_sv[[cfg$diagnosis_var]], md_sv$sex, md_sv$tissue2, drop = TRUE)
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
    control_label  = cfg$control_label %||% cfg$control_level
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
  
  if (qc_plots) {
    for (nm in names(de_results)) {
      tt <- de_results[[nm]]
      plot(tt$logFC, -log10(tt$adj.P.Val), main = nm, xlab = "logFC", ylab = "-log10(adj.P.Val)")
    }
  }
  
  DE_save <- list(
    metadata_res      = md_sv,
    fit_res_contrasts = fit_res_contrasts,
    ebayes_res        = ebayes_res,
    results           = de_results
  )
  
  # 10. Save + upload DE results
  save_and_store(DE_save, cfg$output_filename_de, cfg)
  
  invisible(list(DE_res = DE_res, DE_save = DE_save))
}

## ----------------------------------------------------------------------------
## Cohort-specific metadata prep
## ----------------------------------------------------------------------------

msbb_prep_metadata <- function(md_sv, qc_plots = TRUE) {
  md_sv$ageDeath_invnorm <- rank_inverse_normal(md_sv$ageDeath)
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
  md_sv$PMI_log <- log(md_sv$PMI)
  md_sv$final_batch <- factor(md_sv$final_batch)
  md_sv$diagnosis   <- factor(md_sv$diagnosis, levels = c("CT", "AD", "OTHER"))
  if (qc_plots) print(table(md_sv$diag2, useNA = "always"))
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

de_formula_shared <- "~ 0 + group + apoe4Status + ageDeath_invnorm + (1|individualID)"

msbb_cfg <- list(
  name                 = "MSBB",
  input_synid          = "syn77540009",  # MSBB_md_counts_cqn_FINAL.rds - new syn folder
  #output_parent_synid  = "syn77539791", # Staging
  output_parent_synid  = "syn77615867", # amp-ad-rnaseq_reprocessing_intermediate_files
  output_filename_res  = "msbb_md_res-counts_final.rds",
  output_filename_de   = "amp-ad-de-res_MSBB_neuro.path.diag_sex_stratified.rds",
  tissue_levels        = c("FP", "IFG", "PG", "STG"),
  sex_levels           = c("male", "female"),
  diagnosis_var        = "diag2",
  diagnosis_levels     = c("CT2", "AD2", "OTHER2"),
  case_level           = "AD2",
  control_level        = "CT2",
  residual_formula_str = "~ PMI_log + RIN + PC1_metrics + PC2_metrics + PC3_metrics + PC4_metrics + (1|sequencingBatch)",
  residuals_method     = "variancePartition",
  de_formula_str       = de_formula_shared,
  prep_metadata        = msbb_prep_metadata,
  n_workers            = 4
)

rosmap_cfg <- list(
  name                 = "ROSMAP",
  input_synid         = "syn77559937",  # ROSMAP_md_counts_cqn_DLPFC_CN_PCC_FINAL.rds - new syn folder
  #output_parent_synid = "syn77539791", # Staging
  output_parent_synid = "syn77615867", # amp-ad-rnaseq_reprocessing_intermediate_files
  output_filename_res  = "rosmap_md_res-counts_final.rds",
  output_filename_de   = "amp-ad-de-res_ROSMAP_neuro.path.diag_sex_stratified.rds",
  tissue_levels        = c("DLPFC", "PCC", "CN"),
  sex_levels           = c("male", "female"),
  diagnosis_var        = "diag2",
  diagnosis_levels     = c("CT2", "AD2", "OTHER2"),
  case_level           = "AD2",
  control_level        = "CT2",
  residual_formula_str = "~ PMI_log + RIN + PC1_metrics + PC2_metrics + PC3_metrics + PC4_metrics + (1|final_batch)",
  residuals_method     = "variancePartition",
  de_formula_str       = de_formula_shared,
  prep_metadata        = rosmap_prep_metadata,
  n_workers            = 4
)

mayo_cfg <- list(
  name                 = "MAYO",
  input_synid          = "syn77548692",  # MAYO_md_counts_cqn_FINAL.rds - new syn folder
  output_parent_synid  = "syn77615867", # amp-ad-rnaseq_reprocessing_intermediate_files
  output_filename_res  = "mayo_md_res-counts_final.rds",
  output_filename_de   = "amp-ad-de-res_MAYO_neuro.path.diag_sex_stratified.rds",
  tissue_levels        = c("CER", "TCX"),
  sex_levels           = c("male", "female"),
  diagnosis_var        = "diag2",
  diagnosis_levels     = c("CT2", "AD2", "PathAg2", "PSP2"),
  case_level           = "AD2",
  control_level        = "CT2",
  residual_formula_str = "~ RIN + PC1_metrics + PC2_metrics + PC3_metrics + (1|flowcell)",
  residuals_method     = "variancePartition",
  de_formula_str       = de_formula_shared,
  prep_metadata        = mayo_prep_metadata,
  n_workers            = 4
)

## ----------------------------------------------------------------------------
## RUN
## Edit Synapse IDs above (or override a field right before calling, e.g.
##   msbb_cfg$input_synid <- "synNEWID"
##   msbb_cfg$output_parent_synid <- "synNEWPARENT"
## ) then run one of:
## ----------------------------------------------------------------------------

result_msbb   <- run_DE_residuals_pipeline(msbb_cfg)
result_rosmap <- run_DE_residuals_pipeline(rosmap_cfg)
result_mayo   <- run_DE_residuals_pipeline(mayo_cfg)
