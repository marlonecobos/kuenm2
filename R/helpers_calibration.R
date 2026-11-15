eval_stats <- function(cal_res, error_considered, algorithm) {
  if (any(missing(cal_res), missing(error_considered), missing(algorithm))) {
    stop("Arguments 'cal_res', 'error_considered', 'algorithm' must be defined.")
  }
  if (!inherits(cal_res, "data.frame")) {
    stop("Argument 'cal_res' must be a 'dataframe'.")
  }
  if (!inherits(error_considered, "numeric")) {
    stop("Argument 'error_considered' must be 'numeric'.")
  }
  if (!inherits(algorithm, "character") || !(algorithm %in% c("glm", "maxnet"))) {
    stop("Argument algorithm must be 'glm' or 'maxnet'")
  }

  omission_rates <- paste0("Omission_rate_at_", error_considered)
  proc_values <- paste0("Mean_AUC_ratio_at_", error_considered)
  pval_values <- paste0("pval_pROC_at_", error_considered)
  toagg <- c(omission_rates, proc_values, pval_values)

  meta_cols <- if (algorithm == "maxnet") {
    c("ID", "Formulas", "R_multiplier", "Features", "AICc",
      "Parameters", "Is_concave")
  } else {
    c("ID", "Formulas", "Features", "AIC", "Parameters", "Is_concave")
  }

  # Fast path when cal_res holds folds for a single candidate model
  if (length(unique(cal_res$ID)) == 1L) {
    n_k <- nrow(cal_res)
    means <- round(colMeans(cal_res[, toagg, drop = FALSE]), 4L)
    sds <- if (n_k > 1L) {
      round(vapply(toagg, function(cn) stats::sd(cal_res[[cn]]), numeric(1L)), 4L)
    } else {
      rep(NA_real_, length(toagg))
    }

    stat_vals <- unlist(lapply(seq_along(toagg), function(idx) {
      v <- c(means[idx], sds[idx])
      names(v) <- c(paste0(toagg[idx], ".mean"), paste0(toagg[idx], ".sd"))
      v
    }))

    return(cbind(cal_res[1L, meta_cols, drop = FALSE],
                 as.data.frame(as.list(stat_vals)),
                 row.names = NULL))
  }

  # General path for multi-model data.frames
  agg_by <- if (algorithm == "maxnet") {
    c("Formulas", "R_multiplier", "Features")
  } else {
    c("Formulas", "Features")
  }
  agg_formula <- paste("~", paste(agg_by, collapse = " + "))

  xy <- lapply(toagg, function(x) {
    do.call(
      data.frame,
      stats::aggregate(stats::as.formula(paste(x, agg_formula)),
                       data = cal_res, FUN = function(y) {
                         c(mean = round(mean(y), 4L), sd = round(stats::sd(y), 4L))
                       }, na.action = NULL)
    )
  })

  stats <- Reduce(function(x, y) merge(x, y, by = agg_by), xy)
  stats_AICS <- cal_res[!duplicated(cal_res[, meta_cols]), meta_cols, drop = FALSE]
  merge(stats, stats_AICS, by = agg_by)
}


empty_replicates <- function(error_considered, n_row, replicates,
                             is_c, algorithm) {
  if (any(missing(error_considered), missing(n_row), missing(replicates),
          missing(is_c), missing(algorithm))) {
    stop("Arguments 'error_considered', 'n_row', 'replicates', 'is_c', and 'algorithm' must be defined.")
  }
  if (!inherits(error_considered, "numeric")) {
    stop("Argument 'error_considered' must be 'numeric'.")
  }
  if (!is.numeric(n_row)) {
    stop("Argument 'n_row' must be 'numeric'.")
  }
  if (all(!inherits(replicates, "numeric") & !inherits(replicates, "character"))) {
    stop("Argument 'replicates' must be 'numeric' or 'character'.")
  }
  if (!is.na(is_c) && !inherits(is_c, "logical")) {
    stop("Argument 'is_c' must be NA or 'logical'.")
  }
  if (!inherits(algorithm, "character") || !(algorithm %in% c("glm", "maxnet"))) {
    stop("Argument algorithm must be 'glm' or 'maxnet'")
  }

  aic_col <- if (algorithm == "maxnet") "AICc" else "AIC"
  column_names <- c("Partition",
                    paste0("Omission_rate_at_", error_considered),
                    paste0("Mean_AUC_ratio_at_", error_considered),
                    paste0("pval_pROC_at_", error_considered),
                    aic_col, "Parameters", "Is_concave")

  df_eval_q <- data.frame(matrix(NA, nrow = n_row, ncol = length(column_names)))
  colnames(df_eval_q) <- column_names
  df_eval_q$Partition <- replicates
  df_eval_q$Is_concave <- is_c

  return(df_eval_q)
}


empty_summary <- function(error_considered, is_c, algorithm) {
  if (any(missing(error_considered), missing(is_c), missing(algorithm))) {
    stop("Arguments 'error_considered', 'is_c', and 'algorithm' must be defined.")
  }
  if (!inherits(error_considered, "numeric")) {
    stop("Argument 'error_considered' must be 'numeric'.")
  }
  if (!is.na(is_c) && !inherits(is_c, "logical")) {
    stop("Argument 'is_c' must be NA or 'logical'.")
  }
  if (!inherits(algorithm, "character") || !(algorithm %in% c("glm", "maxnet"))) {
    stop("Argument algorithm must be 'glm' or 'maxnet'")
  }

  om_means <- paste0("Omission_rate_at_", error_considered, ".mean")
  om_sd <- paste0("Omission_rate_at_", error_considered, ".sd")
  auc_means <- paste0("Mean_AUC_ratio_at_", error_considered, ".mean")
  auc_sd <- paste0("Mean_AUC_ratio_at_", error_considered, ".sd")
  pval_means <- paste0("pval_pROC_at_", error_considered, ".mean")
  pval_sd <- paste0("pval_pROC_at_", error_considered, ".sd")

  aic_col <- if (algorithm == "maxnet") "AICc" else "AIC"
  column_names <- c(om_means, om_sd,
                    auc_means, auc_sd,
                    pval_means, pval_sd,
                    aic_col, "Parameters", "Is_concave")

  eval_final_q <- data.frame(matrix(NA, nrow = 1, ncol = length(column_names)))
  colnames(eval_final_q) <- column_names
  eval_final_q$Is_concave <- is_c

  return(eval_final_q)
}


# Precompute design matrix and base regularization once per formula
prep_design <- function(p, data, f, addsamplestobackground = TRUE,
                        weights = NULL) {
  iniweigth <- is.null(weights)
  if (iniweigth) {
    weights <- ifelse(p == 1, 1, 100)
  }

  if (addsamplestobackground) {
    is_p <- (p == 1)
    pdata <- data[is_p, , drop = FALSE]
    ndata <- data[!is_p, , drop = FALSE]

    wadd <- !(do.call(paste, pdata) %in% do.call(paste, ndata))
    if (any(wadd)) {
      n_add <- sum(wadd)
      p <- c(p, rep(0, n_add))
      data <- rbind(data, pdata[wadd, , drop = FALSE])

      if (iniweigth) {
        weights <- c(weights, rep(100, n_add))
      } else {
        pweight <- weights[is_p]
        weights <- c(weights, pweight[wadd])
      }
    }
  }

  mm <- stats::model.matrix(f, data)
  if (ncol(mm) < 2L) {
    mm <- cbind(0, mm)
  }
  reg <- maxnet.default.regularization(p, mm)
  y_fac <- structure(as.integer(p) + 1L, levels = c("0", "1"), class = "factor")

  list(p = p, y_fac = y_fac, data = data, weights = weights, mm = mm, reg = reg)
}


# Direct matrix prediction on precomputed model.matrix
predict_mm <- function(mod, mm_new, type = "cloglog") {
  betas <- mod$betas
  if (length(betas) > 0L) {
    link <- as.numeric(mm_new[, names(betas), drop = FALSE] %*% betas) + mod$alpha
  } else {
    link <- rep(mod$alpha, nrow(mm_new))
  }
  if (type == "exponential") return(exp(link))
  if (type == "cloglog") return(1 - exp(-exp(link)))
  stop("Unsupported type.")
}


# Grouped maxnet evaluation: builds design matrices once per formula and fold,
# reusing them across all R_multiplier values
fit_eval_formula_group <- function(group_df, data, error_considered,
                                   omission_rate, write_summary,
                                   addsamplestobackground, weights = NULL,
                                   return_all_results, proc_for_all,
                                   out_dir) {
  formula_x <- stats::as.formula(group_df$Formulas[1L])
  n_rm <- nrow(group_df)
  n_folds <- length(data$part_data)
  bgind <- which(data$calibration_data$pr_bg == 0)
  occ_rows <- which(data$calibration_data$pr_bg == 1)

  # 1. Design matrix on calibration_data for fast prediction across folds
  mm_cal <- stats::model.matrix(formula_x, data$calibration_data)
  if (ncol(mm_cal) < 2L) mm_cal <- cbind(0, mm_cal)

  # 2. Precompute full-data design
  full_pre <- try(
    prep_design(p = data$calibration_data$pr_bg,
                data = data$calibration_data,
                f = formula_x,
                addsamplestobackground = addsamplestobackground,
                weights = weights),
    silent = TRUE
  )

  # 3. Precompute per-fold designs (dropping heavy data copy since only_fit = TRUE)
  if (!inherits(full_pre, "try-error")) {
    fold_pres <- lapply(seq_len(n_folds), function(i) {
      notrain <- -data$part_data[[i]]
      data_i <- data$calibration_data[notrain, , drop = FALSE]
      weights_i <- if (!is.null(weights)) weights[notrain] else NULL
      p_i <- try(
        prep_design(p = data_i$pr_bg, data = data_i, f = formula_x,
                    addsamplestobackground = addsamplestobackground,
                    weights = weights_i),
        silent = TRUE
      )
      if (inherits(p_i, "try-error")) return(NULL)
      p_i$data <- NULL
      p_i$notrain <- notrain
      p_i
    })
  } else {
    fold_pres <- vector("list", n_folds)
  }

  # 4. Fit and evaluate each R_multiplier reusing precomputed designs
  group_res <- lapply(seq_len(n_rm), function(r_idx) {
    grid_x <- group_df[r_idx, , drop = FALSE]
    reg_x <- grid_x$R_multiplier

    if (inherits(full_pre, "try-error")) {
      m_aic <- full_pre
    } else {
      m_aic <- try(
        glmnet_mx(f = formula_x, regmult = reg_x,
                  addsamplestobackground = addsamplestobackground,
                  only_fit = FALSE, precomp = full_pre),
        silent = TRUE
      )
    }

    if (inherits(m_aic, "try-error")) {
      is_c <- NA
      eval_final <- cbind(
        grid_x,
        empty_replicates(error_considered = error_considered,
                         n_row = n_folds,
                         replicates = names(data$part_data),
                         is_c = is_c, algorithm = "maxnet")
      )
      eval_final_summary <- reorder_stats_columns(
        cbind(grid_x, empty_summary(error_considered, is_c, "maxnet")),
        error_considered = error_considered
      )
    } else {
      npar <- length(m_aic$betas)
      vals <- predict_mm(m_aic, mm_cal[occ_rows, , drop = FALSE],
                         type = "exponential")
      AICc <- aic_ws(pred_occs = vals, ncoefs = npar)

      q_betas <- m_aic$betas[grepl("\\^2", names(m_aic$betas))]
      is_c <- if (length(q_betas) == 0L) FALSE else any(q_betas > 0)

      mods <- lapply(seq_len(n_folds), function(i) {
        fp_i <- fold_pres[[i]]
        if (is.null(fp_i)) return(NULL)

        mod_i <- try(
          glmnet_mx(f = formula_x, regmult = reg_x,
                    addsamplestobackground = addsamplestobackground,
                    only_fit = TRUE, precomp = fp_i),
          silent = TRUE
        )
        if (inherits(mod_i, "try-error")) return(NULL)

        pred_i <- predict_mm(mod_i, mm_cal, type = "cloglog")
        notrain <- fp_i$notrain
        suit_val_cal <- pred_i[unique(c(notrain, -bgind))]
        suit_val_eval <- pred_i[which(!-notrain %in% bgind)]
        om_rate <- omrat(threshold = error_considered, pred_train = suit_val_cal,
                         pred_test = suit_val_eval)

        if (proc_for_all) {
          proc_i <- unlist(lapply(error_considered, function(omr) {
            proc_omr <- fpROC::auc_metrics(test_prediction = suit_val_eval,
                                           prediction = pred_i,
                                           threshold = omr)$summary[, 4:5]
            names(proc_omr) <- c(paste0("Mean_AUC_ratio_at_", omr),
                                 paste0("pval_pROC_at_", omr))
            proc_omr
          }))
        } else {
          proc_i <- rep(NA, length(error_considered) * 2L)
          names(proc_i) <- c(paste0("Mean_AUC_ratio_at_", error_considered),
                             paste0("pval_pROC_at_", error_considered))
        }

        df_eval <- data.frame(Partition = i, t(om_rate), t(proc_i), AICc = AICc,
                              Parameters = npar, Is_concave = is_c,
                              row.names = NULL)
        cbind(grid_x, df_eval)
      })

      if (any(sapply(mods, is.null))) {
        is_c <- NA
        eval_final <- cbind(
          grid_x,
          empty_replicates(error_considered = error_considered,
                           n_row = n_folds,
                           replicates = names(data$part_data),
                           is_c = is_c, algorithm = "maxnet")
        )
        eval_final_summary <- reorder_stats_columns(
          cbind(grid_x, empty_summary(error_considered, is_c, "maxnet")),
          error_considered = error_considered
        )
      } else {
        names(mods) <- names(data$part_data)
        eval_final <- do.call("rbind", mods)
        eval_final_summary <- reorder_stats_columns(
          eval_stats(eval_final, error_considered, "maxnet"),
          error_considered = error_considered
        )
      }
    }

    if (write_summary) {
      utils::write.csv(eval_final_summary,
                       file.path(out_dir, paste0("Summary_cand_model_", grid_x$ID, ".csv")),
                       row.names = FALSE)
    }

    if (!return_all_results) eval_final <- NULL
    list(All_results = eval_final, Summary = eval_final_summary)
  })

  list(
    All_results = do.call("rbind", lapply(group_res, function(x) x$All_results)),
    Summary = do.call("rbind", lapply(group_res, function(x) x$Summary))
  )
}


fit_eval_concave <- function(x, q_grids, data, formula_grid, error_considered,
                             omission_rate, write_summary, addsamplestobackground,
                             weights = NULL, return_all_results, algorithm,
                             proc_for_all, out_dir) {
  grid_x <- q_grids[x, , drop = FALSE]

  if (algorithm == "maxnet") {
    formula_x <- stats::as.formula(grid_x$Formulas)
    reg_x <- grid_x$R_multiplier
    m_aic <- try(glmnet_mx(p = data$calibration_data$pr_bg,
                           data = data$calibration_data,
                           f = formula_x, regmult = reg_x,
                           addsamplestobackground = addsamplestobackground,
                           weights = weights),
                 silent = TRUE)
  } else if (algorithm == "glm") {
    formula_x <- stats::as.formula(paste("pr_bg", grid_x$Formulas))
    m_aic <- try(glm_mx(formula = formula_x, family = binomial(link = "cloglog"),
                        data = data$calibration_data, weights = weights),
                 silent = TRUE)
  } else {
    stop("Unsupported model type. Please use 'maxnet' or 'glm'.")
  }

  if (inherits(m_aic, "try-error")) {
    npar <- NA
    AICc <- NA
    is_c <- NA
  } else {
    npar <- if (algorithm == "maxnet") {
      length(m_aic$betas)
    } else {
      length(m_aic$coefficients[-1])
    }

    if (algorithm == "maxnet") {
      mm_cal <- stats::model.matrix(formula_x, data$calibration_data)
      if (ncol(mm_cal) < 2L) mm_cal <- cbind(0, mm_cal)
      occ_rows <- which(data$calibration_data$pr_bg == 1)
      vals <- predict_mm(m_aic, mm_cal[occ_rows, , drop = FALSE],
                         type = "exponential")
      AICc <- aic_ws(pred_occs = vals, ncoefs = npar)
    }

    q_betas <- if (algorithm == "maxnet") {
      m_aic$betas[grepl("\\^2", names(m_aic$betas))]
    } else {
      m_aic$coefficients[-1][grepl("\\^2", names(m_aic$coefficients[-1]))]
    }

    is_c <- if (length(q_betas) == 0L) FALSE else any(q_betas > 0)
  }

  if (isTRUE(is_c) || is.na(is_c)) {
    grid_q <- if (algorithm == "maxnet") {
      all_reg <- unique(formula_grid$R_multiplier)
      do.call("rbind", lapply(seq_along(all_reg), function(k) {
        grid_x_i <- grid_x
        grid_x_i$R_multiplier <- all_reg[k]
        grid_x_i$ID <- formula_grid[formula_grid$Formulas == grid_x$Formulas &
                                      formula_grid$R_multiplier == all_reg[k], "ID"]
        return(grid_x_i)
      }))
    } else {
      grid_x
    }

    n_folds <- length(data$part_data)
    rep_blocks <- lapply(seq_len(nrow(grid_q)), function(r_i) {
      cbind(
        grid_q[rep(r_i, n_folds), , drop = FALSE],
        empty_replicates(error_considered = error_considered,
                         n_row = n_folds,
                         replicates = names(data$part_data),
                         is_c = is_c,
                         algorithm = algorithm)
      )
    })
    eval_final_q <- reorder_stats_columns(do.call("rbind", rep_blocks),
                                          error_considered)
    row.names(eval_final_q) <- NULL

    sum_blocks <- lapply(seq_len(nrow(grid_q)), function(r_i) {
      cbind(
        grid_q[r_i, , drop = FALSE],
        empty_summary(error_considered = error_considered,
                      is_c = is_c,
                      algorithm = algorithm)
      )
    })
    eval_final_q_summary <- reorder_stats_columns(
      do.call("rbind", sum_blocks),
      error_considered = error_considered
    )
  } else {
    bgind <- which(data$calibration_data$pr_bg == 0)
    mods <- lapply(seq_along(data$part_data), function(i) {
      notrain <- -data$part_data[[i]]
      data_i <- data$calibration_data[notrain, , drop = FALSE]
      weights_i <- if (!is.null(weights)) weights[notrain] else NULL

      if (algorithm == "maxnet") {
        mod_i <- try(glmnet_mx(p = data_i$pr_bg, data = data_i,
                               f = formula_x, regmult = reg_x,
                               addsamplestobackground = addsamplestobackground,
                               weights = weights_i, only_fit = TRUE),
                     silent = TRUE)
      } else {
        mod_i <- try(glm_mx(formula = formula_x, family = binomial(link = "cloglog"),
                            data = data_i, weights = weights_i),
                     silent = TRUE)
      }

      if (inherits(mod_i, "try-error")) return(NULL)

      pred_i <- if (algorithm == "maxnet") {
        predict_mm(mod_i, mm_cal, type = "cloglog")
      } else {
        enmpa::predict_glm(model = mod_i, newdata = data$calibration_data,
                           type = "response")
      }

      suit_val_cal <- pred_i[unique(c(notrain, -bgind))]
      suit_val_eval <- pred_i[which(!-notrain %in% bgind)]
      om_rate <- omrat(threshold = error_considered, pred_train = suit_val_cal,
                       pred_test = suit_val_eval)

      if (proc_for_all) {
        proc_i <- unlist(lapply(error_considered, function(omr) {
          proc_omr <- fpROC::auc_metrics(test_prediction = suit_val_eval,
                                         prediction = pred_i,
                                         threshold = omr)$summary[, 4:5]
          names(proc_omr) <- c(paste0("Mean_AUC_ratio_at_", omr),
                               paste0("pval_pROC_at_", omr))
          proc_omr
        }))
      } else {
        proc_i <- rep(NA, length(error_considered) * 2L)
        names(proc_i) <- c(paste0("Mean_AUC_ratio_at_", error_considered),
                           paste0("pval_pROC_at_", error_considered))
      }

      df_eval_q <- if (algorithm == "maxnet") {
        data.frame(Partition = i, t(om_rate), t(proc_i), AICc = AICc,
                   Parameters = npar, Is_concave = is_c, row.names = NULL)
      } else {
        data.frame(Partition = i, t(om_rate), t(proc_i), AIC = m_aic$aic,
                   Parameters = npar, Is_concave = is_c, row.names = NULL)
      }
      cbind(grid_x, df_eval_q)
    })

    if (any(sapply(mods, is.null))) {
      is_c <- NA
      eval_final_q <- cbind(
        grid_x,
        empty_replicates(error_considered = error_considered,
                         n_row = length(data$part_data),
                         replicates = names(data$part_data),
                         is_c = is_c, algorithm = algorithm)
      )
      eval_final_q_summary <- reorder_stats_columns(
        cbind(grid_x, empty_summary(error_considered, is_c, algorithm)),
        error_considered = error_considered
      )
    } else {
      names(mods) <- names(data$part_data)
      eval_final_q <- do.call("rbind", mods)
      eval_final_q_summary <- reorder_stats_columns(
        eval_stats(eval_final_q, error_considered, algorithm),
        error_considered = error_considered
      )
    }
  }

  if (write_summary) {
    utils::write.csv(eval_final_q_summary,
                     file.path(out_dir, paste0("Summary_cand_model_", grid_x$ID, ".csv")),
                     row.names = FALSE)
  }

  if (!return_all_results) eval_final_q <- NULL
  return(list(All_results = eval_final_q, Summary = eval_final_q_summary))
}


fit_eval_models <- function(x, formula_grid, data, error_considered, omission_rate,
                            write_summary, addsamplestobackground, weights = NULL,
                            return_all_results, algorithm,
                            proc_for_all, out_dir) {
  grid_x <- formula_grid[x, , drop = FALSE]

  if (algorithm == "maxnet") {
    reg_x <- grid_x$R_multiplier
    formula_x <- stats::as.formula(grid_x$Formulas)
    m_aic <- try(glmnet_mx(p = data$calibration_data$pr_bg,
                           data = data$calibration_data,
                           f = formula_x, regmult = reg_x,
                           addsamplestobackground = addsamplestobackground,
                           weights = weights),
                 silent = TRUE)
  } else if (algorithm == "glm") {
    formula_x <- stats::as.formula(paste("pr_bg", grid_x$Formulas))
    m_aic <- try(glm_mx(formula = formula_x,
                        family = binomial(link = "cloglog"),
                        data = data$calibration_data,
                        weights = weights),
                 silent = TRUE)
  } else {
    stop("Unsupported model type. Please use 'maxnet' or 'glm'.")
  }

  if (inherits(m_aic, "try-error")) {
    npar <- NA
    AICc <- NA
    is_c <- NA
    mods <- NA
    class(mods) <- "try-error"
  } else {
    npar <- if (algorithm == "maxnet") {
      length(m_aic$betas)
    } else {
      length(m_aic$coefficients[-1])
    }

    if (algorithm == "maxnet") {
      mm_cal <- stats::model.matrix(formula_x, data$calibration_data)
      if (ncol(mm_cal) < 2L) mm_cal <- cbind(0, mm_cal)
      occ_rows <- which(data$calibration_data$pr_bg == 1)
      vals <- predict_mm(m_aic, mm_cal[occ_rows, , drop = FALSE],
                         type = "exponential")
      AICc <- aic_ws(pred_occs = vals, ncoefs = npar)
    }

    q_betas <- if (algorithm == "maxnet") {
      m_aic$betas[grepl("\\^2", names(m_aic$betas))]
    } else {
      m_aic$coefficients[-1][grepl("\\^2", names(m_aic$coefficients[-1]))]
    }

    is_c <- if (length(q_betas) == 0L) FALSE else any(q_betas > 0)
    bgind <- which(data$calibration_data$pr_bg == 0)

    mods <- try(lapply(seq_along(data$part_data), function(i) {
      notrain <- -data$part_data[[i]]
      data_i <- data$calibration_data[notrain, , drop = FALSE]
      weights_i <- if (!is.null(weights)) weights[notrain] else NULL

      if (algorithm == "maxnet") {
        mod_i <- try(glmnet_mx(p = data_i$pr_bg, data = data_i,
                               f = formula_x, regmult = reg_x,
                               addsamplestobackground = addsamplestobackground,
                               weights = weights_i, only_fit = TRUE),
                     silent = TRUE)
      } else {
        mod_i <- try(glm_mx(formula = formula_x,
                            family = binomial(link = "cloglog"),
                            data = data_i, weights = weights_i),
                     silent = TRUE)
      }

      if (inherits(mod_i, "try-error")) return(NULL)

      pred_i <- if (algorithm == "maxnet" && inherits(mod_i, "glmnet_mx")) {
        predict_mm(mod_i, mm_cal, type = "cloglog")
      } else {
        enmpa::predict_glm(model = mod_i, newdata = data$calibration_data,
                           type = "response")
      }

      suit_val_cal <- pred_i[unique(c(notrain, -bgind))]
      suit_val_eval <- pred_i[which(!-notrain %in% bgind)]
      om_rate <- omrat(threshold = error_considered, pred_train = suit_val_cal,
                       pred_test = suit_val_eval)

      if (proc_for_all) {
        proc_i <- unlist(lapply(error_considered, function(omr) {
          proc_omr <- fpROC::auc_metrics(test_prediction = suit_val_eval,
                                         prediction = pred_i,
                                         threshold = omr)$summary[, 4:5]
          names(proc_omr) <- c(paste0("Mean_AUC_ratio_at_", omr),
                               paste0("pval_pROC_at_", omr))
          proc_omr
        }))
      } else {
        proc_i <- rep(NA, length(error_considered) * 2L)
        names(proc_i) <- c(paste0("Mean_AUC_ratio_at_", error_considered),
                           paste0("pval_pROC_at_", error_considered))
      }

      df_eval <- if (algorithm == "maxnet") {
        data.frame(Partition = i, t(om_rate), t(proc_i), AICc = AICc,
                   Parameters = npar, Is_concave = is_c, row.names = NULL)
      } else {
        data.frame(Partition = i, t(om_rate), t(proc_i), AIC = m_aic$aic,
                   Parameters = npar, Is_concave = is_c, row.names = NULL)
      }
      cbind(grid_x, df_eval)
    }), silent = TRUE)
  }

  if (inherits(mods, "try-error") || any(sapply(mods, is.null))) {
    is_c <- NA
    eval_final <- cbind(grid_x,
                        empty_replicates(error_considered = error_considered,
                                         n_row = length(data$part_data),
                                         replicates = names(data$part_data),
                                         is_c = is_c, algorithm = algorithm))
    eval_final_summary <- reorder_stats_columns(
      cbind(grid_x, empty_summary(error_considered, is_c, algorithm)),
      error_considered = error_considered
    )
  } else {
    names(mods) <- names(data$part_data)
    eval_final <- do.call("rbind", mods)
    eval_final_summary <- reorder_stats_columns(
      eval_stats(eval_final, error_considered, algorithm),
      error_considered = error_considered
    )
  }

  if (write_summary) {
    utils::write.csv(eval_final_summary,
                     file.path(out_dir, paste0("Summary_cand_model_", grid_x$ID, ".csv")),
                     row.names = FALSE)
  }

  if (!return_all_results) eval_final <- NULL
  return(list(All_results = eval_final, Summary = eval_final_summary))
}


fit_best_model <- function(x, dfgrid, cal_res, n_replicates = 1,
                           rep_data = NULL, algorithm = "maxnet") {
  if (missing(x)) stop("Argument 'x' must be defined.")
  if (missing(dfgrid)) stop("Argument 'dfgrid' must be defined.")
  if (missing(cal_res)) stop("Argument 'cal_res' must be defined.")
  if (!inherits(dfgrid, "data.frame")) stop("Argument 'dfgrid' must be a 'data.frame'.")
  if (!inherits(cal_res, "calibration_results")) stop("Argument 'cal_res' must be a 'calibration_results' object.")
  if (!is.numeric(n_replicates)) stop("Argument 'n_replicates' must be 'numeric'.")
  if (n_replicates > 1 && !inherits(rep_data, "list")) stop("Argument 'rep_data' must be a 'list'.")
  if (!inherits(algorithm, "character") || !(algorithm %in% c("glm", "maxnet"))) {
    stop("Argument 'algorithm' must be 'glm' or 'maxnet'.")
  }

  grid_x <- dfgrid[x, ]
  m_id <- grid_x$models
  rep_x <- grid_x$replicates

  best_model <- cal_res$selected_models[cal_res$selected_models$ID == m_id, ]
  best_formula <- best_model$Formulas

  if (algorithm == "maxnet") {
    best_regm <- best_model$R_multiplier
  }

  if (n_replicates > 1) {
    rep_i <- rep_data[[rep_x]]
    data_x <- cal_res$calibration_data[-rep_i, , drop = FALSE]
  } else {
    data_x <- cal_res$calibration_data
  }

  if (algorithm == "maxnet") {
    mod_x <- glmnet_mx(p = data_x$pr_bg, data = data_x,
                       f = stats::as.formula(best_formula),
                       regmult = best_regm,
                       addsamplestobackground = cal_res$addsamplestobackground,
                       weights = cal_res$weights,
                       calculate_AIC = FALSE)
  } else {
    mod_x <- glm_mx(formula = stats::as.formula(paste("pr_bg", best_formula)),
                    family = binomial(link = "cloglog"),
                    data = data_x,
                    weights = cal_res$weights)
  }

  mod_x$checkModel <- m_id
  mod_x$checkPartition <- rep_x

  return(mod_x)
}


bind_rows_projection <- function(data_frames) {
  all_columns <- c("Time", "Period", "Scenario", "ssp", "GCM",
                   "input_path", "output_path")
  result <- NULL

  for (df in data_frames) {
    if (!is.null(df)) {
      missing_columns <- setdiff(all_columns, colnames(df))
      for (col in missing_columns) {
        df[[col]] <- NA
      }
      result <- rbind(result, df[, all_columns])
    }
  }

  result <- result[, colSums(is.na(result)) < nrow(result)]
  return(result)
}


reorder_stats_columns <- function(stats_final, error_considered) {
  if (!inherits(stats_final, "data.frame")) {
    stop("Argument 'stats_final' must be a 'data.frame'.")
  }
  if (!inherits(error_considered, "numeric")) {
    stop("Argument 'error_considered' must be 'numeric'.")
  }

  first_cols <- intersect(c("ID", "Formulas", "R_multiplier", "Features"),
                          colnames(stats_final))
  metric_cols <- c(paste0("Omission_rate_at_", error_considered, ".mean"),
                   paste0("Omission_rate_at_", error_considered, ".sd"),
                   paste0("Mean_AUC_ratio_at_", error_considered, ".mean"),
                   paste0("Mean_AUC_ratio_at_", error_considered, ".sd"),
                   paste0("pval_pROC_at_", error_considered, ".mean"),
                   paste0("pval_pROC_at_", error_considered, ".sd"))
  ordered_metric_cols <- unlist(lapply(error_considered, function(rate) {
    grep(paste0("_", rate, "\\."), metric_cols, value = TRUE)
  }))
  last_cols <- setdiff(colnames(stats_final), c(first_cols, metric_cols))
  orders_cols <- c(first_cols, ordered_metric_cols, last_cols)
  return(stats_final[, orders_cols])
}


#### PROC ####
proc <- function(x, formula_grid, data, error_considered = 10,
                 addsamplestobackground = TRUE, weights = NULL,
                 algorithm) {
  if (missing(x)) stop("Argument 'x' must be defined.")
  if (missing(formula_grid)) stop("Argument 'formula_grid' must be defined.")
  if (missing(data)) stop("Argument 'data' must be defined.")
  if (!inherits(formula_grid, "data.frame")) stop("Argument 'formula_grid' must be a 'data.frame'.")
  if (!inherits(error_considered, "numeric")) stop("Argument 'error_considered' must be 'numeric'.")
  if (!inherits(addsamplestobackground, "logical")) stop("Argument 'addsamplestobackground' must be 'logical'.")
  if (!is.null(weights) && !inherits(weights, "numeric")) stop("Argument 'weights' must be NULL or 'numeric'.")
  if (!inherits(algorithm, "character") || !(algorithm %in% c("glm", "maxnet"))) {
    stop("Argument 'algorithm' must be 'glm' or 'maxnet'.")
  }

  grid_x <- formula_grid[x, ]
  formula_x <- stats::as.formula(grid_x$Formulas)
  reg_x <- grid_x$R_multiplier

  if (algorithm == "glm") {
    formula_x <- stats::as.formula(paste("pr_bg", grid_x$Formulas))
  } else {
    mm_cal <- stats::model.matrix(formula_x, data$calibration_data)
    if (ncol(mm_cal) < 2L) mm_cal <- cbind(0, mm_cal)
  }

  bgind <- which(data$calibration_data$pr_bg == 0)

  mods <- try(lapply(seq_along(data$part_data), function(i) {
    notrain <- -data$part_data[[i]]
    data_i <- data$calibration_data[notrain, , drop = FALSE]
    weights_i <- if (!is.null(weights)) weights[notrain] else NULL

    if (algorithm == "maxnet") {
      mod_i <- glmnet_mx(p = data_i$pr_bg, data = data_i,
                         f = formula_x, regmult = reg_x,
                         addsamplestobackground = addsamplestobackground,
                         weights = weights_i, only_fit = TRUE)
      pred_i <- predict_mm(mod_i, mm_cal, type = "cloglog")
    } else {
      mod_i <- glm_mx(formula = formula_x,
                      family = binomial(link = "cloglog"),
                      data = data_i, weights = weights_i)
      pred_i <- enmpa::predict_glm(model = mod_i,
                                   newdata = data$calibration_data,
                                   type = "response")
    }

    suit_val_cal <- pred_i[unique(c(notrain, -bgind))]
    suit_val_eval <- pred_i[which(!-notrain %in% bgind)]

    proc_i <- unlist(lapply(error_considered, function(omr) {
      proc_omr <- suppressWarnings(
        fpROC::auc_metrics(test_prediction = suit_val_eval,
                           prediction = pred_i,
                           threshold = omr)[[1]]
      )

      if (all(is.na(proc_omr))) {
        proc_omr <- data.frame(proc_omr[1], proc_omr[2])
      } else {
        proc_omr <- proc_omr[, 4:5]
      }

      names(proc_omr) <- c(paste0("Mean_AUC_ratio_at_", omr),
                           paste0("pval_pROC_at_", omr))
      proc_omr
    }))

    data.frame(Partition = i,
               t(proc_i),
               row.names = NULL)
  }), silent = TRUE)

  proc_df <- do.call("rbind", mods)

  means <- colMeans(proc_df[, -1, drop = FALSE], na.rm = TRUE)
  sds <- apply(proc_df[, -1, drop = FALSE], 2, stats::sd, na.rm = TRUE)

  proc_df <- data.frame(t(c(means, sds)))
  names(proc_df) <- c(paste0(names(means), ".mean"), paste0(names(sds), ".sd"))
  proc_df$ID <- grid_x$ID

  return(proc_df)
}
