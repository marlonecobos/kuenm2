####Helpers from eval_m ####

#Function to summarize data
eval_stats <- function(calib_results, toagg = NULL, agg_by = NULL,
                       to_keep = NULL){
  if(is.null(toagg)){
    toagg <- c("Omission_rate_at_5", "Omission_rate_at_10",
               "proc_auc_ratio", "proc_pval") }
  if(is.null(agg_by)){
    agg_by <- c("Formulas", "regm", "Features")
  }
  if(is.null(to_keep)){
    to_keep <- c("ID", "Formulas", "regm", "Features", "AIC_nk",
                 "AIC_ws", "npar", "is_convex")
  }
  agg_formula <- paste("~", paste(agg_by, collapse = " + "))

  #Get summary
  xy <- lapply(toagg, function(x) {
    do.call(data.frame, stats::aggregate(as.formula(paste(x, agg_formula)),
                                         data = calib_results, FUN = function(y) c(mean = round(mean(y), 4), sd = round(sd(y), 4)), na.action=NULL))
  })

  #Summarize stats
  stats <- Reduce(function(x, y) merge(x, y,
                                       by = agg_by),
                  xy)

  stats_AICS <- calib_results[!duplicated(calib_results[,to_keep]),][,to_keep]
  stats_final <- merge(stats, stats_AICS, by = agg_by)
  return(stats_final)
}


####Create empty dataframes####
empty_replicates <- function(omrat_thr = c(5, 10),
                             n_row = 4, replicates = 1:4,
                             is_c = NA) {
  column_names <- c("Replicate", paste0("Omission_rate_at_", omrat_thr),
                    "proc_auc_ratio", "proc_pval", "AIC_nk", "AIC_ws", "npar", "is_convex")
  df_eval_q <- data.frame(matrix(NA, nrow = n_row, ncol = length(column_names)))
  colnames(df_eval_q) <- column_names
  df_eval_q$Replicate <- replicates
  df_eval_q$is_convex = is_c
  return(df_eval_q)
}

empty_summary <- function(omrat_thr, is_c){
  om_means <- paste0("Omission_rate_at_", omrat_thr, ".mean")
  om_sd <- paste0("Omission_rate_at_", omrat_thr, ".sd")
  column_names <- c(om_means, om_sd,
                    "proc_auc_ratio.mean", "proc_auc_ratio.sd", "proc_pval.mean",
                    "proc_pval.sd", "AIC_nk", "AIC_ws", "npar", "is_convex")
  eval_final_q  <- data.frame(matrix(NA, nrow = 1, ncol = length(column_names)))
  colnames(eval_final_q) <- column_names
  eval_final_q$is_convex = is_c
  return(eval_final_q)
}
