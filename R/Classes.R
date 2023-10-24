new_explore_back <- function(hist_m, hist_b, hist_o, mean_m, mean_b, mean_o,
                             cl_m, cl_b, cl_o, range_m, range_b, range_o) {
  # error checking

  val <- list(hist_m = hist_m, hist_b = hist_b, hist_o = hist_o,
              mean_m = mean_m, mean_b = mean_b, mean_o = mean_o,
              cl_m = cl_m, cl_b = cl_b, cl_o = cl_o, range_m = range_m,
              range_b = range_b, range_o = range_o)

  class(val) <- "explore_back"
  return(val)
}


new_explore_list <- function(summary, exploration_stats) {
  # error checking

  val <- list(summary = summary, exploration_stats = exploration_stats)

  class(val) <- "explore_list"
  return(val)
}


new_occ_cal <- function(summary, occurrences, block_limits = NULL,
                        calibration_regime = NULL) {
  # error checking

  val <- list(summary = summary, occurrences = occurrences,
              block_limits = block_limits,
              calibration_regime = calibration_regime)
  class(val) <- "occ_cal"
  return(val)
}



new_back_cal <- function(summary, background, background_bias = NULL,
                         training_regime = NULL) {
  # error checking

  val <- list(summary = summary, background = background,
              background_bias = background_bias,
              training_regime = training_regime)
  class(val) <- "back_cal"
  return(val)
}

new_lambdas <- function(lambdas_path){
  lambdas <- readLines(lambdas_path)
  maxcalcsV <- c("linearPredictorNormalizer",
                 "densityNormalizer",
                 "numBackgroundPoints",
                 "entropy")
  maxids <-  sapply(seq_along(maxcalcsV), function(x){
    ids <- stringr::str_detect(string = lambdas,maxcalcsV[x])
    return(which(ids))
  })
  if(length(maxids)==0L){
    stop("Please provide a valid lambdas file")
  }
  maxcalcs <- lapply(maxids, function(x){
    val <- stringr::str_split(lambdas[x],", ")[[1]]
    v <- as.numeric(val[2])
    return(v)
  })

  names(maxcalcs) <- maxcalcsV
  lambvals <- data.frame(stringr::str_split(string = lambdas[-maxids],
                                            pattern = "[,]",simplify = T))
  lambvals[,-1] <- as.numeric(as.matrix(lambvals[,-1]))
  names(lambvals) <- c("variable","lambda","min","max")
  feature <- sapply(lambvals$variable, function(x){
    qq <- which(stringr::str_detect(string =x ,"\\^"))
    pd <- which(stringr::str_detect(string =x ,"\\*"))
    fh <- which(stringr::str_detect(string =x ,"\\'"))
    rh <- which(stringr::str_detect(string =x ,"\\`"))
    th <- which(stringr::str_detect(string =x ,"\\<"))
    if(length(qq)){
      return("quadratic")
    } else if(length(pd)){
      return("product")
    } else if(length(fh)){
      return("forward_hinge")
    } else if(length(rh)){
      return("reverse_hinge")
    } else if(length(rh)){
      return("reverse_hinge")
    } else if(length(th)){
      return("threshold")
    } else{
      return("linear")
    }
  })
  lambvals$feature <- feature
  res <- list(lambdas_df=lambvals,maxmeta=maxcalcs)
  class(res) <- c("lambdas")
  return(res)
}
