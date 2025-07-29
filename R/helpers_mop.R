#Function to run mop on new records
mop_with_records <- function(train_data, test_data, variables = NULL,
                             categorical_variables = NULL,
                             mop_type = "detailed", calculate_distance = TRUE,
                             where_distance = "all",
                             progress_bar = TRUE, ...){

  if (missing(train_data)) {
    stop("Argument 'train_data' must be defined.")
  }
  if (!inherits(train_data, "data.frame")) {
    stop("Argument 'train_data' must be a 'data.frame' object.")
  }
  if (missing(test_data)) {
    stop("Argument 'test_data' must be defined.")
  }
  if (!inherits(test_data, "data.frame")) {
    stop("Argument 'test_data' must be a 'data.frame' object.")
  }

  if (!inherits(mop_type, "character")) {
    stop("Argument 'mop_type' must be a 'character'.")
  }

  mop_type_out <- setdiff(mop_type, c("basic", "simple", "detailed"))
  if (length(mop_type_out) > 0) {
    stop("Invalid 'mop_type' provided.",
         "\nAvailable options are: 'basic', 'simple', or 'detailed'.")
  }

  if(!inherits(calculate_distance, "logical")){
    stop("Argument 'calculate_distance' must be 'logical'.")
  }

  distance_out <- setdiff(where_distance, c("in_range", "out_range", "all"))
  if (length(distance_out) > 0) {
    stop("Invalid 'where_distance' provided.",
         "\nAvailable options are: 'in_range', 'out_range', and 'all'.")
  }

  if(!inherits(progress_bar, "logical")){
    stop("Argument 'progress_bar' must be 'logical'.")
  }

  if(is.null(variables)){
    v <- colnames(train_data)
    v <- setdiff(v, c("pr_bg", "x", "y"))
  } else{
    v <- variables
  }

  # remove categorical variables
  if (!is.null(categorical_variables)) {
    v <- setdiff(v, categorical_variables)
  }

  #Identify different categories...
  mop_res <- mop::mop(m = as.matrix(train_data[,v]),
                      g = as.matrix(test_data[,v]),
                      type = mop_type,
                      calculate_distance = calculate_distance,
                      where_distance = where_distance,
                      progress_bar = progress_bar, ...)

  #Get results in data.frame
  df <- list()
  #Distance
  if(calculate_distance){
    df[["mop_distance"]] <- mop_res$mop_distances
  }
  #Mop basic
  df[["inside_range"]] <- is.na(mop_res$mop_basic)
  #Mop simple
  if(mop_type %in% c("simple", "detailed")){
    n_var <- mop_res$mop_simple
    n_var[is.na(n_var)] <- 0
    df[["n_var_out"]] <- n_var
  }
  #Mop detailed
  if(mop_type == "detailed"){
    #Create map to get variables names
    interpretation_map <- stats::setNames(
      mop_res[["mop_detailed"]][["interpretation_combined"]]$extrapolation_variables,
      mop_res[["mop_detailed"]][["interpretation_combined"]]$values
    )
    #Replace values in toward high and low
    towards_low <- interpretation_map[as.character(
      mop_res[["mop_detailed"]][["towards_low_combined"]])]
    towards_high <- interpretation_map[as.character(
      mop_res[["mop_detailed"]][["towards_high_combined"]])]
    #Append results
    df[["towards_low"]] <- towards_low
    df[["towards_high"]] <- towards_high
  } #End of detailed
  #Convert to dataframe
  df <- as.data.frame(df)

  #Get results with categorical
  if(!is.null(categorical_variables)){
    cat_train <- sapply(categorical_variables, function(ct){
      unique(train_data[[ct]])}, simplify = FALSE)
    cat_test <- sapply(categorical_variables, function(ct){
      unique(test_data[[ct]])}, simplify = FALSE)
    #Get differences
    cat_diff <- sapply(categorical_variables, function(ct){
      setdiff(cat_test[[ct]], cat_train[[ct]])
    }, simplify = FALSE)

    #Append results
    for(ct in categorical_variables){
      #Create column
      df[[ct]] <- NA
      ct_x <- cat_diff[[ct]]
      #Get occurrences with these values
      for(x in ct_x){
        df[test_data[[ct]] == x, ct] <- x
      }

      #Update inside_range, if necessary
      df[!is.na(df[[ct]]), "inside_range"] <- FALSE
    }

    if(mop_type %in% c("simple", "detailed")){
      #Compute number of categorical variables outside the range
      n_cat_out <- rowSums(1 - is.na(as.data.frame(df[,categorical_variables])))
      df[["n_var_out"]] <- df[["n_var_out"]] + n_cat_out
    }
  } #End of !is.null(categorical_variables)


  return(list("mop_results" = mop_res,
              "mop_records" = df))
}

