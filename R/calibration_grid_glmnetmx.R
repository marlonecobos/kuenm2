####Get grid of formulas####
calibration_grid_glmnetmx <- function(var_names = NULL, swd = NULL, x = NULL, y = NULL,
                                      fold_column = "folds",
                                      min_number = 2,
                                      categorical_var = NULL,
                                      features = c("l", "q", "lq", "lqp", "p"),
                                      min_continuous = NULL,
                                      regm = c(0.1, 1, 2, 3, 5)){
  if(is.null(var_names) & !is.null(swd)) {
    ignore_columns <- intersect(c(x, y, "pr_bg", "species", fold_column),
                                colnames(swd))
    var_names <- colnames(swd[, !names(swd) %in% ignore_columns])}

  #Get var combinations
  var_comb <- enmpa:::aux_var_comb(var_names = var_names,
                                   minvar = min_number)

  #Split features
  formula_x <- unlist(lapply(seq_along(features), function(x){
    f_x <- features[x]
    if(grepl("p", f_x) & !is.null(categorical_var)) {
      var_comb_new <- var_comb[sapply(var_comb,
                                      function(x)
                                        sum(!grepl(categorical_var, x))) >= 2]
    } else {
      var_comb_new <- var_comb}
    ff_x <- unlist(lapply(seq_along(var_comb_new), function(i){
      #If type = p, get only combinations with 2 or more combinations

      f_l <- prepare_formulas_glmnetmx(independent = var_comb_new[[i]], type = f_x,
                                       categorical_var = categorical_var)
      names(f_l) <- f_x
      return(f_l) }))
  }))
  #Create dataframe with formulas
  formula_d <- data.frame(formula = paste(formula_x, -1),
                          features = names(formula_x))
  #Expand grid
  f_grid <- unique(expand.grid(formula_d$formula, regm))
  colnames(f_grid) <- c("Formulas", "regm")
  f_grid <- merge(f_grid, formula_d, by.x = "Formulas", by.y = "formula")
  f_grid$ID <- 1:nrow(f_grid)
  colnames(f_grid) <- c("Formulas", "regm", "Features", "ID")
  #Reorder columns
  f_grid <- f_grid[, c("ID", "Formulas", "regm", "Features")]
  #Class formulas
  f_grid$Formulas <- as.character(f_grid$Formulas)
  return(f_grid)
}

# #Test
# my_grid <- calibration_grid_glmnetmx(var_names = c("bio_1", "bio_7", "bio_12",
#                                                   "bio_15", "soil"),
#                       min_number = 2, categorical_var = "soil",
#                       features =  c("l", "q", "lq", "lqp", "p"),
#                       regm = c(0.1, 1, 2, 3, 5))

####Get formulas####
prepare_formulas_glmnetmx <- function(independent, type = "lqpht",
                                      categorical_var = NULL,
                                      minvar = 1, maxvar = NULL)
{
  if (is.character(type)) {
    if (!all(unlist(strsplit(type, "")) %in% c("l", "p",
                                               "q", "h", "t"))) {
      stop("'type' must be: 'l', 'p', 'q', 'h', 't', or a combination of those three.")
    }
  }  else {
    stop("'type' must be of class character.")
  }

  #Check if categorical variable is in independente variables. If not, set null
  if(!is.null(categorical_var)) {
    if(!(categorical_var %in% independent)) {
      categorical_var <- NULL
    }  }

  if(!is.null(categorical_var)) {
    independent <- setdiff(independent, categorical_var)
  }

  predictors <- independent
  npred <- length(predictors)
  aux <- " "
  #Linear
  if (grepl("l", type)) {
    aux <- paste(aux, paste(predictors, collapse = " + "),
                 sep = " + ")
  }
  #Quadratic
  if (grepl("q", type)) {
    for (v in 1:length(predictors)) {
      aux <- paste(aux, paste0("I(", predictors[v], "^2)"),
                   sep = " + ")
    }
  }
  #Product
  if (grepl("p", type)) {
    if (npred > 1) {
      inter_tab <- utils::combn(predictors, 2)
      aux_inter <- paste0(" + ", paste(apply(inter_tab,
                                             2, paste, collapse = ":"), collapse = " + "))
      if (length(aux_inter) > 0) {
        aux <- paste0(aux, aux_inter)
      }
    }
    else {
      if (grepl("l", type) | grepl("q", type)) {
        message("'p' is is only possible with 2 or more independent variables.",
                "\nReturning other combinations.")
      }
      else {
        stop("'p' is is only possible with 2 or more independent variables.",
             "\nTry other combinations of type.")
      }
    }
  }
  #Hinge
  if (grepl("h", type)) {
    for (v in 1:length(predictors)) {
      aux <- paste(aux, paste0("hinge(", predictors[v], ")"),
                   sep = " + ")
    }
  }
  #Threshold
  if (grepl("t", type)) {
    for (v in 1:length(predictors)) {
      aux <- paste(aux, paste0("thresholds(", predictors[v], ")"),
                   sep = " + ")
    }
  }
  out <- paste0("~",
                gsub(pattern = "  \\+ ", x = aux, replacement = ""))

  if(!is.null(categorical_var)) {
    out <- paste0(out, " + categorical(", categorical_var, ")")
  }
  return(out)
}

#my_formula <- prepare_formulas_glmnetmx(independent = c("PC1", "PC2"), type = "lq")
