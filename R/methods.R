predict_maxent <- function(object,newdata,type="logistic",...){
  if(!inherits(object,"lambdas")){
    lambvals <- new_lambdas(lambdas_path = object)
  } else{
    lambvals <- object
  }
  eval_layers <- FALSE
  if(methods::is(newdata,"RasterStack") || methods::is(newdata,"RasterBrick")){
    eval_layers <- TRUE
    varnames <- names(newdata)
    rmask <- newdata[[1]]

    newdata <- lapply(1:raster::nlayers(newdata), function(x){
      rval <- newdata[[x]][]
      #rval <- rval[nona]
      return(rval)
    })
    newdata <- do.call(cbind,newdata)
    nona <- which(complete.cases(newdata))
    newdata <- newdata[nona,]
    colnames(newdata) <- varnames
  }
  # Function to estimate maxents raw output
  eval_lamb <- function(lambvals,bg_vars){
    ll <- lambvals[which(lambvals$lambda !=0),]
    lls <- ll %>% split(.$feature)
    fxL <- seq_along(lls) %>% purrr::map(function(j){
      fx <- 0
      if(names(lls)[j] %in% "forward_hinge"){
        vars <- gsub("\\'",replacement = "",lls[[j]]$variable)
        for (i in seq_along(vars)){
          x <- bg_vars[,which(colnames(bg_vars) %in% vars[i])]
          lambda <- lls[[j]]$lambda[i]
          minn <- lls[[j]]$min[i]
          maxx <- lls[[j]]$max[i]
          fx <- fx + ifelse(x < minn,0, lambda*(x - minn)/(maxx - minn))
        }
      }
      if(names(lls)[j] %in% "reverse_hinge"){
        vars <- gsub(pattern = "\\`",replacement = "",lls[[j]]$variable)
        for (i in seq_along(vars)) {
          x <- bg_vars[,which(colnames(bg_vars) %in% vars[i])]
          lambda <- lls[[j]]$lambda[i]
          minn <- lls[[j]]$min[i]
          maxx <- lls[[j]]$max[i]
          fx <- fx + ifelse(x < maxx,lambda*(maxx - x)/(maxx - minn) , 0)
        }
      }
      if(names(lls)[j] %in% "product"){
        vars <- stringr::str_split(string = lls[[j]]$variable,pattern = "\\*")
        for(i in seq_along(vars)){
          x <- bg_vars[,vars[[i]][1]]*bg_vars[,vars[[i]][2]]
          lambda <- lls[[j]]$lambda[i]
          minn <- lls[[j]]$min[i]
          maxx <- lls[[j]]$max[i]
          fx <- fx +  lambda*(x - minn)/(maxx - minn)
        }
      }
      if(names(lls)[j] %in% "quadratic"){
        vars <- gsub(pattern = "\\^2",replacement = "",lls[[j]]$variable)
        for (i in seq_along(vars)) {
          x <- bg_vars[,which(colnames(bg_vars) %in% vars[i])]
          lambda <- lls[[j]]$lambda[i]
          minn <- lls[[j]]$min[i]
          maxx <- lls[[j]]$max[i]
          fx <- fx +  lambda*(x*x - minn)/(maxx - minn)
        }
      }
      if(names(lls)[j] %in% "threshold"){
        vars1 <- stringr::str_split(string = lls[[j]]$variable,
                                    pattern = "\\<",simplify = T)
        vars <- gsub("\\)","",vars1[,2])

        thvals <- as.numeric(gsub("\\(","",vars1[,1]))
        for (i in seq_along(vars)) {
          x <- bg_vars[,which(colnames(bg_vars) %in% vars[i])]
          lambda <- lls[[j]]$lambda[i]
          minn <- lls[[j]]$min[i]
          maxx <- lls[[j]]$max[i]
          fx <- fx + ifelse(x < thvals[i],0,lambda)
        }
      }
      if(names(lls)[j] %in% "linear"){
        vars <- lls[[j]]$variable
        for (i in seq_along(vars)) {
          x <- bg_vars[,which(colnames(bg_vars) %in% vars[i])]
          lambda <- lls[[j]]$lambda[i]
          minn <- lls[[j]]$min[i]
          maxx <- lls[[j]]$max[i]
          fx <- fx + lambda*(x - minn)/(maxx - minn)
        }
      }
      return(fx)
    })
    Reduce("+", fxL)
  }
  maxcalcs <- lambvals$maxmeta
  maxv <- eval_lamb(lambvals = lambvals$lambdas_df,bg_vars = newdata) -
    maxcalcs$linearPredictorNormalizer
  rm(newdata)
  raw_max <- exp(maxv)/maxcalcs$densityNormalizer
  vals <- switch(
    type,
    "raw" = raw_max,
    "logistic" = raw_max*exp(maxcalcs$entropy)/
      (1 + raw_max*exp(maxcalcs$entropy)),
    "cloglog" =  1 - exp(-exp(maxcalcs$entropy) * raw_max)
  )
  rm(raw_max)
  if(eval_layers){
    predv <- rep(NA,raster::ncell(rmask))
    predv[nona] <- vals
    rmask[] <- predv
    vals <- rmask
    rm(rmask)
  }

  return(vals)
}

#' Print Method for fitted_models Class
#' @export
print.fitted_models <- function(x, ...){
  cat("fitted_models object summary\n")
  cat("==========================\n")
  cat("Species:", x$species, "\n")

  cat("Number of fitted models:", length(x$Models), "\n")
  #Get number of replicates
  nr <- sum(grepl("Rep", names(x$Models[[1]])))
  cat("Models fitted with", nr, "replicates")
}
