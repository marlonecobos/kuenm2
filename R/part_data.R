# data <- read.csv("Models/Piper_fuligineum/occ_bg.csv")
# data$folds <- NULL
part_data <- function(data,
                      pr_bg = "pr_bg",
                      train_portion = 0.7,
                      n_replicates = 5,
                      method = "subsample") {
  #Get data
  d <- data[pr_bg]
  #Split presence and absence
  pre <- which(d[, pr_bg] == 1)
  aus <- which(d[, pr_bg] == 0)

  if(method == "kfold") {
  foldp <- sample(cut(seq(1, length(pre)), breaks = n_replicates, labels = FALSE))
  folda <- sample(cut(seq(1, length(aus)), breaks = n_replicates, labels = FALSE))
  #Join data
  d$folds <- NA
  d$folds[which(d[,pr_bg] == 1)] <- foldp
  d$folds[which(d[,pr_bg] == 0)] <- folda

  rep_data <- lapply(unique(d$folds), function(f) {
    which(d$folds != f) })
  names(rep_data) <- paste0("Rep_", 1:n_replicates)
  }

  if(method == "subsample" | method == "bootstrap") {
      if(method == "subsample"){
      replacement = FALSE
    } else if(method == "bootstrap"){
      replacement = TRUE
    }


  #Create list of replicates
  rep_data <- lapply(1:n_replicates, function(i){
  set.seed(42*i)
  foldp <- sample(pre,
                  size = floor(train_portion * length(pre)),
                  replace = FALSE) #Always false
  folda <- sample(aus,
                  size = floor(train_portion * length(aus)),
                  replace = FALSE) #Always false
  foldpa <- c(foldp, folda)
  return(foldpa)
  })
  names(rep_data) <- paste0("Rep_", 1:n_replicates)
  }
return(rep_data)
} #End of function



# #Test function
data <- read.csv("Models/Piper_fuligineum/occ_bg.csv")
a <- part_data(data = data, train_portion = 0.75,
               n_replicates = 5,
               method = "subsample")
