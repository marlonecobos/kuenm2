# 5. Project Models to a Single Scenario

- [Introduction](#introduction)
- [Predict selected models for a single
  scenario](#predict-selected-models-for-a-single-scenario)
  - [Predict to SpatRaster](#predict-to-spatraster)
  - [Predict to data.frame](#predict-to-data.frame)
  - [Binarize models](#binarize-models)
  - [Clamping variables](#clamping-variables)
  - [No Extrapolation](#no-extrapolation)
  - [Output Type](#output-type)
- [Saving Predictions](#saving-predictions)

------------------------------------------------------------------------

### Introduction

Once selected models have been fit using
[`fit_selected()`](https://marlonecobos.github.io/kuenm2/reference/fit_selected.md),
projections to single or multiple scenarios can be performed. The
[`predict_selected()`](https://marlonecobos.github.io/kuenm2/reference/predict_selected.md)
function is designed for projections to single scenarios.

To predict using the selected models, a `fitted_models` object is
required. For detailed information on model fitting, please consult the
vignette [Fit and Explore Selected
Models](https://marlonecobos.github.io/kuenm2/articles/model_exploration.md).
The `fitted_models` object generated in that vignette is included as an
example dataset within the package. Let’s load it.

``` r
#Load packages
library(kuenm2)
library(terra)
#> terra 1.8.86

#Import calib_results_maxnet
data("fitted_model_maxnet", package = "kuenm2")
#Print calibration result
fitted_model_maxnet
#> fitted_models object summary
#> ============================
#> Species: Myrcia hatschbachii 
#> Algortihm: maxnet 
#> Number of fitted models: 2 
#> Only full models fitted, no replicates
```

To compare the results, let’s import a `fitted_models` object generated
using the GLM algorithm:

``` r
#Import calib_results_maxnet
data("fitted_model_glm", package = "kuenm2")
#Print calibration result
fitted_model_glm
#> fitted_models object summary
#> ============================
#> Species: Myrcia hatschbachii 
#> Algortihm: glm 
#> Number of fitted models: 1 
#> Only full models fitted, no replicates
```

  

## Predict Selected Models for a Single Scenario

To predict selected models for a single scenario, you need a
`fitted_models` object and the corresponding predictor variables. These
predictor variables can be provided as either a `SpatRaster` or a
`data.frame`. The names of the variables (or columns in the
`data.frame`) must precisely match those used for model calibration or
those used when running PCA (if `do_pca = TRUE` was set in the
[`prepare_data()`](https://marlonecobos.github.io/kuenm2/reference/prepare_data.md)
function; see [Prepare Data for Model
Calibration](https://marlonecobos.github.io/kuenm2/articles/prepare_data.md)
for more details).

### Predict to SpatRaster

Let’s use the same raster variables that were used to prepare the data
and calibrate the models. These are included as example data within the
package:

``` r
# Import raster layers
var <- rast(system.file("extdata", "Current_variables.tif", package = "kuenm2"))
#Plot raster layers
plot(var)
```

![](model_predictions_files/figure-html/Import%20raster%20layers-1.png)

Let’s check which variables were used to calibrate our models. They are
available in the `calibration_data` element of the object:

``` r
# Variables used to calibrate maxnet models
colnames(fitted_model_maxnet$calibration_data)
#> [1] "pr_bg"    "bio_1"    "bio_7"    "bio_12"   "bio_15"   "SoilType"

#Variables used to calibrate glm models
colnames(fitted_model_glm$calibration_data)
#> [1] "pr_bg"    "bio_1"    "bio_7"    "bio_12"   "bio_15"   "SoilType"
```

The first column, “pr_bg”, indicates the presence (1) and absence (0)
records, while the other columns represent the environmental variables.
In this case, the variables are `bio_1`, `bio_7`, `bio_12`, `bio_15`,
and `SoilType`. All these variables are present in the `SpatRaster`
(`var`) we imported. Therefore, we can now predict our models to this
raster. Let’s begin by predicting the maxnet model:

``` r
p_maxnet <- predict_selected(models = fitted_model_maxnet, 
                             new_variables = var,
                             progress_bar = FALSE)
```

By default, the function computes consensus metrics (mean, median,
range, and standard deviation) for each model across its replicates (if
more than one model was selected), as well as a general consensus across
all models. In this case, the output is a `list` containing `SpatRaster`
predictions for each replicate, along with the consensus results for
each model and the overall general consensus:

``` r
#See objects in the output of predict_selected
names(p_maxnet)
#> [1] "Model_192"         "Model_219"         "General_consensus"
```

Let’s plot the general consensus:

``` r
plot(p_maxnet$General_consensus)
```

![](model_predictions_files/figure-html/plot%20maxnet%20general-1.png)

We can also plot the results for each replicate and the consensus for
each model:

``` r
#Predictions for each replicate from model 192
plot(p_maxnet$Model_192$Replicates)
```

![](model_predictions_files/figure-html/plot%20models%20maxnet-1.png)

``` r

#Consensus across each replicate from model 192
plot(p_maxnet$Model_192$Model_consensus)
```

![](model_predictions_files/figure-html/plot%20models%20maxnet-2.png)

  

For comparison, let’s predict the GLM model:

``` r
# Predict glm model
p_glm <- predict_selected(models = fitted_model_glm, 
                          new_variables = var,
                          progress_bar = FALSE)
#See selected models that were predicted
names(p_glm)
#> [1] "Model_85"          "General_consensus"

#Compare general consensus (mean) between maxnet and glm
par(mfrow= c(1, 2)) #Set grid to plot
plot(p_maxnet$General_consensus$mean, main = "Maxnet")
plot(p_glm$General_consensus, main = "GLM")
```

![](model_predictions_files/figure-html/unnamed-chunk-2-1.png)

``` r
on.exit() #Reinitiate grid
```

  

### Predict to data.frame

Instead of a `SpatRaster`, we can also predict the models to a
`data.frame` that stores the variable values. To see an example, let’s
convert the raster variables `var` to a `data.frame`:

``` r
var_df <- as.data.frame(var)
head(var_df)
#>       bio_1    bio_7 bio_12   bio_15 SoilType
#> 11 22.77717 18.12400   1180 48.03594       NA
#> 12 22.76711 17.74400   1191 49.31194       10
#> 13 22.68580 17.46575   1206 51.51922       10
#> 14 22.50121 17.84525   1228 53.90265       10
#> 15 22.07609 18.14125   1254 54.10397       10
#> 16 21.88485 18.80800   1276 54.07279       10
```

Note that each column stores the values for each variable. Let’s predict
our Maxnet models to this `data.frame`:

``` r
p_df <- predict_selected(models = fitted_model_maxnet, 
                         new_variables = var_df, #Now, a data.frame
                         progress_bar = FALSE) 
```

Now, instead of `SpatRaster` objects, the function returns `data.frame`
objects with the predictions:

``` r
#Results by replicate of the model 192
head(p_df$Model_192$Replicates)
#>   Partition_1  Partition_2  Partition_3  Partition_4   Full_model
#> 1 0.006521501 0.0006209852 0.0005883615 9.831561e-05 9.526648e-08
#> 2 0.006446437 0.0005356316 0.0005713501 9.009486e-05 8.993584e-08
#> 3 0.006233583 0.0003396879 0.0004975279 6.967025e-05 8.559934e-08
#> 4 0.005797668 0.0001458500 0.0003605775 4.303576e-05 8.306910e-08
#> 5 0.008513515 0.0002034105 0.0006532983 8.529550e-05 4.220225e-07
#> 6 0.009240381 0.0001753492 0.0006784553 9.035171e-05 6.540127e-07

#Consensus across replicates of the model 192
head(p_df$Model_192$Model_consensus)
#>         median       range        mean       stdev
#> 1 0.0005883615 0.006521406 0.001565852 0.002784420
#> 2 0.0005356316 0.006446347 0.001528721 0.002761027
#> 3 0.0003396879 0.006233497 0.001428111 0.002693874
#> 4 0.0001458500 0.005797585 0.001269443 0.002535186
#> 5 0.0002034105 0.008513093 0.001891188 0.003710540
#> 6 0.0001753492 0.009239727 0.002037038 0.004035351

#General consensus across all models
head(p_df$General_consensus)
#>         median       range        mean       stdev
#> 1 0.0005892047 0.006521444 0.001506370 0.002502583
#> 2 0.0005348599 0.006446378 0.001477821 0.002497636
#> 3 0.0003390597 0.006233519 0.001390363 0.002461362
#> 4 0.0001459093 0.005797595 0.001250202 0.002350678
#> 5 0.0002032855 0.008513101 0.001870336 0.003456477
#> 6 0.0001758720 0.009239727 0.002031620 0.003796839
```

  

### Binarize Models

The `fitted_models` object stores the thresholds that can be used to
binarize the models into suitable and unsuitable areas. These thresholds
correspond to the omission error rate used during model selection (e.g.,
5% or 10%).

You can access the omission error rate used to calculate the thresholds
directly from the object:

``` r
#Get omission error used to select models and calculate the thesholds
## For maxnet model
fitted_model_maxnet$omission_rate
#> [1] 10

## For glm model
fitted_model_glm$omission_rate
#> [1] 10
```

In both models, a 10% omission error rate was used to calculate the
thresholds. This means that when predictions are binarized,
approximately 10% of the presence records used in model calibration will
fall into areas classified as unsuitable.

The thresholds are summarized in two ways: the mean and median across
replicates for each model, and the consensus mean and median across all
selected models (when more than one model is selected). Let’s check the
thresholds for the general consensus:

``` r
#For maxnet
fitted_model_maxnet$thresholds$consensus
#> $mean
#> [1] 0.3095083
#> 
#> $median
#> [1] 0.259534

#For glm
fitted_model_glm$thresholds$consensus
#> $mean
#> [1] 0.1204713
#> 
#> $median
#> [1] 0.1204713
```

Let’s use these thresholds to binarize the models (this functionality is
only available when predicting to a `SpatRaster`):

``` r
#Get the thersholds for models (general consensus)
thr_mean_maxnet <- fitted_model_maxnet$thresholds$consensus$mean #Maxnet
thr_mean_glm <- fitted_model_glm$thresholds$consensus$mean #glm

#Binarize models
mean_maxnet_bin <- (p_maxnet$General_consensus$mean > thr_mean_maxnet) * 1
mean_glm_bin <- (p_glm$General_consensus > thr_mean_glm) * 1

#Compare results
par(mfrow= c(1, 2)) #Set grid to plot
plot(mean_maxnet_bin, main = "Maxnet")
plot(mean_glm_bin, main = "GLM")
```

![](model_predictions_files/figure-html/binarize%20models-1.png)

``` r
on.exit() #Reinitiate grid
```

  

### Clamping Variables

By default, predictions are performed with free extrapolation
(`extrapolation_type = "E"`). This can be problematic when the peak of
suitability occurs at the extremes of a predictor’s range. For example,
let’s examine the response curve of the Maxnet model for `bio_7`
(Temperature Annual Range):

``` r
response_curve(models = fitted_model_maxnet, variable = "bio_7", 
               extrapolation_factor = 1)
```

![](model_predictions_files/figure-html/response%20bio_7-1.png)

Note that higher suitability occurs at low values of the temperature
range. However, the lower limit of the calibration data used to fit the
models (dashed line) is at 15.7ºC. The premise that suitability will
increase and stabilize at lower values of `bio_7` is an extrapolation of
the model (the area to the left of the dashed line). It’s possible that
suitability could begin to decrease at extremely low values, rendering
this extrapolation inaccurate, but the calibration data is insufficient
for the model to predict this.

One way to address this is by clamping the variables. This means that
all values outside the calibration range (both below the lower value and
above the upper value) are set to the respective lower and upper limits
of the calibration range. For example, in the calibration data for the
Maxnet models, the lower and upper limits for `bio_7` are 15.7ºC and
23.3ºC, respectively:

``` r
range(fitted_model_maxnet$calibration_data$bio_7)
#> [1] 15.71120 23.30475
```

To observe the effect of clamping this variable, let’s create a
hypothetical (and extreme) scenario where `bio_7` has extremely low
values:

``` r
#From bio_7, reduce values
new_bio7 <- var$bio_7 - 3
#Create new scenario
new_var <- var
#Replace bio_7 with new_bio7 in this scenario
new_var$bio_7 <- new_bio7

#Plot the differences
par(mfrow = c(1,2))
plot(var$bio_7, main = "Original bio_7", range = c(5, 25))
plot(new_var$bio_7, main = "New bio_7", range = c(5, 25))
```

![](model_predictions_files/figure-html/create%20scenario-1.png)

``` r
on.exit() #Reinitiate grid
```

Let’s predict the Maxnet models for this new scenario with both free
extrapolation (`extrapolation_type = "E"`) and with clamped variables
(`extrapolation_type = "EC"`):

``` r
#Predict to hypothetical scenario with free extrapolation
p_free_extrapolation <- predict_selected(models = fitted_model_maxnet, 
                                         new_variables = new_var, #New scenario
                                         consensus = "mean",
                                         extrapolation_type = "E", #Free extrapolation (Default)
                                         progress_bar = FALSE)

#Predict to hypothetical scenario with clamping
p_clamping <- predict_selected(models = fitted_model_maxnet, 
                               new_variables = new_var, #New scenario
                               consensus = "mean",
                               extrapolation_type = "EC", #Extrapolation with clamping
                               progress_bar = FALSE)

#Get and see differences
p_difference <- p_free_extrapolation$General_consensus$mean - p_clamping$General_consensus$mean

#Plot the differences
par(mfrow = c(2,2))
plot(p_free_extrapolation$General_consensus$mean, main = "Free extrapolation",
     zlim = c(0, 1))
plot(p_clamping$General_consensus$mean, main = "Clamping",
     zlim = c(0, 1))
plot(p_difference, main = "Difference")
plot(new_bio7, main = "Hypothetical bio_7", type = "interval")
```

![](model_predictions_files/figure-html/unnamed-chunk-3-1.png)

``` r
on.exit() #Reinitiate grid
```

Note that when we clamp the variables, regions with extremely low values
of (the hypothetical) `bio_7` exhibit lower predicted suitabilities
compared to when free extrapolation is allowed.

By default, when `extrapolation_type = "EC"` is set, all predictor
variables are clamped. You can specify which variables to clamp using
the `var_to_clamp` argument.

  

### No Extrapolation

A more rigorous approach is to **predict with no extrapolation**, where
regions outside the limits of the calibration data are assigned a
suitability value of 0. Let’s predict the Maxnet models using the
hypothetical scenario we created in the previous step to observe the
difference:

``` r
#Predict to hypothetical scenario with no extrapolation
p_no_extrapolation <- predict_selected(models = fitted_model_maxnet, 
                                       new_variables = new_var, #New scenario
                                       consensus = "mean",
                                       extrapolation_type = "NE", #No extrapolation
                                       progress_bar = FALSE)
#Plot the differences
par(mfrow = c(2,2))
plot(p_free_extrapolation$General_consensus$mean, main = "Free extrapolation",
     zlim = c(0, 1))
plot(p_clamping$General_consensus$mean, main = "Clamping",
     zlim = c(0, 1))
plot(p_no_extrapolation$General_consensus$mean, main = "No extrapolation",
     zlim = c(0, 1))
plot(new_bio7, main = "Hypothetical bio_7", type = "interval")
```

![](model_predictions_files/figure-html/no%20extrapolation-1.png)

``` r
on.exit() #Reinitiate grid
```

n this example, a large portion of the predicted area shows zero
suitability. This is because, in this hypothetical scenario, much of the
region has `bio_7` values lower than those in the calibration data,
which has a minimum of 15ºC. Suitability values greater than zero are
only predicted in areas where `bio_7` falls within the range of the
calibration data.

By default, when `extrapolation_type = "NE"` is set, all predictor
variables are considered for this process. You can specify a subset of
variables to be considered for extrapolation using the `var_to_clamp`
argument.

  

### Output Type

Maximum entropy models (maxnet) produce four different types of output
for their predictions: raw, cumulative, logistic, and cloglog. These are
described in [Merow et
al. 2013](https://nsojournals.onlinelibrary.wiley.com/doi/10.1111/j.1600-0587.2013.07872.x)
and [Phillips et
al. 2017](https://nsojournals.onlinelibrary.wiley.com/doi/10.1111/ecog.03049).

All four output types are monotonically related. Therefore, rank-based
metrics for model fit (e.g., omission rate and partial ROC) will be
identical. However, the output types have different scaling, which leads
to distinct interpretations and visually different prediction maps.

- **Raw (or exponential) output** is interpreted as a Relative
  Occurrence Rate (ROR). The ROR sums to 1 when predicted to the
  calibration data.
- **Cumulative output** assigns to a location the sum of all raw values
  less than or equal to the raw value for that location, and then
  rescales this to range between 0 and 100. Cumulative output can be
  interpreted in terms of an omission rate because thresholding at a
  value of *c* to predict a suitable/unsuitable cell will omit
  approximately *c*% of presences.
- **Cloglog output (Default)** transforms the raw values into a scale of
  relative suitability ranging between 0 and 1, using a logistic
  transformation based on a user-specified parameter ‘$\tau$’, which
  represents the probability of presence at ‘average’ presence
  locations. In this context, the tau value defaults to
  $\tau \approx 0.632$.
- **Logistic output** is similar to Cloglog, but it assumes that
  $\tau = 0.5$.

Let’s examine the differences between these four output types for Maxnet
models:

``` r
p_cloglog <- predict_selected(models = fitted_model_maxnet, new_variables = var, 
                              type = "cloglog", progress_bar = FALSE)
p_logistic <- predict_selected(models = fitted_model_maxnet, new_variables = var, 
                              type = "logistic", progress_bar = FALSE)
p_cumulative <- predict_selected(models = fitted_model_maxnet, new_variables = var, 
                              type = "cumulative", progress_bar = FALSE)
p_raw <- predict_selected(models = fitted_model_maxnet, new_variables = var, 
                              type = "raw", progress_bar = FALSE)

#Plot the differences
par(mfrow = c(2,2))
plot(p_cloglog$General_consensus$mean, main = "Cloglog (Default)",
     zlim = c(0, 1))
plot(p_logistic$General_consensus$mean, main = "Logistic",
     zlim = c(0, 1))
plot(p_cumulative$General_consensus$mean, main = "Cumulative",
     zlim = c(0, 1))
plot(p_raw$General_consensus$mean, main = "Raw",
     zlim = c(0, 1))
```

![](model_predictions_files/figure-html/output%20types-1.png)

``` r
on.exit() #Reinitiate grid
```

  

## Saving Predictions

We can save the predictions to the disk by setting `write_files = TRUE`.
When this option is enabled, you must provide a directory path in the
`out_dir` argument.

If `new_variables` is a `SpatRaster`, the function will save the output
files as GeoTIFF (.tif) files. If `new_variables` is a `data.frame`, the
function will save the output files as Comma Separated Value (.csv)
files.

``` r
p_save <- predict_selected(models = fitted_model_maxnet, 
                           new_variables = var, 
                           write_files = TRUE, #To save to the disk
                           write_replicates = TRUE, #To save predictions for each replicate
                           out_dir = tempdir(), #A path to save the resuls (here, the temporary directory)
                           progress_bar = FALSE)
```

Alternatively, we can use
[`writeRaster()`](https://rspatial.github.io/terra/reference/writeRaster.html)
to save specific output predictions manually. For example, to save only
the mean layer from the general consensus results:

``` r
writeRaster(p_maxnet$General_consensus$mean, 
            filename = file.path(tempdir(), "Mean_consensus.tif"))
```
