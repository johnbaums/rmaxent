<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Travis-CI Build Status](https://travis-ci.org/johnbaums/rmaxent.svg?branch=master)](https://travis-ci.org/johnbaums/rmaxent)

rmaxent: working with Maxent Species Distribution Models in R
=============================================================

(*Please note that some characters and equations do not render correctly in Github Flavored Markdown. Please refer to the [html version](https://rawgit.com/johnbaums/rmaxent/master/README.html) of this document.*)

Correlative species distribution models (SDMs; Franklin & Miller 2010) are now the most common tool for predicting habitat suitability. Maxent, a machine-learning regression-type approach to fitting SDMs based on the principle of maximum entropy (Phillips *et al.* 2004, 2006; Elith *et al.* 2011), is used in a vast proportion of published SDM studies. The Maxent software is written in Java, and provides a graphical user interface in addition to command line operation. In 2010, the `dismo` R package (Hijmans *et al.* 2016) was added to CRAN, providing, amongst other features, an R interface to Maxent that streamlined the process of preparing data, and fitting, evaluating, and projecting models.

Additional functionality is provided by the `rmaxent` package, which allows Java-free projection of previously-fitted Maxent models, and provides several other convenience functions. The core function of the package is `project`, which builds upon a previous description of the relationship between covariate (i.e., "feature") values and Maxent's fitted values (Wilson 2009). In my test, projection with the `project` function is at least twice as fast as maxent.jar (the standard Java implementation), and there is scope for further gains by taking advantage of C++ libraries (e.g., via the `Rcpp` package---planned for future releases). These speed gains are of particular use when projecting numerous Maxent models, such as when exploring sensitivity of suitability surfaces to model settings, or when projecting models to numerous environmental scenarios, as is increasingly common when considering potential climate change.

The `rmaxent` package also includes function `ic`, which calculates information criteria (AIC, AIC<sub>c</sub>, BIC) for Maxent models as implemented in ENMTools (Warren *et al.* 2010). These quantities can be used to optimise model complexity (e.g., Warren & Seifert 2011), and for highlighting model parsimony. The user should note, though, that this approach uses the number of parameters (Maxent features with non-zero weights) in place of degrees of freedom when calculating model likelihood, and this may underestimate the true degrees of freedom, particularly when hinge and/or threshold features are in use (see Warren *et al.* 2014 for details). However, despite this potential issue, model selection based on this calculation of AIC<sub>c</sub> has been shown to outperform selection based on predictive capacity (i.e., using AUC; Warren & Seifert 2011).

Finally, `rmaxent` also provides functions to:

-   import raster data stored in Maxent’s binary .mxe raster format (`read_mxe`; written in collaboration with Peter D. Wilson);
-   parse Maxent .lambdas files (files that contain information about model features), returning information about feature types, weights, minima and maxima, as well as the model’s entropy and other estimated constants (`parse_lambdas`); and
-   create limiting factor maps (Elith *et al.* 2010) that identify the environmental variable that is least favourable at each point across the landscape (`limiting`).

Installation
------------

We can install the `rmaxent` package from GitHub, using the `devtools` package:

``` r
library(devtools)
install_github('johnbaums/rmaxent')
```

``` r
library(rmaxent)
```

Examples
--------

Projecting a fitted Maxent model requires predictor states for all variables included in the model, and the model's ".lambdas" file---a plain text file containing details of all features considered by the model, including their weights (i.e., coefficients), minima, maxima, and some constants required in the calculation of fitted values.

Below, we use the example data distributed with the `dismo` package. These data include coordinates representing localitions where the brown-throated three-toed sloth, *Bradypus variegatus* has been recorded, and spatial, gridded data giving biome classification and values for a range of current climate variables.

Let's import the *B. variegatus* occurrence and predictor data from the appropriate paths:

``` r
occ_file <- system.file('ex/bradypus.csv', package='dismo')
occ <- read.table(occ_file, header=TRUE, sep=',')[,-1]

library(raster)
pred_files <- list.files(system.file('ex', package='dismo'), '\\.grd$', full.names=TRUE )
predictors <- stack(pred_files)
```

The object `predictors` is a `RasterStack` comprising nine raster layers, one for each predictor used in the model.

We can now fit the model using the `maxent` function from the `dismo` package. Note that this function calls Maxent's Java program, `maxent.jar`. Our objective here is to fit a model in order to demonstrate the functionality of `rmaxent`. For the sake of the exercise we will disable hinge and threshold features.

``` r
library(dismo)
me <- maxent(predictors, occ, factors='biome', args=c('hinge=false', 'threshold=false'))
```

The Maxent model has now been fit, and the resulting object, `me`, which is of class `MaxEnt`, can be passed to various functions in `rmaxent`. For example, `project` takes a trained Maxent model and predicts it to new data. The procedure for calculating fitted values from a Maxent .lambdas file and a vector of predictor values for a given site is as follows:

1.  clamp each untransformed predictor to its training extrema (i.e., the maximum and minimum of the model-fitting data), by setting all values greater than the maximum to maximum, and all values less than the minimum to the minimum;
2.  considering only non-linear features with non-zero weights (see description of `parse_lambdas`), take each and calculate its value. For example, if a quadratic feature has a non-zero weight, the quadratic feature's value is the square of the corresponding linear feature;
3.  clamp each non-hinge feature to its training extrema, as in step 1;
4.  normalise all features so that their values span the range [0, 1]. Maxent's procedure for this depends on the feature type. For each feature \(x_j\), the corresponding normalised feature \(x_j^\ast\) is calculated as

\[
  \begin{equation} \label{eq:normfeat}
    x_j^\ast=
    \begin{cases}
      \frac{\text{max}x_j - x_j}{\text{max}x_j - \text{min}x_j},  & \text{if }x_j\text{ is a reverse hinge feature}\\%[1em]
      \frac{x_j - \text{min}x_j}{\text{max}x_j - \text{min}x_j},  & \text{otherwise}
    \end{cases} 
  \end{equation}
  \]

1.  calculate \(X^\ast\cdot\beta\), the dot product of the vector of normalised feature values, and the corresponding vector of feature weights;
2.  calculate \(y_{\text{raw}}\) Maxent's "raw" output by subtracting a normalising constant from \(X^\ast\cdot\beta\), exponentiating the result, and dividing by a second normalising constant (these constants are, respectively, the `linearPredictorNormalizer` and `densityNormalizer` returned by `parse_lambdas`); and finally,
3.  calculate Maxent's "logistic" output (often interpreted as habitat suitability, \(HS\)) as follows, where \(H\) is the model entropy (returned by `parse_lambdas`)

\[
\begin{equation} \label{eq:maxentlogistic}
  HS = 1 - \frac{1}{e^H y_{\text{raw}} + 1}.
\end{equation}
\]

Using this procedure, we predict the model to the model-fitting data below:

``` r
prediction <- project(me, predictors)
```

And plot the result:

``` r
library(rasterVis)
library(viridis)
levelplot(prediction$prediction_logistic, margin=FALSE, col.regions=viridis, at=seq(0, 1, len=100)) +
  layer(sp.points(SpatialPoints(occ), pch=20, col=1))
```

<img src="README_files/figure-markdown_github/plot1-1.png" alt="__Figure 1. Maxent habitat suitability prediction for the brown-throated three-toed sloth, _Bradypus variegatus_.__"  />
<p class="caption">
**Figure 1. Maxent habitat suitability prediction for the brown-throated three-toed sloth, *Bradypus variegatus*.**
</p>

We can compare the time taken to project the model to the model-fitting landscape with `project`, versus using the typical `predict.MaxEnt` method shipped with `dismo`.

``` r
library(microbenchmark)
timings <- microbenchmark(
  rmaxent=pred_rmaxent <- project(me, predictors),
  dismo=pred_dismo <- predict(me, predictors), 
  times=10)
```

``` r
print(timings, signif=2)
```

    ## Unit: milliseconds
    ##     expr min  lq     mean median  uq max neval cld
    ##  rmaxent 100 110 141.1127    120 200 200    10  a 
    ##    dismo 270 280 358.5439    320 370 560    10   b

On average, the `dismo` method takes approximately 2.5 times as long as the `rmaxent` method. Here the difference is rather trivial, but when projecting to data with higher spatial resolution and/or larger extent, the gains in efficiency are welcome, particularly if projecting many models to multiple environmental scenarios.

We can check that the predictions are equivalent, at least to machine precision:

``` r
all.equal(values(pred_rmaxent$prediction_logistic), values(pred_dismo))
```

    ## [1] "Mean relative difference: 2.143141e-08"

It is useful to know that `project` returns a list containing Maxent's raw output as well as its logistic output. The raw output can be accessed with `pred_rmaxent$prediction_raw`, and is required for calculating model information criteria, as we will see when demonstrating the use of `ic`, below.

Once a model has been projected, information about the features used in the model can be extracted from the fitted model object, or the .lambdas file, with `parse_lambdas`. For example,

``` r
parse_lambdas(me)
```

    ## $lambdas
    ##          feature        var       lambda    min     max        type
    ## 1   (biome==1.0)      biome   1.49226955      0       1 categorical
    ## 2   (biome==2.0)      biome   1.15527404      0       1 categorical
    ## 3   (biome==9.0)      biome   2.33122618      0       1 categorical
    ## 4  (biome==13.0)      biome   1.96179711      0       1 categorical
    ## 5  (biome==14.0)      biome   0.32250218      0       1 categorical
    ## 6           bio1       bio1   6.60855618    -23     289      linear
    ## 7          bio12      bio12   0.00000000      0    7682      linear
    ## 8          bio16      bio16   0.32171392      0    2458      linear
    ## 9          bio17      bio17  -4.08708925      0    1496      linear
    ## 10          bio5       bio5   0.00000000     61     422      linear
    ## 11          bio6       bio6   0.00000000   -212     240      linear
    ## 12          bio7       bio7 -16.18361978     62     461      linear
    ## 13          bio8       bio8   1.85855298    -66     323      linear
    ## 14        bio5^2       bio5  -2.16295546   3721  178084   quadratic
    ## 15        bio6^2       bio6  -5.26414653      0   57600   quadratic
    ## 16        bio7^2       bio7  -7.49172108   3844  212521   quadratic
    ## 17        bio8^2       bio8   0.12388710      0  104329   quadratic
    ## 18    bio12*bio7 bio12,bio7   4.15357004      0  737464     product
    ## 19    bio12*bio8 bio12,bio8   2.44714442 -27324 2020366     product
    ## 20    bio16*bio8 bio16,bio8   0.01359303 -13332  638038     product
    ## 21    bio17*bio7 bio17,bio7   0.07878229      0  145705     product
    ## 22     bio5*bio7  bio5,bio7  -3.27012487   8235  162812     product
    ## 
    ## $linearPredictorNormalizer
    ## [1] 5.526369
    ## 
    ## $densityNormalizer
    ## [1] 175.0504
    ## 
    ## $numBackgroundPoints
    ## [1] 9766
    ## 
    ## $entropy
    ## [1] 7.932286

This information can be useful, since it shows how many, and which, features have non-zero weights. However, the function is perhaps more useful in its role as a helper function for other functions in the package. For example, the values returned by `parse_lambdas` are required for calculating fitted values (shown above).

To identify which variable is most responsible for decreasing suitability in a given environment, we can use the `limiting` function. This is an R implementation of an approach described previously and now incorporated into Maxent. The limiting variable at a given location is identified by calculating the decrease in suitability, for each predictor in turn, relative to the suitability (logistic prediction) that would be achieved if that predictor took the value equal to the mean at occurrence sites (median for categorical variables). The predictor associated with the largest decrease in suitability is the most limiting factor.

``` r
lim <- limiting(predictors, me)
levelplot(lim, col.regions=rainbow) +
  layer(sp.points(SpatialPoints(occ), pch=20, col=1))
```

<img src="README_files/figure-markdown_github/limiting-1.png" alt="__Figure 2. The variable that most limits the suitability of habitat for the brown-throated three-toed sloth, _Bradypus variegatus_. Black points indicate occurrence localities.__"  />
<p class="caption">
**Figure 2. The variable that most limits the suitability of habitat for the brown-throated three-toed sloth, *Bradypus variegatus*. Black points indicate occurrence localities.**
</p>

Figure 2 shows that for much of the Americas, the BIOCLIM variable 7 (annual temperature range) is most limiting for *B. variegatus*.

Finally, we can calculate information criteria describing the balance of complexity and fit of the Maxent model. The calculation of these criteria follows that of Warren *et al*. (2010). In the context of Maxent, it has been suggested that likelihood may not be calculated correctly by this approach, since the number of parameters may not be equal to the number of features with non-zero weights (Hastie *et al.* 2009). This may lead to underparameterised models , but relative to other common approaches to SDM selection (e.g., AUC), AIC<sub>c</sub>-based model selection, in particular, has been shown to lead to models with improved transferability, accuracy, and ecological relevance (Warren & Seifert 2011).

Information criteria are typically used as relative measures of model support, thus we will fit and project a second model for comparison to the existing model. The new model will have higher beta-regularisation, permitting a smoother fit to the training data that may be less prone to being locally overfit, but is otherwise identical to the first model.

``` r
me2 <- maxent(predictors, occ, factors='biome', args=c('hinge=false', 'threshold=false', 'betamultiplier=5'))
pred2 <- project(me2, predictors)
```

We can now calculate and compare these quantities using `ic`.

``` r
ic(stack(pred_rmaxent$prediction_raw, pred2$prediction_raw), 
   occ, list(me, me2))
```

    ##          n  k        ll      AIC     AICc      BIC
    ## layer.1 94 19 -736.2354 1510.471 1520.741 1558.793
    ## layer.2 94 10 -748.3635 1516.727 1519.378 1542.160

We see above that AIC<sub>c</sub>, which converges to AIC as \(n\) gets large (Burnham & Anderson 2004), is marginally lower for the simpler model. Difference in the two models could be interrogated further by comparing their results for `parse_lambdas`, and by examining response curves in the standard Maxent output.

References
----------

Burnham, K.P. & Anderson, D.R. (2004). Multimodel inference: Understanding AIC and BIC in model selection. *Sociological Methods & Research*, **33**, 261–304.

Elith, J., Kearney, M. & Phillips, S. (2010). The art of modelling range-shifting species. *Methods in Ecology and Evolution*, **1**, 330–342.

Elith, J., Phillips, S.J., Hastie, T., Dudík, M., Chee, Y.E. & Yates, C.J. (2011). A statistical explanation of MaxEnt for ecologists. *Diversity and Distributions*, **17**, 43–57.

Franklin, J. & Miller, J.A. (2010). *Mapping Species Distributions*. Cambridge University Press, New York.

Hastie, T.J., Tibshirani, R.J. & Friedman, J.H. (2009). *The Elements of Statistical Learning: Data Mining, Inference, and Prediction*. Springer.

Hijmans, R.J., Phillips, S., Leathwick, J. & Elith, J. (2016). *dismo: Species Distribution Modeling*.

Phillips, S., Anderson, R. & Schapire, R. (2006). Maximum entropy modeling of species geographic distributions. *Ecological Modelling*, **190**, 231–259.

Phillips, S., Dudík, M. & Schapire, R. (2004). A maximum entropy approach to species distribution modeling. *Proceedings of the Twenty-First International Conference on Machine Learning*.

Warren, D.L. & Seifert, S.N. (2011). Ecological niche modeling in Maxent: The importance of model complexity and the performance of model selection criteria. *Ecological Applications*, **21**, 335–342.

Warren, D.L., Glor, R.E. & Turelli, M. (2010). ENMTools: A toolbox for comparative studies of environmental niche models. *Ecography*, **33**, 607–611.

Warren, D.L., Wright, A.N., Seifert, S.N. & Shaffer, H.B. (2014). Incorporating model complexity and spatial sampling bias into ecological niche models of climate change risks faced by 90 California vertebrate species of concern. *Diversity and distributions*, **20**, 334–343.

Wilson, P.D. (2009). *Guidelines for computing MaxEnt model output values from a lambdas file*.
