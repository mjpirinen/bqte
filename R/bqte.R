# These functions compute
# Quantile Treatment Effect (QTE) and
# Back-transformed Quantile Treatment Effect (BQTE).
# BQTE is a transformation of QTE as a function of the outcome value in control group.
# Uses bootstrap for estimating uncertainty.

# Matti Pirinen & Harri Hemila, version Aug 19, 2025.

# The source code in this file is published under
# the MIT license listed at the end of this file.

# Contains 9 functions:

# bqte(Treatment, Control, at = NULL,
#      bqte.conf = 0.95, B = 2000,
#      K = length(Control), bagging = TRUE,
#      tails = FALSE,
#      discrete.range = FALSE,
#      exact.sample.quantiles = TRUE,
#      interpolation.method = "linear",
#      spline.df = 5,
#      verbose = TRUE)

# plot_bqte(x, plot.ci = TRUE, plot.rbqte = FALSE,
#           plot.range = FALSE,
#           pch = 19, col = "black",
#           col.ci = "red",
#           col.range = "gray",
#           xlim = NULL, ylim = NULL,
#           xlab = NULL, ylab = NULL,
#           xaxs = "r", yaxs = "r",
#           xaxp = NULL, yaxp = NULL,
#           cex = 1, cex.axis = 1, cex.lab = 1,
#           main ="", lwd = 1)

# plot_tbqte(x, utbqte = TRUE, plot.ci = TRUE,
#            plot.relative = FALSE, pch = NULL, col = "black",
#            col.ci = "gray", xlim = NULL, ylim = NULL,
#            xlab = NULL, ylab = NULL, xaxs = "r", yaxs = "r",
#            xaxp = NULL, yaxp = NULL,
#            cex = 1, cex.axis = 1, cex.lab = 1,
#            main ="", lwd = 1)
#
# plot_tbqte_both(
#     x,
#     x.utbqte = NULL, x.ltbqte = NULL,
#     plot.ci = TRUE, plot.relative = FALSE,
#     pch.utbqte = 25, pch.ltbqte = 24,
#     col.utbqte = "blue",
#     col.ltbqte = "purple",
#     col.utbqte.ci = "blue",
#     col.ltbqte.ci = "purple",
#     jitter = c(0,0),
#     xlim = NULL, ylim = NULL,
#     xlab = NULL, ylab = NULL,
#     xaxs = "r", yaxs = "r",
#     xaxp = NULL, yaxp = NULL,
#     cex = 1, cex.axis = 1, cex.lab = 1,
#     main ="", lwd = 1)

# qte(Treatment, Control, at = NULL,
#     qte.conf = 0.95, B = 2000,
#     bagging = TRUE,
#     exact.sample.quantiles = TRUE,
#     verbose = TRUE)

# plot_qte(x, plot.ci = TRUE, plot.rqte = FALSE,
#          pch = 19, col = "black",
#          col.ci = "gray",
#          xlim = NULL, ylim = NULL,
#          xlab = NULL, ylab = NULL,
#          xaxs = "r", yaxs = "r",
#          xaxp = NULL, yaxp = NULL,
#          lwd = 1)

# quantiles_from_discrete(p, pr)

# exact_bqte_discrete(at, x, pr.con, pr.trt, verbose = TRUE)

# bqte_doksum(at, Treatment, Control)

# For examples of usage, see file 'bqte_examples.R'.


################################################################################
############## FUNCTIONS #######################################################
################################################################################


#' Estimate backtransformed quantile treatment effects
#'
#' We have observed survival/recovery times, or other numerical outcome values,
#' in the treatment group and in the control group.
#' Quantile treatment effect (QTE) is the difference between
#' the two groups as a function of quantile level in (0,1).
#' We want to transform QTE to the original scale of the outcome value,
#' at a given grid of values ('at') in the control group.
#' These transformed values are called
#' backtransformed quantile treatment effects' (BQTE).
#'
#' Approach:
#' Fix a grid of outcome values in the control group at which
#' BQTEs are computed (parameter 'at').
#' Recommendation:
#'   Do not try to estimate BQTEs at values that are outside
#'   empirical quantile levels (10/N, 1 - (10/N)),
#'   where N is the minimum sample size of the treatment and control groups
#'   because accuracy of quantile estimates is low at the tails.
#'
#' Estimate the quantile treatment effects as difference between the quantiles
#' of Treatment group and Control group at 'K' quantile levels
#' 1/(K+1),...,K/(K+1).
#' Make a piece-wise linear approximation of BQTE
#' as a function of outcome value in controls,
#' based on the observed K quantiles.
#' If several quantiles are equal, then use average of the corresponding
#' Treatment group values as the value of BQTE at that point.
#' (See parameter 'discrete.range' to get more details of range of tied values.)
#' Use that approximation
#' to estimate the BQTE at the chosen grid of outcome values in controls.
#' NOTE: values are returned to the user on the grid defined by 'at',
#'       with a sensible default range, and not all K values included in the computation
#'       are meant to be returned to the user.
#'
#' Uses bootstrapping to estimate the uncertainty of BQTE at each grid point,
#' by resampling with replacement a set of Treatment and Control group values that
#' have the same sample size as the original Treatment and Control groups, respectively,
#' and applying the procedure explained above to each bootstrap sample.
#' Note that bootstrapping accounts for uncertainty of quantiles of both groups.
#'
#' The final estimate of BQTE at a grid point is the mean over bootstrap samples ("bagging estimate")
#' and the confidence interval for BQTE is estimated from the quantiles of the bootstrap sample.
#' NOTE: In empirical evaluation it was observed that bagging estimate tend to have smaller
#'       mean square error in discrete distributions than the direct estimate from the data and
#'       therefore the bagging estimate is used by default.
#'       This can be switched by argument "bagging".
#'
#' For discrete data, one can also compute the range of Treatment values that correspond
#' to a single Control group value and return that range (See argument 'discrete.range'.)
#' This range is typically wider than the confidence interval of the BQTE estimate
#' as BQTE is defined as the expected value over those Treatment group values that
#' correspond to the quantile levels of a particular Control roup value.
#'
#' Can also compute upper tail BQTE (UTBQTE) and lower tail BQTE (LTBQTE) at each point in 'at'.
#' UTBQTE, at outcome value 'x' is the difference in the mean outcome values
#'  of the highest 1-F(x) of the two populations, where F is the CDF of the Control population.
#' LTBQTE, at outcome value 'x' is the difference in the mean outcome values
#'  of the lowest F(x) of the two populations, where F is the CDF of Control population.
#' UTBQTE(x) is an upper bound for the average treatment effect of that part of
#' the Control population who have outcome value >= x.
#' LTBQTE(x) is a lower bound for the average treatment effect of that part of
#' the Control population who have outcome value <= x.
#' These bounds are general in the sense that they do not require an assumption
#' about how the treatment operates, for example, they do not assume anything about
#' whether the treatment preserves order.
#' Note that the outcome value 'x' is included in the definition of both UTBQTE(x)
#' and in LTBQTE(x).
#' The confidence intervals of UTBQTE and LTBQTE are estimated using bootstrap
#' and the point estimates are computed either using bagging (if bagging = TRUE)
#' or without bagging (if bagging = FALSE).
#'
#' @param Treatment vector of outcome values in Treatment group.
#' @param Control vector of outcome values in Control group.
#' @param at vector of outcome values in controls at which BQTE is estimated.
#'     By default, between 10 and 20 points between the empirical quantile levels
#'     of 10/N and 1 - (10/N) in controls where N = min(|Treatment|, |Control|).
#' @param bqte.conf (default 0.95), confidence level for BQTE.
#' @param B (default 2000), number of bootstrap samples.
#' @param K (default 'length(Control)), number of quantiles where BQTE function is estimated.
#' @param bagging (default TRUE), if TRUE, reports bagging estimate (mean over bootstrap samples)
#'                           otherwise reports the direct estimate from the data.
#' @param tails (default FALSE), if TRUE, returns also estimates and confidence intervals for
#'                         UTBQTE and LTBQTE.
#' @param discrete.range (default FALSE), if TRUE, returns, for every Control group value,
#'                                the range of Treatment group values that correspond
#'                                through the shared quantile levels.
#' @param exact.sample.quantiles (default TRUE), if TRUE, uses empirical sample quantiles,
#'                         that are called type = 1 in R's quantile function,
#'                         if FALSE, uses R's default quantiles (type = 7).
#' @param interpolation.method either "linear" for piecewise-linear model or "spline"
#'                       for natural splines of df 'spline.df'.
#' @param spline.df (default 5) degrees of freedom of natural spline.
#'                        Relevant only if interpolation.method = "spline".
#' @param verbose (default TRUE) if TRUE, prints parameters and range recommendations on console.
#'
#' @return list with 1,2, or 3 data frames
#'
#' $res
#'
#' 'at', outcome value in control group;
#'
#' 'bqte', backtransformed quantile treatment effect (averaged in case of ties);
#'
#' 'bqte.low', lower CI end point for 'bqte';
#'
#' 'bqte.up', upper CI end point for 'bqte';
#'
#' 'rbqte', relative backtransformed quantile treatment effect with respect to control outcome value;
#'
#' 'rbqte.low', lower CI end point for 'rbqte';
#'
#' 'rbqte.up', upper CI end point for 'rbqte';
#'
#' $discrete.range, returned only if argument discrete.range = TRUE;
#'
#' 'at', outcome value in control group;
#'
#' 'range.low', lower CI end point of raw, non-averaged BQTE values corresponding to this Control group value;
#'
#' 'range.up' upper CI end point;
#'
#' 'rrange.low' lower CI end point of relative BQTE;
#'
#' 'rrange.up' upper CI end point of relative BQTE;
#'
#' $tails, returned only if parameter 'tails = TRUE'
#'
#' 'utbqte, upper tail average treatment effect estimate;
#'
#' 'utbqte.low, lower CI end point for 'utbqte';
#'
#' 'utbqte.up, upper CI end point for 'utbqte';
#'
#' 'ltbqte, lower tail average treatment effect estimate;
#'
#' 'ltbqte.low, lower CI end point for 'ltbqte';
#'
#' 'ltbqte.up, upper CI end point for 'ltbqte';
#'
#' 'rutbqte, relative upper tail average treatment effect estimate with respect to control value;
#'
#' 'rutbqte.low, lower CI end point for 'rutbqte';
#'
#' 'rutbqte.up, upper CI end point for 'rutbqte';
#'
#' 'rltbqte, relative lower tail average treatment effect estimate with respect to control value;
#'
#' 'rltbqte.low, lower CI end point for 'rltbqte';
#'
#' 'rltbqte.up, upper CI end point for 'rltbqte'
#' @examples
#' set.seed(12); Tr = rweibull(100,1.2,1); Co = rnorm(100,3,5)
#' bqte(Tr, Co)
#' @export
bqte <- function(Treatment, Control, at = NULL,
                 bqte.conf = 0.95, B = 2000,
                 K = length(Control), bagging = TRUE,
                 tails = FALSE,
                 discrete.range = FALSE,
                 exact.sample.quantiles = TRUE,
                 interpolation.method = "linear",
                 spline.df = 5,
                 verbose = TRUE){

  if( K < 2 ) stop("K should be at least 2.")
  if( bqte.conf <= 0 | bqte.conf >= 1 ) stop("Value \'bqte.conf\' not valid.")
  if( B < 1 ) stop("B is not positive integer.")
  if(!(interpolation.method %in% c("linear","spline"))) stop("interpolation.method must be either 'linear' or 'spline'")
  C.n = length(Control)
  T.n = length(Treatment)
  min.n = min(c(C.n, T.n))

  if(exact.sample.quantiles) q.type = 1 else q.type = 7
  #recommended to keep value in 'at' within the range 'x':
  x = as.numeric(quantile(Control, prob = c(10/min.n, 1 - 10/min.n), type = q.type))
  if(is.null(at)){ # Set 'at' to its default value
    # that has between 10 and 20 values from interval between empirical quantiles
    # of Control vector of levels 10/C.n and 1 - 10/C.n.
    # The step size is set to the appropriate power of 2.
    d = 2^(floor(log2((x[2] - x[1])/10))) # step size
    at = d*seq(ceiling(x[1]/d), floor(x[2]/d), 1)
    # Check the endpoints for rounding errors and add possibly missing endpoints:
    tolerance = 1e-10 #tolerance for rounding error
    if(at[1] - d >= x[1] - tolerance){
      at = c(at[1] - d, at)} # add left endpoint
    if(at[length(at)] + d <= x[2] + tolerance) {
      at = c(at, at[length(at)] + d)} # add right end point
  }

  if( verbose ){
    cat("Running bqte() with the following parameters.\n")
    cat(paste("length(Treatment):", T.n),"\n")
    cat(paste("length(Control):", C.n),"\n")
    cat(paste("K:", K),"\n")
    cat(paste("B:", B),"\n")
    cat(paste("bqte.conf:", bqte.conf),"\n")
    cat(paste("at:", paste(signif(at,3),collapse=", ")),"\n")
    cat(paste("bagging:", bagging),"\n")
    cat(paste("tails:", tails),"\n")
    cat(paste("exact.sample.quantiles:",exact.sample.quantiles),"\n")
    cat(paste("discrete.range:",discrete.range),"\n")
    cat(paste("interpolation.method:",interpolation.method),"\n")
    cat(paste0("If you choose to change 'at' values, a recommended interval is [",
               signif(x[1],3),", ",signif(x[2],3),"].\n"))
  }
  probs = (1:K)/(K+1) # quantile levels at which the BQTE-function is estimated
  bqte.b = matrix(NA, nrow = B, ncol = length(at)) # bootstrapped  BQTEs at grid 'at'
  low.range.b = matrix(NA, nrow = B, ncol = length(at)) # bootstrapped low end of discrete values at grid 'at'
  up.range.b = matrix(NA, nrow = B, ncol = length(at)) # bootstrapped upper end of discrete values at grid 'at'
  utbqte.b = matrix(NA, nrow = B, ncol = length(at)) # bootstrapped UTbqtes at grid 'at'
  rutbqte.b = matrix(NA, nrow = B, ncol = length(at)) # bootstrapped relative UTbqtes at grid 'at'
  ltbqte.b = matrix(NA, nrow = B, ncol = length(at)) # bootstrapped LTbqtes at grid 'at'
  rltbqte.b = matrix(NA, nrow = B, ncol = length(at)) # bootstrapped relative LTbqtes at grid 'at'

  for (ii in 1:B){

    #1. Make bootstrap sample
    C.b = sample(Control, size = C.n, replace = TRUE)
    T.b = sample(Treatment, size = T.n, replace = TRUE)

    #2. Pick the fixed set of quantiles from the bootstrap sample
    C.q = as.numeric(quantile(C.b, prob = probs, type = q.type))
    T.q = as.numeric(quantile(T.b, prob = probs, type = q.type))

    #3. Estimate BQTE at grid 'at' using linear approximation
    # NOTE: In case of ties in C.q, the value is taken to be the mean.
    #       Thus we estimate an _average_ difference for each grid point,
    #       and do not consider the full uncertainty of values
    #       corresponding to each point but only uncertainty related to the mean!
    yy = T.q - C.q
    if(interpolation.method == "linear"){
      bqte.b[ii,] = as.numeric(approx(C.q, yy, xout = at, rule = 2, ties = mean)$y)}
    else{
      bqte.b[ii,] = as.numeric(predict(lm(yy ~ splines::ns(C.q, df = spline.df)),
                                       newdata = data.frame(C.q = at)))
    }

    #4. Estimate lower and upper end points of ranges corresponding to discrete values
    #  (if discrete.range = TRUE)
    if(discrete.range){
       ranges = lapply(at, function(x){
         kk = which(C.q == x)
         if(length(kk) > 0){return(range(T.q[kk]-x))}
         else{return(c(NA,NA))}
         })
      A = matrix(unlist(ranges), byrow = TRUE, ncol = 2)
      low.range.b [ii,] = as.numeric(A[,1])
      up.range.b [ii,] = as.numeric(A[,2])
    }

    #5. Estimate utbqte and ltbqte at grid 'at' using linear approximation (if 'tails = TRUE')
    if(tails){

      ut.C.mean.q = sapply(1:K, function(x){mean(C.q[as.integer(x):K])})
      ut.T.mean.q = sapply(1:K, function(x){mean(T.q[as.integer(x):K])})

      #UTBQTE is used as an upper bound for the average treatment effect of the upper tail.
      # Thus, we take the whole right tail including the target value.

      utbqte.q = ut.T.mean.q  - ut.C.mean.q
      utbqte.b[ii,] = approx(C.q, utbqte.q, xout = at, rule = 2, ties = function(x){x[1]})$y
      rutbqte.q = (ut.T.mean.q  - ut.C.mean.q)/ut.C.mean.q
      rutbqte.b[ii,] = approx(C.q, rutbqte.q, xout = at, rule = 2, ties = function(x){x[1]})$y

      lt.C.mean.q = sapply(1:K, function(x){mean(C.q[1:as.integer(x)])})
      lt.T.mean.q = sapply(1:K, function(x){mean(T.q[1:as.integer(x)])})

      #LTBQTE is used as a lower bound for the average treatment effect of the
      # lower tail. Thus, we take the whole left tail including the target value.

      ltbqte.q = lt.T.mean.q  - lt.C.mean.q
      ltbqte.b[ii,] = approx(C.q, ltbqte.q, xout = at, rule = 2, ties = function(x){x[length(x)]})$y
      rltbqte.q = (lt.T.mean.q  - lt.C.mean.q)/lt.C.mean.q
      rltbqte.b[ii,] = approx(C.q, rltbqte.q, xout = at, rule = 2, ties = function(x){x[length(x)]})$y
    }
  }# ends bootstrap

  bqte.low = as.numeric(apply(bqte.b, 2,
                    function(X){quantile(X, probs = (1 - bqte.conf)/2, na.rm = TRUE)}))
  bqte.up = as.numeric(apply(bqte.b, 2,
                    function(X){quantile(X, probs = (1 + bqte.conf)/2, na.rm = TRUE)}))
  if(discrete.range){
    range.low = as.numeric(apply(low.range.b, 2,
                    function(X){quantile(X, probs = (1 - bqte.conf)/2, na.rm = TRUE)}))
    range.up = as.numeric(apply(up.range.b, 2,
                    function(X){quantile(X, probs = (1 + bqte.conf)/2, na.rm = TRUE)}))
  }

  if(tails){
    ltbqte.low = as.numeric(apply(ltbqte.b, 2, function(X){quantile(X, (1 - bqte.conf)/2)}))
    ltbqte.up = as.numeric(apply(ltbqte.b, 2, function(X){quantile(X, (1 + bqte.conf)/2)}))
    utbqte.low = as.numeric(apply(utbqte.b, 2, function(X){quantile(X, (1 - bqte.conf)/2)}))
    utbqte.up = as.numeric(apply(utbqte.b, 2, function(X){quantile(X, (1 + bqte.conf)/2)}))
    rltbqte.low = as.numeric(apply(rltbqte.b, 2, function(X){quantile(X, (1 - bqte.conf)/2)}))
    rltbqte.up = as.numeric(apply(rltbqte.b, 2, function(X){quantile(X, (1 + bqte.conf)/2)}))
    rutbqte.low = as.numeric(apply(rutbqte.b, 2, function(X){quantile(X, (1 - bqte.conf)/2)}))
    rutbqte.up = as.numeric(apply(rutbqte.b, 2, function(X){quantile(X, (1 + bqte.conf)/2)}))
  }

  if(bagging) { # Mean as the final estimate
    bqte = apply(bqte.b, 2, mean)
    if(tails){
      utbqte = apply(utbqte.b, 2, mean)
      ltbqte = apply(ltbqte.b, 2, mean)
      rutbqte = apply(rutbqte.b, 2, mean)
      rltbqte = apply(rltbqte.b, 2, mean)
    }
  }
  else{ # Compute estimate directly from the data (not bagging)
    C.q = as.numeric(quantile(Control, prob = probs, type = q.type))
    T.q = as.numeric(quantile(Treatment, prob = probs, type = q.type))
    bqte = as.numeric(approx(C.q, T.q - C.q, xout = at, rule = 2, ties = mean)$y)
    if(tails){
      ut.T.mean.q = sapply(1:K, function(x){mean(T.q[as.integer(x):K])})
      ut.C.mean.q = sapply(1:K, function(x){mean(C.q[as.integer(x):K])})
      utbqte.q = ut.T.mean.q  - ut.C.mean.q

      utbqte = as.numeric(approx(C.q, utbqte.q, xout = at, rule = 2, ties = function(x){x[1]})$y)
      rutbqte.q = ut.T.mean.q/ut.C.mean.q - 1
      rutbqte = as.numeric(approx(C.q, utbqte.q, xout = at, rule = 2, ties = function(x){x[1]})$y)

      lt.T.mean.q = sapply(1:K, function(x){mean(T.q[1:as.integer(x)])})
      lt.C.mean.q = sapply(1:K, function(x){mean(C.q[1:as.integer(x)])})
      ltbqte.q = lt.T.mean.q  - lt.C.mean.q
      ltbqte = as.numeric(approx(C.q, ltbqte.q, xout = at, rule = 2, ties = function(x){x[length(x)]})$y)
      rltbqte.q = lt.T.mean.q/lt.C.mean.q - 1
      rltbqte = as.numeric(approx(C.q, ltbqte.q, xout = at, rule = 2, ties = function(x){x[length(x)]})$y)
    }
  }

  res = data.frame(
    at = at,
    bqte = bqte,
    bqte.low = bqte.low,
    bqte.up = bqte.up,
    rbqte = bqte/at,
    rbqte.low = bqte.low/at,
    rbqte.up = bqte.up/at)

  res.range = NULL
  if(discrete.range){
    res.range = data.frame(
      at = at,
      range.low = range.low,
      range.up = range.up,
      rrange.low = range.low/at,
      rrange.up = range.up/at
    )
  }

  res.tails = NULL
  if(tails){
    res.tails = data.frame(
      at = at,
      utbqte = utbqte,
      utbqte.low = utbqte.low,
      utbqte.up = utbqte.up,
      rutbqte = rutbqte,
      rutbqte.low = rutbqte.low,
      rutbqte.up = rutbqte.up,
      ltbqte = ltbqte,
      ltbqte.low = ltbqte.low,
      ltbqte.up = ltbqte.up,
      rltbqte = rltbqte,
      rltbqte.low = rltbqte.low,
      rltbqte.up = rltbqte.up)}

  return(list(res = res, range = res.range, tails = res.tails))
}

#' Plots BQTE
#'
#' Plots backtransformed quantile treatment effect (BQTE) on direct or relative scale.
#'
#' @param x data.frame returned by function bqte( )
#' @param plot.ci if TRUE, plots confidence intervals around estimates using color col.ci
#' @param plot.rbqte if TRUE, plots relative BQTEs, otherwise plots actual BQTEs
#' @param plot.range if TRUE, plots bqte$range results (requires discrete data)
#' @param pch pointstyle, if NULL, then lines are drawn instead of separate points
#' @param col color
#' @param col.ci color of confidence interval
#' @param col.range color of range
#' @param xlim range of x-coordinates
#' @param ylim range of y-coordinates
#' @param xlab label of x-axis
#' @param ylab label of y-axis
#' @param xaxs as R's standard plotting parameter
#' @param yaxs as R's standard plotting parameter
#' @param xaxp as R's standard plotting parameter
#' @param yaxp as R's standard plotting parameter
#' @param cex as R's standard plotting parameter
#' @param cex.axis as R's standard plotting parameter
#' @param cex.lab as R's standard plotting parameter
#' @param main as R's standard plotting parameter
#' @param lwd as R's standard plotting parameter
#' @return none
#' @examples
#' set.seed(12); Tr = rweibull(30,1.2,1); Co = rnorm(30,3,5)
#' plot_bqte(bqte(Tr, Co))
#' @export
plot_bqte <- function(x, plot.ci = TRUE,
                      plot.rbqte = FALSE,
                      plot.range = FALSE,
                      pch = 19, col = "black",
                      col.ci = "red",
                      col.range = "gray",
                      xlim = NULL, ylim = NULL,
                      xlab = NULL, ylab = NULL,
                      xaxs = "r", yaxs = "r",
                      xaxp = NULL, yaxp = NULL,
                      cex = 1, cex.axis = 1, cex.lab = 1,
                      main ="", lwd = 1){

  at = x$res$at
  if(plot.range & is.null(x$range)) stop("Input bqte object misses 'range'. Rerun bqte with 'discrete.range = TRUE'.")
  #Use percentages for relative effects
  if(plot.rbqte) {
    x$res[,c("rbqte","rbqte.low","rbqte.up")] = 100*x$res[,c("rbqte","rbqte.low","rbqte.up")]
    y = x$res[,"rbqte"]; y.low = x$res[,"rbqte.low"]; y.up = x$res[,"rbqte.up"]
    if(plot.range) {
      y.low.range = 100*x$range[,"rrange.low"]
      y.up.range = 100*x$range[,"rrange.up"]
    }
  }else{
    y = x$res[,"bqte"]; y.low = x$res[,"bqte.low"]; y.up = x$res[,"bqte.up"]
    if(plot.range) {
      y.low.range = x$range[,"range.low"]
      y.up.range = x$range[,"range.up"]
    }
  }

  if(is.null(xlim)) {
    xlim = range(at) + c(-1,1) * 0.05 * diff(range(at))}
  if(is.null(ylim)) {
    ran = range(c(y.low, y.up))
    if(plot.range) {
      ran = range(c(ran, y.low.range, y.up.range))}
    ylim = ran + c(-1,1) * 0.05 * diff(ran)}
  if(is.null(xlab)) xlab = "Outcome in controls"
  if(is.null(ylab)){
    if(!plot.rbqte) ylab = "BQTE"
    else ylab = "Relative BQTE (%)"
  }

  plot(NULL,
       xlim = xlim, ylim = ylim,
       xlab = xlab, ylab = ylab,
       xaxs = xaxs, yaxs = yaxs,
       xaxp = xaxp, yaxp = yaxp,
       main = main, cex.lab = cex.lab,
       cex.axis = cex.axis)
  grid()
  abline(h = 0, lty = 3)

  #add intervals
  if(plot.range){
    arrows(at, y.low.range,
           at, y.up.range,
           col = col.range, code = 3, angle = 90, length = 0.0, lty = 2, lwd = lwd)}

  if(plot.ci){
    if(!is.null(pch)) {
      arrows(at, y.low,
             at, y.up,
             col = col.ci, code = 3, angle = 90, length = 0.0, lwd = lwd*1.5)}
    else{
      polygon(c(at,rev(at)),
              c(y.low, rev(y.up)),
              col = col.ci, border = NA)
    }
  }

  #add points
  if(!is.null(pch)) {
    points(at, y, pch = pch, col = col, cex = cex)}
  else {
    lines(at, y, lty = 1, col = col, lwd = lwd)}
}

#' Plot tail BQTE
#'
#' Plot either lower (LTBQTE) or upper (UTBQTE) tail BQTE on direct or relative scale.
#' For definitions of these measures, see help of bqte().
#' @param x data.frame returned by function bqte( ,tails = TRUE).
#' @param utbqte if TRUE, prints UTBQTE else prints LTBQTE.
#' @param plot.ci if TRUE, plots confidence intervals around estimates using color col.ci.
#' @param plot.relative if TRUE, plots relative TBQTEs, otherwise plots direct TBQTEs.
#' @param pch pointstyle if NULL, then lines are drawn instead of separate points.
#' @param col color.
#' @param col.ci color of confidence interval.
#' @param xlim range of x-coordinates.
#' @param ylim range of y-coordinates.
#' @param xlab label of x-axis.
#' @param ylab label of y-axis.
#' @param xaxs as R's standard plotting parameter.
#' @param yaxs as R's standard plotting parameter.
#' @param xaxp as R's standard plotting parameter.
#' @param yaxp as R's standard plotting parameter.
#' @param cex as R's standard plotting parameter.
#' @param cex.axis as R's standard plotting parameter.
#' @param cex.lab as R's standard plotting parameter.
#' @param main as R's standard plotting parameter.
#' @param lwd as R's standard plotting parameter.
#' @return none
#' @examples
#' set.seed(12); Tr = rweibull(30,1.2,1); Co = rnorm(30,3,5)
#' plot_tbqte(bqte(Tr, Co, tails = TRUE), utbqte = FALSE, pch = 19)
#' @export
plot_tbqte <- function(x, utbqte = TRUE, plot.ci = TRUE,
                       plot.relative = FALSE,
                      pch = NULL, col = "black",
                      col.ci = "gray",
                      xlim = NULL, ylim = NULL,
                      xlab = NULL, ylab = NULL,
                      xaxs = "r", yaxs = "r",
                      xaxp = NULL, yaxp = NULL,
                      cex = 1, cex.axis = 1, cex.lab = 1,
                      main ="", lwd = 1){

  at = x$res$at
  if(utbqte) {
    if(plot.relative){
      X = x$tails[,c("rutbqte","rutbqte.low","rutbqte.up")]}
    else{
      X = x$tails[,c("utbqte","utbqte.low","utbqte.up")]}
  }else {
    if(plot.relative){
      X = x$tails[,c("rltbqte","rltbqte.low","rltbqte.up")]}
    else{
      X = x$tails[,c("ltbqte","ltbqte.low","ltbqte.up")]}
  }

  #Use percentages for relative effects
  if(plot.relative) X = 100*X

  if(is.null(xlim)) {
    xlim = range(at) + c(-1,1) * 0.05 * diff(range(at))}
  if(is.null(ylim)) {
      ylim = range(c(X[,2], X[,3])) + c(-1,1) * 0.05 * diff(range(c(X[,2], X[,3])))}

  if(is.null(xlab)) xlab = "Outcome in controls"
  if(is.null(ylab)){
    lab = c("LTBQTE","UTBQTE")[1 + as.numeric(utbqte)]
    if(!plot.relative) ylab = lab
    else ylab = paste("Relative",lab,"(%)")
  }

  plot(NULL,
       xlim = xlim, ylim = ylim,
       xlab = xlab, ylab = ylab,
       xaxs = xaxs, yaxs = yaxs,
       xaxp = xaxp, yaxp = yaxp,
       main = main, cex.lab = cex.lab,
       cex.axis = cex.axis)
  grid()
  abline(h = 0, lty = 3, lwd = lwd)

  #add intervals
  if(plot.ci){
    if(!is.null(pch)) {
      arrows(at, X[,2],
             at, X[,3],
             col = col.ci, code = 3, angle = 90, length = 0.0, lwd = lwd)}
    else{
      polygon(c(at,rev(at)),
              c(X[,2], rev(X[,3])),
              col = col.ci, border = NA)
    }
  }

  #add points
  if(!is.null(pch)) {
    points(at, X[,1], pch = pch, col = col, cex = cex)}
  else{
    lines(at, X[,1], col = col, lty = 1, lwd = lwd)}
}

#' Plot simultaneously both tail BQTEs
#'
#' Plot in the same figure both UTBQTE and LTBQTE on direct or relative scale.
#' For definitions of these measures, see help of bqte()
#'
#' @param x data.frame returned by function bqte( ,tails = TRUE).
#' @param x.utbqte control group's outcome value from which utbqte is plotted towards the upper tail.
#'                  If NULL, then the lowest point of data in 'x' is used.
#' @param x.ltbqte control group's outcome value from which ltbqte is plotted towards the lower tail.
#'                   If NULL, then highest point of data in 'x' is used.
#' @param plot.relative If TRUE, plots relative TBQTEs, otherwise plots direct TBQTEs.
#' @param plot.ci  If TRUE, plots confidence intervals around estimates using colors.
#'                col.utbqte.ci and col.ltbqte.ci
#' @param jitter vector of length 2, gives the horizontal offsets for points to avoid overlaps.
#'                 between plotted LTBQTE and UTBQTE estimates and intervals.
#'                 1st value is offset for UTBQTE, 2nd value is offset for LTBQTE
#'                 default is (0,0), i.e., no offset.
#' @param pch.utbqte pointstyle for utbqte. If NULL then the corresponding estimates are shown by a line
#  and the confidence intervals' endpoints are connected by lines.
#' @param pch.ltbqte pointstyle for utbqte. If NULL then the corresponding estimates are shown by a line
#  and the confidence intervals' endpoints are connected by lines.
#' @param col.utbqte color of utbqte.
#' @param col.ltbqte color of ltbqte.
#' @param col.utbqte.ci color of confidence interval for utbqte.
#' @param col.ltbqte.ci color of confidence interval for ltbqte.
#' @param xlim range of x-coordinates.
#' @param ylim range of y-coordinates.
#' @param xlab label of x-axis.
#' @param ylab label of y-axis.
#' @param xaxs as R's standard plotting parameter.
#' @param yaxs as R's standard plotting parameter.
#' @param xaxp as R's standard plotting parameter.
#' @param yaxp as R's standard plotting parameter.
#' @param cex as R's standard plotting parameter.
#' @param cex.axis as R's standard plotting parameter.
#' @param cex.lab as R's standard plotting parameter.
#' @param main as R's standard plotting parameter.
#' @param lwd as R's standard plotting parameter.
#' @return none
#' @examples
#' set.seed(12); Tr = rweibull(30,1.2,1); Co = rnorm(30,3,5)
#' plot_tbqte_both(bqte(Tr, Co, tails = TRUE))
#' @export
plot_tbqte_both <- function(
    x,
    x.utbqte = NULL, x.ltbqte = NULL,
    plot.ci = TRUE, plot.relative = FALSE,
    pch.utbqte = 25, pch.ltbqte = 24,
    col.utbqte = "blue",
    col.ltbqte = "purple",
    col.utbqte.ci = "blue",
    col.ltbqte.ci = "purple",
    jitter = c(0,0),
    xlim = NULL, ylim = NULL,
    xlab = NULL, ylab = NULL,
    xaxs = "r", yaxs = "r",
    xaxp = NULL, yaxp = NULL,
    cex = 1, cex.axis = 1, cex.lab = 1,
    main ="", lwd = 1){

  at = x$res$at
  if(plot.relative){
      X = x$tails[,c("rutbqte","rutbqte.low","rutbqte.up","rltbqte","rltbqte.low","rltbqte.up")]}
    else{
      X = x$tails[,c("utbqte","utbqte.low","utbqte.up","ltbqte","ltbqte.low","ltbqte.up")]}

  #Use percentages for relative effects
  if(plot.relative) X = 100*X

  if(is.null(xlim)) {
    xlim = range(at) + c(-1,1) * 0.05 * diff(range(at))}
  if(is.null(ylim)) {
    ylim = range(as.vector(X[,c(2,3,5,6)])) + c(-1,1) * 0.05 * diff(range(as.vector(X[,c(2,3,5,6)])))}

  if(is.null(xlab)) xlab = "Outcome in controls"
  if(is.null(ylab)){
    lab = c("L/UTBQTE")
    if(!plot.relative) ylab = lab
    else ylab = paste("Relative",lab,"(%)")
  }

  plot(NULL,
       xlim = xlim, ylim = ylim,
       xlab = xlab, ylab = ylab,
       xaxs = xaxs, yaxs = yaxs,
       xaxp = xaxp, yaxp = yaxp,
       main = main, cex.lab = cex.lab,
       cex.axis = cex.axis)
  grid()
  abline(h = 0, lty = 3, lwd = lwd)

  if(is.null(x.utbqte)) x.utbqte = min(at)
  if(is.null(x.ltbqte)) x.ltbqte = max(at)
  u.ind = (at >= x.utbqte)
  l.ind = (at <= x.ltbqte)

  #add intervals
  if(plot.ci){
    if(!is.null(pch.utbqte)) {
      arrows(at[u.ind]+jitter[1], X[u.ind,2],
             at[u.ind]+jitter[1], X[u.ind,3],
             col = col.utbqte.ci, code = 3, angle = 90, length = 0.0, lwd = lwd)}
    else{
      polygon(c(at[u.ind],rev(at[u.ind])),
              c(X[u.ind,2], rev(X[u.ind,3])),
              col = col.utbqte.ci, border = NA)
    }
    if(!is.null(pch.ltbqte)) {
      arrows(at[l.ind]+jitter[2], X[l.ind,5],
             at[l.ind]+jitter[2], X[l.ind,6],
             col = col.ltbqte.ci, code = 3, angle = 90, length = 0.0, lwd = lwd)}
    else{
      polygon(c(at[l.ind],rev(at[l.ind])),
              c(X[l.ind,5], rev(X[l.ind,6])),
              col = col.ltbqte.ci, border = NA)
    }
  }

  #add points
  if(!is.null(pch.utbqte)) {
    points(at[u.ind]+jitter[1], X[u.ind,1], cex = cex,
           pch = pch.utbqte, col = col.utbqte, bg = col.utbqte)}
  else{
    lines(at[u.ind], X[u.ind,1], col = col.utbqte, lty = 1, lwd = lwd)}

  if(!is.null(pch.ltbqte)) {
    points(at[l.ind]+jitter[2], X[l.ind,4], cex = cex,
           pch = pch.ltbqte, col = col.ltbqte, bg = col.ltbqte)}
  else{
    lines(at[l.ind], X[l.ind,4], col = col.ltbqte, lty = 1, lwd = lwd)}

}


#' Estimate quantile treatment effects.
#'
#' We have observed outcome values, such as survival or recovery times,
#' in the treatment group and in the control group.
#' Quantile treatment effect (QTE) is the difference between the two groups
#' as a function of the quantile levels.
#'
#' Fix a grid of quantile levels where QTEs are computed (parameter 'at').
#'  Recommendation:
#'   Do not try to estimate QTEs at values that are outside
#'   quantile levels (10/N, 1 - (10/N)),
#'   where N is the minimum of sample sizes of treated and controls because of low accuracy at the tails.
#'
#' Estimate the quantile treatment effects as difference between the quantiles
#'  of Treatment group and Control group
#'
#' Use bootstrapping to estimate the uncertainty of QTE,
#'  by resampling with replacement a set of Treatment and Control group values that
#'  have the same sample size as the original Treatment and Control groups, respectively,
#'  and by applying the procedure explained above to each bootstrap sample.
#'  Note that bootstrapping accounts for uncertainty of quantiles of both groups.
#'
#' The final estimate of QTE is the mean over bootstrap samples ("bagging estimate")
#'  and the confidence interval for QTE is estimated from the quantiles of the bootstrap sample.
#' NOTE: In empirical evaluation it was observed that bagging estimate tend to have smaller
#'       mean square error in discrete distributions than the direct estimate from the data and
#'      therefore the bagging estimate is used by default. This can be switched by parameter "bagging".
#' @param Treatment vector of outcome values (e.g. survival times) in Treatment group.
#' @param Control vector of outcome values in Control group.
#' @param at vector of quantiles at which QTE is estimated.
#'       By default, between 10 and 20 points between the quantile levels
#'       of 10/N and 1 - (10/N) where N = min(c(length(Control),length(Treatment))).
#' @param qte.conf (default 0.95), confidence level for QTE estimate.
#' @param B (default 2000), number of bootstrap samples.
#' @param bagging (default TRUE), if TRUE, reports bagging estimate (mean over bootstrap samples)
#'                           otherwise reports the direct estimate from the data.
#' @param exact.sample.quantiles (default TRUE), if TRUE, uses empirical sample quantiles (type = 1)
#'                             if FALSE, uses R's default quantiles (type = 7), where 'type' is
#'                              argument of R's quantile() function.
#' @param verbose (default TRUE) if TRUE, prints parameters and recommendations on console.
#' @return data frame with 7 columns
#'
#' at, quantile level;
#'
#' qte, quantile treatment effect;
#'
#' qte.low, lower CI end point for 'qte';
#'
#' qte.up, upper CI end point for 'qte';
#'
#' rqte, relative quantile treatment effect with respect to control outcome value;
#'
#' rqte.low, lower CI end point for 'rqte';
#'
#' rqte.up, upper CI end point for 'rqte';
#' @examples
#' set.seed(12); Tr = rweibull(100,1.2,1); Co = rnorm(100,3,5)
#' qte(Tr, Co)
#' @export
qte <- function(Treatment, Control, at = NULL,
                qte.conf = 0.95, B = 2000,
                bagging = TRUE,
                exact.sample.quantiles = TRUE,
                verbose = TRUE){

  if( qte.conf <= 0 | qte.conf >= 1 ) stop("Value \'qte.conf\' not valid.")
  if( B < 1 ) stop("B is not positive integer.")
  C.n = length(Control)
  T.n = length(Treatment)

  #recommended to keep value in 'at' within the range 'x'
  min.n = min(c(C.n,T.n))
  if(min.n < 5) stop("Minimum allowed sample size is 5 in both groups.")
  if( is.null(at) ){ # Set 'at' to its default value
    # that has at most 19 values from interval between quantile levels 10/min.n and 1 - 10/min.n.
    K = max(which(1/(2:20) >= 10/min.n))
    at = (1:K)/(K+1)
  }
  at = sort(at)
  if(any(at <= 0 | at >= 1)) stop("Values in 'at' should be in interval (0,1).")

  if( verbose ){
    cat("Running qte() with the following parameters.\n")
    cat(paste("length(Treatment):", T.n),"\n")
    cat(paste("length(Control):", C.n),"\n")
    cat(paste("B:", B),"\n")
    cat(paste("qte.conf:", qte.conf),"\n")
    cat(paste("at:", paste(signif(at, 3),collapse=", ")),"\n")
    cat(paste("bagging:", bagging),"\n")
    x = c(10/min.n, 1 - 10/min.n)
    cat(paste0("If you choose to change 'at' values, a recommended interval is [",
               signif(x[1],3),", ",signif(x[2],3),"].\n"))
  }

  if(exact.sample.quantiles) q.type = 1 else q.type = 7

  C.b = matrix(NA, nrow = B, ncol = length(at)) # bootstrapped control quantiles
  T.b = matrix(NA, nrow = B, ncol = length(at)) # bootstrapped treatment quantiles
  QTE.b = matrix(NA, nrow = B, ncol = length(at)) # bootstrapped approximations of QTEs at grid 'at'
  for (ii in 1:B){
    C.b[ii,] = as.numeric(quantile(sample(Control, size = C.n, replace = T), prob = at, type = q.type))
    T.b[ii,] = as.numeric(quantile(sample(Treatment, size = T.n, replace = T), prob = at, type = q.type))
    QTE.b[ii,] = as.numeric(T.b[ii,] - C.b[ii,])
  }

  qte.low = as.numeric(apply(QTE.b, 2, function(X){quantile(X, (1 - qte.conf)/2)}))
  qte.up = as.numeric(apply(QTE.b, 2, function(X){quantile(X, (1 + qte.conf)/2)}))

  if(bagging) { # Mean as the final estimate
    C.q = apply(C.b, 2, mean)
    qte = apply(QTE.b, 2, mean)
  }
  else{ #Compute QTE estimate directly from the data (not bagging)
    C.q = as.numeric(quantile(Control, prob = at, type = q.type))
    T.q = as.numeric(quantile(Treatment, prob = at, type = q.type))
    qte = as.numeric(T.q - C.q)
  }

  res = data.frame(
    at = at,
    qte = qte,
    qte.low = qte.low,
    qte.up = qte.up,
    rqte = qte/C.q,
    rqte.low = qte.low/C.q,
    rqte.up = qte.up/C.q)

  return(res)
}

#' Plot quantile treatment effect
#'
#' Plot quantile treatment effect (QTE) on direct or relative scale
#'
#' @param x data.frame returned by the qte () function.
#' @param plot.ci if TRUE, plots confidence intervals around estimates.
#' @param plot.rqte if TRUE, plots relative QTEs, otherwise plots actual QTEs.
#' @param pch pointstyle, if NULL, then lines are drawn instead of separate points
#' @param col color
#' @param col.ci color of confidence interval
#' @param xlim range of x-coordinates
#' @param ylim range of y-coordinates
#' @param xlab label of x-axis
#' @param ylab label of y-axis
#' @param xaxs as R's standard plotting parameter
#' @param yaxs as R's standard plotting parameter
#' @param xaxp as R's standard plotting parameter
#' @param yaxp as R's standard plotting parameter
#' @param lwd as R's standard plotting parameter
#' @return none
#' @examples
#' set.seed(12); Tr = rweibull(30,1.2,1); Co = rnorm(30,3,5)
#' plot_qte(qte(Tr, Co))
#' @export
plot_qte <- function(x, plot.ci = TRUE, plot.rqte = FALSE,
                     pch = 19, col = "black",
                     col.ci = "gray",
                     xlim = NULL, ylim = NULL,
                     xlab = NULL, ylab = NULL,
                     xaxs = "r", yaxs = "r",
                     xaxp = NULL, yaxp = NULL,
                     lwd = 1){

  #Use percentages for relative effects
  if(plot.rqte) x[,c("rqte","rqte.low","rqte.up")] = 100*x[,c("rqte","rqte.low","rqte.up")]

  if(is.null(xlim)) {
    xlim = c(0,1)}
  if(is.null(ylim)) {
    if(!plot.rqte){
      ylim = range(c(x$qte.low, x$qte.up)) + c(-1,1) * 0.05 * diff(range(c(x$qte.low, x$qte.up)))}
    else{
      ylim = range(c(x$rqte.low, x$rqte.up)) + c(-1,1) * 0.05 * diff(range(c(x$rqte.low, x$rqte.up)))}
  }
  if(is.null(xlab)) xlab = "Quantile level"
  if(is.null(ylab)){
    if(!plot.rqte) ylab = "QTE"
    else ylab = "Relative QTE (%)"
  }

  plot(NULL,
       xlim = xlim, ylim = ylim,
       xlab = xlab, ylab = ylab,
       xaxs = xaxs, yaxs = yaxs,
       xaxp = xaxp, yaxp = yaxp)
  grid()
  abline(h = 0, lty = 2, lwd = lwd)

  #add intervals
  if(plot.ci){
    if(!plot.rqte){y.low = x$qte.low; y.up = x$qte.up}
    else{y.low = x$rqte.low; y.up = x$rqte.up}
    arrows(x$at, y.low,
           x$at, y.up,
           col = col.ci, code = 3, angle = 90, length = 0.0, lwd = lwd*1.5)
  }

  #add points
  if(!plot.rqte) {
    points(x$at, x$qte, pch = pch, col = col)}
  else {
    points(x$at, x$rqte, pch = pch, col = col)}
}

#' Doksum estimator
#'
#' Estimator of BQTE published by Doksum (1974), Annals of Statistics Vol.2 No.2 267-277.
#' @param at Control values where estimation is done.
#' @param Treatment Observed values of the treatment group
#' @param Control Observed values of the control group
#' @return vector of estimates at the values given in 'at'.
#' @examples
#' set.seed(12); Tr = rweibull(30,1.2,1); Co = rnorm(30,3,5)
#' bqte_doksum(at = c(3,5), Tr, Co)
#' @export
bqte_doksum <- function(at, Treatment, Control){
  y = sapply(at, function(x){floor(mean(Control <= x) * length(Treatment)) + 1})
  y[y > length(Treatment)] = length(Treatment)
  sort(Treatment)[y] - at
}

#' Quantile from a discrete distribution
#'
#' @param p quantile levels.
#' @param pr probability mass vector of non-negative values.
#' @return list with two vectors: (1) quantile.level and (2) indexes.
#'         'quantile level' equals input vector 'p' and
#'         indexes contain the indexes of 'pr' corresponding to
#'         the quantile levels in 'p'.
#' @examples
#' quantiles_from_discrete(p = c(0.3,0.5), pr = c(0.2,0.2,0.3,0.3))
#' @export
quantiles_from_discrete <- function(p, pr){
  if(any(pr < 0)) stop("'pr' must have non-negative values")
  pr = pr/sum(pr)
  y = sapply(p,function(x){min(which(cumsum(pr) >= x))})
  return(list(quantile.level = p, index = y))
}

#' BQTE for discrete distribution
#'
#' @param at value at which bqte is computed.
#' @param x outcome values corresponding to probabilities in
#'          pr.con and pr.trt. Value 'at' must be found among 'x'.
#' @param pr.con (unnormalized) probability mass function of the control group.
#' @param pr.trt (unnormalized) probability mass function of the treatment group.
#' @param verbose (default TRUE), if TRUE, prints info on the console.
#' @return BQTE at value 'at'
#' @examples
#' exact_bqte_discrete(3, x = 1:4, pr.con = c(0.1,0.2,0.3,0.4), pr.trt = c(0.3,0.5,0.15,0.05))
#' @export
exact_bqte_discrete <- function(at, x, pr.con, pr.trt, verbose = TRUE){

  ii = which(x == at)
  if(length(ii) != 1) stop("Didn't find 'at' from 'x'.")
  if(cumsum(pr.con[ii]) < 1e-10) stop("Value of 'at' must be in support of 'pr.con'.")
  if(length(x) != length(pr.con) | length(x) != length(pr.trt)) {
    stop("'x', 'pr,con' and 'pr.trt' must have same length.")}
  if(any(pr.con < 0)) stop("pr.con must have non-negative values")
  if(any(pr.trt < 0)) stop("pr.con must have non-negative values")
  #normalize probabilities
  pr.con = pr.con/sum(pr.con)
  pr.trt = pr.trt /sum(pr.trt)
  #order values
  i.ordered = order(x)
  x = x[i.ordered]
  pr.con = pr.con[i.ordered]
  pr.trt = pr.trt[i.ordered]

  if(verbose){
    cat("Running with distributions:\n")
    print(rbind(x = x, pr.con = pr.con, pr.trt = pr.trt))
    cat("\n")}

  low = 0 #lower probability mass where value 'at' starts in controls
  if(ii > 1) low = sum(pr.con[1:(ii-1)])
  up = low + pr.con[ii] #upper probability mass where value 'at' ends in controls
  # Probability region(low,up) in the treated may include many values
  # Find them all and weight them according to the probabilities
  trt.low.i = min(which(cumsum(pr.trt) >= low-1e-10)) #first trt index in the interval
  trt.up.i = min(which(cumsum(pr.trt) >= up-1e-10)) #last trt index in the interval
  if(verbose){
    cat(paste0("At value",at,".\n"))
    cat(paste0("Quantile levels from controls: (",low,",",up,"].\n"))
    cat(paste0("Treatment values in that quantile range: ",x[trt.low.i],",...,",x[trt.up.i],".\n"))
  }
  if(trt.low.i == trt.up.i) {
    trt.mean = x[trt.low.i] # interval includes only one value in trt
  }else{
    #add contribution from first value
    trt.mean = (sum(pr.trt[1:trt.low.i]) - low)*x[trt.low.i]
    trt.mean = trt.mean + (pr.trt[trt.up.i]-(sum(pr.trt[1:trt.up.i]) - up))*x[trt.up.i]
    if(trt.up.i-trt.low.i > 1){
      ind = (trt.low.i+1):(trt.up.i-1)
      trt.mean = trt.mean + sum(pr.trt[ind]*x[ind])
    }
    trt.mean = trt.mean/(up-low)
  }
  return(trt.mean - at)
}

#MIT License

#Copyright (c) 2022-2025 Matti Pirinen, Harri Hemila

#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:

#  The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.
