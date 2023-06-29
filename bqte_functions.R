# These functions compute
# Quantile Treatment Effect (QTE) and
# Back-transformed Quantile Treatment Effect (BQTE). 
# BQTE is a transformation of QTE as a function of the outcome value in control group. 
# Uses bootstrap for estimating uncertainty.

# Matti Pirinen & Harri Hemila, version June 26, 2023.

# The source code in this file is published under 
# the MIT license listed at the end of this file.

# Contains 8 functions:
#
# bqte(Treatment, Control, at = NULL, 
#      bqte.conf = 0.95, B = 2000, 
#      K = length(Control), bagging = TRUE,
#      tails = FALSE, verbose = TRUE)
#
# plot_bqte(x, plot.ci = TRUE, plot.rbqte = FALSE,
#           pch = 19, col = "black",
#           col.ci = "gray",
#           xlim = NULL, ylim = NULL, 
#           xlab = NULL, ylab = NULL,
#           xaxs = "r", yaxs = "r",
#           xaxp = NULL, yaxp = NULL,
#           cex = 1, cex.axis = 1, cex.lab = 1, 
#           main ="", lwd = 1)
#
# plot_tbqte(x, utbqte = TRUE, plot.ci = TRUE, 
#            plot.relative = FALSE,pch = NULL, col = "black", 
#            col.ci = "gray", xlim = NULL, ylim = NULL, 
#            xlab = NULL, ylab = NULL, xaxs = "r", yaxs = "r",
#            xaxp = NULL, yaxp = NULL,
#            cex = 1, cex.axis = 1, cex.lab = 1, 
#            main ="", lwd = 1)
#
# qte(Treatment, Control, at = NULL, 
#     qte.conf = 0.95, B = 2000, 
#     bagging = TRUE, verbose = TRUE)
#  
# plot_qte(x, plot.ci = TRUE, plot.rqte = FALSE,
#          pch = 19, col = "black",
#          col.ci = "gray",
#          xlim = NULL, ylim = NULL, 
#          xlab = NULL, ylab = NULL,
#          xaxs = "r", yaxs = "r",
#          xaxp = NULL, yaxp = NULL)
#  
# quantiles.from.density(p, pr)
#
# exact.bqte.discrete(at, x, pr.con, pr.trt)
#
# bqte.doksum(at, Treatment, Control)
#             
#
# For examples of usage, see file 'bqte_examples.R'.  

bqte <- function(Treatment, Control, at = NULL, 
                 bqte.conf = 0.95, B = 2000, 
                 K = length(Control), bagging = TRUE,
                 tails = FALSE, verbose = TRUE){
  
  #Goal: 
  # We have observed survival/recovery times, or other numerical outcome values, 
  # in the treatment group and in the control group. 
  # Sample sizes can differ between the groups.
  # Quantile treatment effect (QTE) is the difference between the two groups as function of quantiles.
  # We want to transform QTE to the original scale of the outcome value,
  # at a given grid of values ('at') in the control group.
  # These transformed values are called 'backtransformed quantile treatment effects' (BQTE).
  
  #Approach:
  # Fix a grid of outcome values in the control group at which
  #  BQTEs are computed (parameter 'at').
  #  Recommendation: 
  #   Do not try to estimate BQTEs at values that are outside
  #   empirical quantile levels (5/N, 1 - (5/N)), 
  #   where N is the sample size of controls because of low accuracy at the tails. 
  #
  # Estimate the quantile treatment effects as difference between the quantiles 
  #  of Treatment group and Control group at 'K' quantile levels 
  #  1/(K+1),...,K/(K+1). 
  #  Make a piece-wise linear approximation of BQTE
  #  as a function of outcome value in controls, 
  #  based on the observed K quantiles. 
  #  Use that approximation
  #  to estimate the BQTE at the chosen grid of outcome values in controls.
  #
  # Use bootstrapping to estimate the uncertainty of BQTE at each grid point,
  #  by resampling with replacement a set of Treatment and Control group values that
  #  have the same sample size as the original Treatment and Control groups, respectively,
  #  and applying the procedure explained above to each bootstrap sample.
  #  Note that bootstrapping accounts for uncertainty of quantiles of both groups.
  # 
  # The final estimate of BQTE at a grid point is the mean over bootstrap samples ("bagging estimate")
  #  and the confidence interval for BQTE is estimated from the quantiles of the bootstrap sample.
  # NOTE: In empirical evaluation it was observed that bagging estimate tend to have smaller
  #       mean square error in discrete distributions than the direct estimate from the data and
  #       therefore the bagging estimate is used by default. This can be switched by parameter "bagging".
  #
  # Can also compute upper tail BQTE (UTBQTE) and lower tail BQTE (LTBQTE) at each point in 'at'.
  # UTBQTE, at outcome value 'x' is the difference in the mean outcome values 
  #  of the highest 1-F(x) of the two populations, where F is the CDF of the Control population.
  # LTBQTE, at outcome value 'x' is the difference in the mean outcome values 
  #  of the lowest F(x) of the two populations, where F is the CDF of Control population.
  # UTBQTE(x) is an upper bound for the average treatment effect of that part of the Control population who
  # have outcome value >= x.
  # LTBQTE(x) is a lower bound for the average treatment effect of that part of the Control population who
  # have outcome value <= x.
  # These bounds are general in the sense that they do not require an assumption about how the treatment operates,
  # for example, they do not assume anything about whether the treatment preserves order.
  # Note that the outcome value 'x' is included in the definition of both UTBQTE(x) and in LTBQTE(x).
  # The confidence intervals of UTBQTE and LTBQTE are estimated using bootstrap
  #  and the point estimates are computed either using bagging (if bagging = TRUE) or not (if bagging = FALSE). 

  
  #INPUT:  
  # 'Treatment', vector of outcome values (e.g. survival times) in Treatment group. 
  # 'Control', vector of outcome values in Control group.
  # 'at', vector of outcome values in controls at which BQTE is estimated.
  #       By default, between 10 and 20 points between the empirical quantile levels 
  #       of 5/N and 1 - (5/N) in controls where N = length(Control).
  # 'bqte.conf' (default 0.95), confidence level for BQTE
  # 'B' (default 2000), number of bootstrap samples
  # 'K' (default 'length(Control)), number of quantiles where BQTE function is estimated.
  # 'bagging' (default TRUE), if TRUE, reports bagging estimate (mean over bootstrap samples)
  #                           otherwise reports the direct estimate from the data.
  # 'tails' (default FALSE), if TRUE, returns also estimates and confidence intervals for 
  #                         UTBQTE and LTBQTE.
  # 'verbose' (default TRUE) if TRUE, prints parameters and recommendations on console.
  
  #OUTPUT:
  #  data frame with 7 or 13 columns
  #1  'at', outcome value in control group
  #2  'bqte', backtransformed quantile treatment effect
  #3  'bqte.low', lower CI end point for 'bqte' 
  #4  'bqte.up', upper CI end point for 'bqte'
  #5  'rbqte', relative backtransformed quantile treatment effect with respect to control outcome value
  #6  'rbqte.low', lower CI end point for 'rbqte' 
  #7  'rbqte.up', upper CI end point for 'rbqte'
  #   columns 8-19 will be returned only if parameter 'tails = TRUE'
  #8  'utbqte, upper tail average treatment effect estimate
  #9  'utbqte.low, lower CI end point for 'utbqte'
  #10 'utbqte.up, upper CI end point for 'utbqte'
  #11 'ltbqte, lower tail average treatment effect estimate
  #12 'ltbqte.low, lower CI end point for 'ltbqte'
  #13 'ltbqte.up, upper CI end point for 'ltbqte'
  #14 'rutbqte, relative upper tail average treatment effect estimate with respect to control value
  #15 'rutbqte.low, lower CI end point for 'rutbqte'
  #16 'rutbqte.up, upper CI end point for 'rutbqte'
  #17 'rltbqte, relative lower tail average treatment effect estimate with respect to control value
  #18 'rltbqte.low, lower CI end point for 'rltbqte'
  #19 'rltbqte.up, upper CI end point for 'rltbqte'
  
  if( K < 2 ) stop("K should be at least 2.")
  if( bqte.conf <= 0 | bqte.conf >= 1 ) stop("Value \'bqte.conf\' not valid.")
  if( B < 1 ) stop("B is not positive integer.")
  C.n = length(Control)
  T.n = length(Treatment)
  
  #recommended to keep value in 'at' within the range 'x'
  x = as.numeric(quantile(Control, prob = c(5/C.n, 1 - 5/C.n)))
  if( is.null(at) ){ # Set 'at' to its default value 
    # that has between 10 and 20 values from interval between empirical quantiles
    # of Control vector of levels 5/C.n and 1 - 5/C.n.
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
    cat(paste0("If you choose to change 'at' values, a recommended interval is [",
               signif(x[1],3),", ",signif(x[2],3),"].\n"))
  }
  probs = (1:K)/(K+1) # quantile levels at which the BQTE-function is estimated
  bqte.b = matrix(NA, nrow = B, ncol = length(at)) # bootstrapped  BQTEs at grid 'at'
  utbqte.b = matrix(NA, nrow = B, ncol = length(at)) # bootstrapped UTbqtes at grid 'at'
  rutbqte.b = matrix(NA, nrow = B, ncol = length(at)) # bootstrapped relative UTbqtes at grid 'at'
  ltbqte.b = matrix(NA, nrow = B, ncol = length(at)) # bootstrapped LTbqtes at grid 'at'
  rltbqte.b = matrix(NA, nrow = B, ncol = length(at)) # bootstrapped relative LTbqtes at grid 'at'
  
  for (ii in 1:B){
    
    #1. Make bootstrap sample
    C.b = sample(Control, size = C.n, replace = T)
    T.b = sample(Treatment, size = T.n, replace = T)
    
    #2. Pick the fixed set of quantiles from the bootstrap sample
    C.q = as.numeric(quantile(C.b, prob = probs))
    T.q = as.numeric(quantile(T.b, prob = probs))
    
    #3. Estimate BQTE at grid 'at' using linear approximation
    bqte.b[ii,] = as.numeric(approx(C.q, T.q - C.q, xout = at, rule = 2, ties = mean)$y)
    
    #4. Estimate utbqte and ltbqte at grid 'at' using linear approximation (if 'tails = TRUE')
    if(tails){
      ut.C.mean.q = sapply(C.q, function(x){mean(C.b[C.b >= x])})
      ut.T.mean.q = sapply(T.q, function(x){mean(T.b[T.b >= x])})

      utbqte.q = ut.T.mean.q  - ut.C.mean.q
      utbqte.b[ii,] = approx(C.q, utbqte.q, xout = at, rule = 2, ties = mean)$y
      rutbqte.q = (ut.T.mean.q  - ut.C.mean.q)/ut.C.mean.q
      rutbqte.b[ii,] = approx(ut.C.mean.q, rutbqte.q, xout = at, rule = 2, ties = mean)$y
      
      lt.C.mean.q = sapply(C.q, function(x){mean(C.b[C.b <= x])})
      lt.T.mean.q = sapply(T.q, function(x){mean(T.b[T.b <= x])})

      ltbqte.q = lt.T.mean.q  - lt.C.mean.q
      ltbqte.b[ii,] = approx(C.q, ltbqte.q, xout = at, rule = 2, ties = mean)$y
      rltbqte.q = (lt.T.mean.q  - lt.C.mean.q)/lt.C.mean.q
      rltbqte.b[ii,] = approx(lt.C.mean.q, rltbqte.q, xout = at, rule = 2, ties = mean)$y
    }
  }# ends bootstrap

  bqte.low = as.numeric(apply(bqte.b, 2, function(X){quantile(X, (1 - bqte.conf)/2)}))
  bqte.up = as.numeric(apply(bqte.b, 2, function(X){quantile(X, (1 + bqte.conf)/2)}))
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
    C.q = as.numeric(quantile(Control, prob = probs))
    T.q = as.numeric(quantile(Treatment, prob = probs))
    bqte = as.numeric(approx(C.q, T.q - C.q, xout = at, rule = 2, ties = mean)$y)
    if(tails){
      ut.T.mean.q = sapply(C.q, function(x){mean(Treatment[Treatment >= x])})
      ut.C.mean.q = sapply(T.q, function(x){mean(Control[Control >= x])})
      utbqte.q = ut.T.mean.q  - ut.C.mean.q
      utbqte = as.numeric(approx(C.q, utbqte.q, xout = at, rule = 2, ties = mean)$y)
      rutbqte.q = ut.T.mean.q/ut.C.mean.q - 1
      rutbqte = as.numeric(approx(C.q, utbqte.q, xout = at, rule = 2, ties = mean)$y)
      
      lt.T.mean.q = sapply(C.q, function(x){mean(Treatment[Treatment <= x])})
      lt.C.mean.q = sapply(T.q, function(x){mean(Control[Control <= x])})
      ltbqte.q = lt.T.mean.q  - lt.C.mean.q
      ltbqte = as.numeric(approx(C.q, ltbqte.q, xout = at, rule = 2, ties = mean)$y)
      rltbqte.q = lt.T.mean.q/lt.C.mean.q - 1
      rltbqte = as.numeric(approx(C.q, ltbqte.q, xout = at, rule = 2, ties = mean)$y)
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

  if(tails){
    res = data.frame(res,
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
  
  return(res)
}


plot_bqte <- function(x, plot.ci = TRUE, plot.rbqte = FALSE,
                      pch = 19, col = "black", 
                      col.ci = "gray",
                      xlim = NULL, ylim = NULL, 
                      xlab = NULL, ylab = NULL,
                      xaxs = "r", yaxs = "r",
                      xaxp = NULL, yaxp = NULL,
                      cex = 1, cex.axis = 1, cex.lab = 1, 
                      main ="", lwd = 1){
  
  #Print backtransformed quantile treatment effect (BQTE) on direct or relative scale
  #INPUT:
  # 'x' data.frame returned by function bqte( )
  # plot.ci, if TRUE, plots confidence intervals around estimates using color col.ci
  # plot.rbqte, if TRUE, plots relative BQTEs, otherwise plots actual BQTEs
  # pch, col, col.ci, xlim, ylim, xlab, ylab, xaxs, yaxs, xaxp, yaxp: 
  #      standard plotting parameters with sensible defaults.
  # if pch = NULL then the estimates are shown by a line and confidence intervals by a band.
  
  #Use percentages for relative effects
  if(plot.rbqte) {
    x[,c("rbqte","rbqte.low","rbqte.up")] = 100*x[,c("rbqte","rbqte.low","rbqte.up")]
    y = x$rbqte; y.low = x$rbqte.low; y.up = x$rbqte.up}
  else{y = x$bqte; y.low = x$bqte.low; y.up = x$bqte.up}

  if(is.null(xlim)) {
    xlim = range(x$at) + c(-1,1) * 0.05 * diff(range(x$at))}
  if(is.null(ylim)) {
      ylim = range(c(y.low,y.up)) + c(-1,1) * 0.05 * diff(range(c(y.low,y.up)))}
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
  if(plot.ci){
    if(!is.null(pch)) {
      arrows(x$at, y.low, 
             x$at, y.up, 
             col = col.ci, code = 3, angle = 90, length = 0.0, lwd = lwd)}
    else{
      polygon(c(x$at,rev(x$at)), 
              c(y.low, rev(y.up)), 
              col = col.ci, border = NA)
    }
  }
  
  #add points
  if(!is.null(pch)) {
    points(x$at, y, pch = pch, col = col, cex = cex)}
  else {
    lines(x$at, y, lty = 1, col = col, lwd = lwd)}
}


plot_tbqte <- function(x, utbqte = TRUE, plot.ci = TRUE, plot.relative = FALSE,
                      pch = NULL, col = "black", 
                      col.ci = "gray",
                      xlim = NULL, ylim = NULL, 
                      xlab = NULL, ylab = NULL,
                      xaxs = "r", yaxs = "r",
                      xaxp = NULL, yaxp = NULL,
                      cex = 1, cex.axis = 1, cex.lab = 1, 
                      main ="", lwd = 1){
  
  # Plot tail BQTEs on direct or relative scale.
  # For definitions of these measures, see "bqte_supplement.pdf"
  # Either upper (UTBQTE) or lower (LTBQTE) tail can be considered.
  #INPUT:
  # 'x' data.frame returned by function bqte( ,tails = TRUE)
  # 'utbqte', if TRUE, prints UTBQTE else prints LTBQTE
  # plot.ci, if TRUE, plots confidence intervals around estimates using color col.ci.
  # plot.relative, if TRUE, plots relative TBQTEs, otherwise plots direct TBQTEs
  # pch, col, col.ci, xlim, ylim, xlab, ylab, xaxs, yaxs, xaxp, yaxp: 
  #      standard plotting parameters with sensible defaults.
  # If pch = NULL then the estimates are shown by a line and confidence intervals by a band.
  
  
  if(utbqte) {
    if(plot.relative){
      X = x[,c("rutbqte","rutbqte.low","rutbqte.up")]}
    else{
      X = x[,c("utbqte","utbqte.low","utbqte.up")]}
  }else {
    if(plot.relative){
      X = x[,c("rltbqte","rltbqte.low","rltbqte.up")]}
    else{
      X = x[,c("ltbqte","ltbqte.low","ltbqte.up")]}
  }
  
  #Use percentages for relative effects
  if(plot.relative) X = 100*X
  
  if(is.null(xlim)) {
    xlim = range(x$at) + c(-1,1) * 0.05 * diff(range(x$at))}
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
      arrows(x$at, X[,2],
             x$at, X[,3], 
             col = col.ci, code = 3, angle = 90, length = 0.0, lwd = lwd)}
    else{
      polygon(c(x$at,rev(x$at)), 
              c(X[,2], rev(X[,3])), 
              col = col.ci, border = NA)
    }
  }
  
  #add points
  if(!is.null(pch)) {
    points(x$at, X[,1], pch = pch, col = col, cex = cex)}
  else{
    lines(x$at, X[,1], col = col, lty = 1, lwd = lwd)}
}


qte <- function(Treatment, Control, at = NULL, 
                qte.conf = 0.95, B = 2000, 
                bagging = TRUE, verbose = TRUE){
  
  #Goal: 
  # We have observed survival/recovery times, or other continuous outcome values, 
  # in the treatment group and in the control group. 
  # Sample sizes can differ between the groups.
  # Quantile treatment effect (QTE) is the difference between the two groups as functions of quantiles.

  #Approach:
  # Fix a grid of quantiles where QTEs are computed (parameter 'at').
  #  Recommendation: 
  #   Do not try to estimate QTEs at values that are outside
  #   quantile levels (2.5/N, 1 - (2.5/N)), 
  #   where N is the minimum of sample sizes of treated and controls because of low accuracy at the tails. 
  #
  # Estimate the quantile treatment effects as difference between the quantiles 
  #  of Treatment group and Control group 
  #
  # Use bootstrapping to estimate the uncertainty of QTE,
  #  by resampling with replacement a set of Treatment and Control group values that
  #  have the same sample size as the original Treatment and Control groups, respectively,
  #  and applying the procedure explained above to each bootstrap sample.
  #  Note that bootstrapping accounts for uncertainty of quantiles of both groups.
  # 
  # The final estimate of QTE is the mean over bootstrap samples ("bagging estimate")
  #  and the confidence interval for QTE is estimated from the quantiles of the bootstrap sample.
  # NOTE: In empirical evaluation it was observed that bagging estimate tend to have smaller
  #       mean square error in discrete distributions than the direct estimate from the data and
  #       therefore the bagging estimate is used by default. This can be switched by parameter "bagging".
  
  #INPUT:  
  # 'Treatment', vector of outcome values (e.g. survival times) in Treatment group. 
  # 'Control', vector of outcome values in Control group.
  # 'at', vector of quantiles at which QTE is estimated.
  #       By default, between 10 and 20 points between the quantile levels 
  #       of 2.5/N and 1 - (2.5/N) where N = min(c(length(Control),length(Treatment))).
  # 'qte.conf' (default 0.95), confidence level for QTE
  # 'B' (default 2000), number of bootstrap samples
  # 'bagging' (default TRUE), if TRUE, reports bagging estimate (mean over bootstrap samples)
  #                           otherwise reports the direct estimate from the data.
  # 'verbose' (default TRUE) if TRUE, prints parameters and recommendations on console.
  
  #OUTPUT:
  #  data frame with 7 columns
  #1 'at', quantile level
  #2 'qte', quantile treatment effect
  #3 'qte.low', lower CI end point for 'qte' 
  #4 'qte.up', upper CI end point for 'qte'
  #5 'rqte', relative quantile treatment effect with respect to control outcome value
  #6 'rqte.low', lower CI end point for 'rqte' 
  #7 'rqte.up', upper CI end point for 'rqte'
  
  if( qte.conf <= 0 | qte.conf >= 1 ) stop("Value \'qte.conf\' not valid.")
  if( B < 1 ) stop("B is not positive integer.")
  C.n = length(Control)
  T.n = length(Treatment)
  
  #recommended to keep value in 'at' within the range 'x'
  min.n = min(c(C.n,T.n)) 
  if(min.n < 5) stop("Minimum allowed sample size is 5 in both groups.")
  if( is.null(at) ){ # Set 'at' to its default value 
    # that has at most 19 values from interval between quantile levels 2.5/min.n and 1 - 2.5/min.n.
    K = max(which(1/(2:20) >= 2.5/min.n))
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
    x = c(2.5/min.n, 1 - 2.5/min.n)
    cat(paste0("If you choose to change 'at' values, a recommended interval is [",
               signif(x[1],3),", ",signif(x[2],3),"].\n"))
  }
  C.b = matrix(NA, nrow = B, ncol = length(at)) # bootstrapped control quantiles
  T.b = matrix(NA, nrow = B, ncol = length(at)) # bootstrapped treatment quantiles
  QTE.b = matrix(NA, nrow = B, ncol = length(at)) # bootstrapped approximations of QTEs at grid 'at'
  for (ii in 1:B){
    C.b[ii,] = as.numeric(quantile(sample(Control, size = C.n, replace = T), prob = at))
    T.b[ii,] = as.numeric(quantile(sample(Treatment, size = T.n, replace = T), prob = at))
    QTE.b[ii,] = as.numeric(T.b[ii,] - C.b[ii,])
  }
  
  qte.low = as.numeric(apply(QTE.b, 2, function(X){quantile(X, (1 - qte.conf)/2)}))
  qte.up = as.numeric(apply(QTE.b, 2, function(X){quantile(X, (1 + qte.conf)/2)}))
  
  if(bagging) { # Mean as the final estimate 
    C.q = apply(C.b, 2, mean)
    qte = apply(QTE.b, 2, mean)
  } 
  else{ #Compute QTE estimate directly from the data (not bagging)
    C.q = as.numeric(quantile(Control, prob = at))
    T.q = as.numeric(quantile(Treatment, prob = at))
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


plot_qte <- function(x, plot.ci = TRUE, plot.rqte = FALSE,
                     pch = 19, col = "black",
                     col.ci = "gray",
                     xlim = NULL, ylim = NULL, 
                     xlab = NULL, ylab = NULL,
                     xaxs = "r", yaxs = "r",
                     xaxp = NULL, yaxp = NULL){
  
  #Print quantile treatment effect (QTE) on direct or relative scale
  #INPUT:
  # 'x' data.frame returned by function qte( )
  # plot.ci, if TRUE, plots confidence intervals around estimates
  # plot.rqte, if TRUE, plots relative QTEs, otherwise plots actual QTEs
  # pch, col, col.ci, xlim, ylim, xlab, ylab, xaxs, yaxs, xaxp, yaxp: 
  #      standard plotting parameters with sensible defaults.
  
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
  abline(h = 0, lty = 2)
  
  #add intervals
  if(plot.ci){
    if(!plot.rqte){y.low = x$qte.low; y.up = x$qte.up}
    else{y.low = x$rqte.low; y.up = x$rqte.up}
    arrows(x$at, y.low, 
           x$at, y.up, 
           col = col.ci, code = 3, angle = 90, length = 0.0)
  }
  
  #add points
  if(!plot.rqte) {
    points(x$at, x$qte, pch = pch, col = col)}
  else {
    points(x$at, x$rqte, pch = pch, col = col)}
}


bqte.doksum <- function(at, Treatment, Control){
  # Estimator of BQTE given by Doksum (1974) Annals of Statistics Vol.2 No.2 267-277.
  y = sapply(at, function(x){floor(mean(Control <= x) * length(Treatment)) + 1})
  y[y > length(Treatment)] = length(Treatment)
  sort(Treatment)[y] - at
}
            

quantiles.from.density <- function(p, pr){ 
  # Returns the index of a discrete probability mass function 'pr' corresponding to 
  # the quantile of probability level 'p'
  min(which(cumsum(pr) >= p)) }


exact.bqte.discrete <- function(at, x, pr.con, pr.trt){
  # Computes true BQTE for discrete distribution
  # INPUT
  # at, value at which bqte is computed
  # x, values corresponding to elements of pr.con and pr.trt
  #    'at' must be found among 'x'
  # pr.con, probability mass function of the control group
  # pr.trt, probability mass function of the treatment group
  # OUTPUT
  # estimate of BQTE at value 'at'
  
  ii = which(x == at)
  if(length(ii) != 1) stop("Didn't find 'at' from 'x'.")
  if(length(x) != length(pr.con) | length(x) != length(pr.trt)) {
    stop("'x', 'pr,con' and 'pr.trt' must have same length.")}
  
  low = 0 #lower probability mass where value 'at' starts in controls
  if(ii > 1) low = sum(pr.con[1:(ii-1)])
  up = low + pr.con[ii] #upper probability mass where value 'at' ends in controls
  # Probability region(low,up) in the treated may include many values
  # Finds the expectation over the region empirically using 'n' unit masses
  n = 10000
  y = rep(x, round(pr.trt*n)) # Each unit mass will be attached to a value
  if(length(y) > n) y = y[-sample(1:n, size = length(y) - n)] #remove random elements to get 'n' elements
  if(length(y) < n) y = sort(c(y, sample(y, size = n - length(y), replace = TRUE)))# add random elements to get 'n' elements
  return(mean(y[round(low*n):round(up*n)]) - at) #compute BQTE
}

#MIT License

#Copyright (c) 2022,2023 Matti Pirinen, Harri Hemila

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
