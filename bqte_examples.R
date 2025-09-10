# Example analyses using bqte(), qte() and plotting the results using
# plot_bqte(), plot_qte(), plot_tbqte() plot_tbqte_both() functions.
# Harri Hemila, Matti Pirinen
# Aug 19, 2025.

# Contents:
# 1. Read in the functions
# 2. Examples of the functions
# 3. Figures 1 - 4 of the main paper
# 4. Compare BQTE and QTE

#
##
### 1. Read in the functions
##
#

# To install bqte from GitHub,
#  you need to install 'devtools' package (in case you don't have it).
install.packages("devtools")
library(devtools)
install_github("mjpirinen/bqte")
library(bqte)
#  If you can't get the above working, all functions are in a single R file,
#  which you can also read in directly to R using command below
# source("https://raw.githubusercontent.com/mjpirinen/bqte/main/R/bqte.R")

# Data files can be found from GitHub
# https://github.com/mjpirinen/bqte folder "data"
# Below, direct web links are provided.



#
##
### 2. Examples of the functions
##
#



# Analysis of Mossad data on effect of zinc gluconate lozenges
# on duration of common cold.
# Data taken from pages 2-4 of the additional file 2 of
# https://doi.org/10.1186/s12874-017-0356-y

filename.1 = "https://raw.githubusercontent.com/mjpirinen/bqte/main/data/Mossad.csv"
Mossad = read.csv(filename.1)
str(Mossad) # Check data looks OK.
table(Mossad$Duration, Mossad$Zinc)

Treatment = subset(Mossad, Zinc == 1)$Duration
table(Treatment)
Control = subset(Mossad, Zinc == 0)$Duration
table(Control)

set.seed(1) #Set random seed for reproducible results
res.bqte = bqte(Treatment = Treatment, Control = Control, tails = TRUE) #Estimate BQTE

#Check output of bqte() function
res.bqte

ate = mean(Treatment) - mean(Control) #average treatment effect
rate = 100*ate/mean(Control) #relative ATE
c(ate = ate, rel.ate = rate)

par(mfrow = c(1,2))

#Plot results on outcome scale ('plot.rte = FALSE')
plot_bqte(res.bqte, plot.rbqte = FALSE,
          xlim = c(0,20), xaxs = "i", col = "red",
          col.ci = "black",
          ylab = "BQTE (Days)",
          xlab = "Duration in the control group (Days)")
abline(ate, 0, col = "blue", lt = 3,  lw = 2) #Add ATE line
text(18, 0, "A", pos = 1, cex = 3) #Label this one as panel "A"

#Plot results on relative scale ('plot.rte = T')
plot_bqte(res.bqte, plot.rbqte = TRUE,
          xlim = c(0,20), xaxs = "i", col = "red",
          col.ci = "black",
          xlab="Duration in the control group (Days)")
abline(rate, 0, col = "blue", lt = 3, lw = 2) #Add relative ATE  line
text(18, 0, "B", pos = 1, cex = 3) #Label this one as panel "B"

#
# Note: Same figures but with confidence intervals connected can be drawn by setting 'pch = NULL'
#       These are still pointwise intervals, not joint confidence regions.

plot_bqte(res.bqte, plot.rbqte = FALSE, pch = NULL,
          xlim = c(0,20), xaxs = "i", col = "red",
          col.ci = "gray",
          ylab = "BQTE (Days)",
          xlab = "Duration in the control group (Days)")
abline(ate, 0, col = "blue", lt = 3,  lw = 2) #Add ATE line
text(18, 0, "A", pos = 1, cex = 3) #Label this one as panel "A"

#Plot results on relative scale ('plot.rte = TRUE')
plot_bqte(res.bqte, plot.rbqte = TRUE, pch = NULL,
          xlim = c(0,20), xaxs = "i", col = "red",
          col.ci = "gray",
          xlab="Duration in the control group (Days)")
abline(rate, 0, col = "blue", lt = 3, lw = 2) #Add relative ATE  line
text(18, 0, "B", pos = 1, cex = 3) #Label this one as panel "B"


# Since bqte() was run with parameter 'tails = TRUE',
# we can also plot the upper tail BQTE (UTBQTE) estimate with 95%CI
# which is an upper bound for the upper tail average treatment effect (UTATE),
# and the lower tail BQTE (LTBQTE) which is a lower bound for
# the lower tail average treatment effect (LTATE).
# Let's demonstrate these on both direct and relative scales.

par(mfrow = c(2,2))
#PLOT UTBQTEs
plot_tbqte(res.bqte, utbqte = TRUE, plot.relative = F, pch = 25,
          xlim = c(0,20), xaxs = "i", col = "blue",
          col.ci = "blue")
abline(ate, 0, col = "black", lt = 3,  lw = 2) #Add known mean effect line
plot_tbqte(res.bqte, utbqte = TRUE, plot.relative = T, pch = 25,
           xlim = c(0,20), xaxs = "i", col = "blue",
           col.ci = "blue")
abline(rate, 0, col = "black", lt = 3, lw = 2) #Add known mean effect line

#PLOT LTBQTEs
plot_tbqte(res.bqte, utbqte = FALSE, plot.relative = FALSE, pch = 24,
           xlim = c(0,20), xaxs = "i", col = "purple",
           col.ci = "purple")
abline(ate, 0, col = "black", lt = 3,  lw = 2) #Add known mean effect line
plot_tbqte(res.bqte, utbqte = FALSE, plot.relative = TRUE, pch = 24,
           xlim = c(0,20), xaxs = "i", col = "purple",
           col.ci = "purple")
abline(rate, 0, col = "black", lt = 3, lw = 2) #Add known mean effect line

#PLOT LTBQTE and UTBQTE in the same plot
# On the original scale:
col.utbqte = "blue"
col.ltbqte = "purple"
par(mfrow = c(1,2))
plot_tbqte_both(
  res.bqte,
  plot.ci = TRUE, plot.relative = FALSE,
  pch.utbqte = 25, pch.ltbqte = 24,
  col.utbqte = col.utbqte,
  col.ltbqte = col.ltbqte,
  col.utbqte.ci = col.utbqte,
  col.ltbqte.ci = col.ltbqte,
  jitter = c(0.2, -0.2),
  xlim = NULL, ylim = NULL,
  xlab = NULL, ylab = NULL,
  xaxs = "r", yaxs = "r",
  xaxp = NULL, yaxp = NULL,
  cex = 1, cex.axis = 1, cex.lab = 1,
  main ="", lwd = 1.2)
abline(ate, 0, col = "black", lt = 3,  lw = 2) #Add known mean effect line

# On the relative scale:
plot_tbqte_both(
  res.bqte,
  plot.ci = TRUE, plot.relative = TRUE,
  pch.utbqte = 25, pch.ltbqte = 24,
  col.utbqte = col.utbqte,
  col.ltbqte = col.ltbqte,
  col.utbqte.ci = col.utbqte,
  col.ltbqte.ci = col.ltbqte,
  jitter = c(0.2, -0.2),
  xlim = NULL, ylim = NULL,
  xlab = NULL, ylab = NULL,
  xaxs = "r", yaxs = "r",
  xaxp = NULL, yaxp = NULL,
  cex = 1, cex.axis = 1, cex.lab = 1,
  main ="", lwd = 1.2)
abline(rate, 0, col = "black", lt = 3,  lw = 2) #Add known mean effect line


# We can also estimate quantile treatment effects (QTEs) and plot them as follows

res.qte = qte(Treatment = Treatment, Control = Control) #Estimate QTE

# Check output of qte() function
res.qte

par(mfrow = c(1,2))
plot_qte(res.qte, plot.rqte = FALSE, col = "limegreen", col.ci = "limegreen")
abline(ate, 0, col = "black", lt = 3,  lw = 2) #Add known mean effect line

plot_qte(res.qte, plot.rqte = TRUE, col = "limegreen", col.ci = "limegreen")
abline(rate, 0, col = "black", lt = 3,  lw = 2) #Add known mean effect line



####

# Analysis of data of the three randomized trials
# on zinc acetate lozenges on duration of common cold.
# The data are from page 8 of supplementary file 2 of
# https://doi.org/10.1093/ofid/ofx059

filename.2 = "https://raw.githubusercontent.com/mjpirinen/bqte/main/data/ZnAcet.csv"
ZnAcet = read.csv(filename.2)
str(ZnAcet)
table(ZnAcet$Duration, ZnAcet$Zinc)

Treatment = subset(ZnAcet, Zinc == 1)$Duration
Control = subset(ZnAcet, Zinc == 0)$Duration

ate = mean(Treatment) - mean(Control) #average treatment effect
rate = 100*ate/mean(Control)
c(ate = ate, rel.ate = rate)

set.seed(2)
res.bqte = bqte(Treatment = Treatment, Control = Control, tails = TRUE)

res.bqte

par(mfrow = c(1,2))
plot_bqte(res.bqte, plot.rbqte = FALSE,
          xaxs = "i", xlim = c(0,15), col = "red",
          col.ci = "black",
          ylab = "BQTE (Days)",
          xlab="Duration in the control group (Days)")
abline(ate, 0, col = "blue", lt = 3,  lw = 2)
text(14, -0.6, "A", pos = 1, cex = 3)

plot_bqte(res.bqte, plot.rbqte = TRUE,
          xaxs = "i", xlim = c(0,15), col = "red",
          col.ci = "black",
          xlab="Duration in the control group (Days)")
abline(rate, 0, col = "blue", lt = 3, lw = 2)
text(14, -13, "B", pos = 1, cex = 3)




####

# Analysis of data of the two nasal iota-carrageenan trials
# on duration of common cold.
# The data are from page 3 of the supplementary file of
# https://doi.org/10.1002/prp2.810
# sensored data on day 20 are imputed with day 20
# (6 in carrageenan, 21 in placebo).


filename.3 = "https://raw.githubusercontent.com/mjpirinen/bqte/main/data/Carrageenan.csv"
Carrageenan = read.csv(filename.3)
str(Carrageenan)
table(Carrageenan$Duration,Carrageenan$Carr)

Treatment = subset(Carrageenan, Carr == 1)$Duration
Control = subset(Carrageenan, Carr == 0)$Duration

ate = mean(Treatment) - mean(Control) #average treatment effect
rate = 100*ate/mean(Control)
c(ate = ate, rel.ate = rate)


set.seed(3)
res.bqte = bqte(Treatment = Treatment, Control = Control, tails = TRUE,
                at = c(3:19)) #at specifies the estimation points
res.bqte

par(mfrow = c(1,2))
plot_bqte(res.bqte, plot.rbqte = FALSE,
          xlim = c(0,20), ylim = c(-9,1),
          yaxp = c(-8, 0, 4),
          col = "red", xaxs = "i",
          col.ci = "black",
          ylab = "BQTE (Days)",
          xlab="Duration in the control group (Days)")
abline(ate, 0, col = "blue", lt = 3, lw = 2)
text(18, 0, "A", pos = 1, cex = 3)

plot_bqte(res.bqte, plot.rbqte = TRUE,
          xlim = c(0,20),  ylim = c(-60,25),
          yaxp = c(-60, 20, 4),
          col = "red", xaxs = "i",
          col.ci = "black",
          xlab="Duration in the control group (Days)")
abline(rate, 0, col = "blue", lt = 3, lw = 2)
text(18, 18, "B", pos = 1, cex = 3)



# Comparison of the Carrageenan data with imputation at 20 days
# (continuous) against imputation at 30 days (dashed line).
# Only difference to above carrageenan data is that
# the censored observations on day 20 are now imputed to be 30 days


filename.4 = "https://raw.githubusercontent.com/mjpirinen/bqte/main/data/Carrageenan30.csv"
Carrageenan30 = read.csv(filename.4)
table(Carrageenan30$Duration, Carrageenan30$Carr)
str(Carrageenan30)

Treatment30 = subset(Carrageenan30, Carr == 1)$Duration
Control30 = subset(Carrageenan30, Carr == 0)$Duration

set.seed(4)
#Use more bootstrap samples to ensure max accuracy (20000 while default = 2000)
res.bqte = bqte(Treatment = Treatment, Control = Control,
                B = 20000, at = c(3:19))
res30.bqte = bqte(Treatment = Treatment30, Control = Control30,
                  B = 20000, at = c(3:19))

res30.bqte

par(mfrow = c(1,1))
plot (res.bqte$res$at, res.bqte$res$bqte, type ="l", xlim = c(0,20),
      yaxp = c(-10, 0, 5), xaxs = "i",
      ylab = "Mean treatment effect (Days)",
      xlab="Duration in the control group (Days)")
lines(res30.bqte$res$at, res30.bqte$res$bqte, lty = 2, col = "red")
legend("bottomleft", col = c("black","red"), lty = c(1,2),
       leg = c("cens. to 20", "cens. to 30"))
# Conclusion:
# No difference observed due to imputation in the range considered here.



#
##
###  3. Figures 1 - 4 of the main paper
##
#


#
### Figure 1
#

plot.pdf = FALSE #set this to TRUE if you want a PDF file, otherwise plotted on screen
if(plot.pdf) pdf("Fig1.pdf", width = 6.6, height = 4.6, pointsize = 9)

filename.1 = "https://raw.githubusercontent.com/mjpirinen/bqte/main/data/Mossad.csv"
Mossad = read.csv(filename.1)
str(Mossad) # Check data looks OK.
table(Mossad$Duration, Mossad$Zinc)

Treatment = subset(Mossad, Zinc == 1)$Duration
table(Treatment)
Control = subset(Mossad, Zinc == 0)$Duration
table(Control)

tt = table(Treatment)
tc = table(Control)
ran = range(c(as.numeric(names(tt)),as.numeric(names(tc))))
pr = matrix(0, ncol = ran[2]-ran[1]+1, nrow = 3)
pr[1,] = ran[1]:ran[2]
pr[2,as.numeric(names(tc))] = as.numeric(tc)/sum(tc)
pr[3,as.numeric(names(tt))] = as.numeric(tt)/sum(tt)
pr
rownames(pr) = c("days","control","treat")
pr

chosen.i = 8 #to highlight for BQTE
exact_bqte_discrete(at = chosen.i, x = pr[1,], pr.con = pr[2,], pr.trt = pr[3,], verbose = TRUE)
#Quantile levels from controls: (0.48,0.58].
#Treatment values in the range: 5,...,6.
#Find the probability mass of 5 and 6 from treatment group that corresponds to these quantile levels
cdf.t = cumsum(pr[3,])
pr.t.chosen = rep(0,length(pr[3,])) #This is probability in treated corresponding to the chosen value in controls
pr.t.chosen[5:6] = c(cdf.t[5]-0.48, 0.58-cdf.t[5]) #MUST BE SET MANUALLY
pr.c.chosen = rep(0, length(pr[2,]))
pr.c.chosen[chosen.i] = pr[2,chosen.i]
#sum(pr.t.chosen*1:19)/sum(pr.t.chosen) #average treatment value to be used in BQTE(8) estimate

cexmain = 1.2
cexaxis = 1.2
cexlab = 1.2
cexnames = 1.2
yran = c(0,max(as.vector(pr[2:3,])))
par(mfcol = c(2,2))
bplot = barplot(rbind(pr[3,]-pr.t.chosen, pr.t.chosen), names.arg = as.integer(pr[1,]),
                main = "Treatment group", ylab = "proportion", cex.names = cexnames,
                ylim = yran, beside = FALSE, col= c("gray","black"), xlab = "Days",
                density = c(1000,30), angle = c(0,45),
                cex.lab = cexlab, cex.axis = cexaxis, cex.main = cexmain)
text(-3, 0.21, "A", pos = 1, cex = 2, xpd = NA)

bplot = barplot(rbind(pr[2,]-pr.c.chosen, pr.c.chosen), names.arg = as.integer(pr[1,]),
                main = "Control group", ylab = "proportion", cex.names = cexnames,
                ylim = yran , beside = FALSE, col= c("gray","black"), xlab = "Days",
                density = c(1000,30), angle = c(0,45),
                cex.lab = cexlab, cex.axis = cexaxis, cex.main = cexmain)
text(-3, 0.21, "B", pos = 1, cex = 2, xpd = NA)

Nt = length(Treatment)
Nc = length(Control)
K = min(c(Nt,Nc))
probs = (1:K)/(K+1) # quantile levels at which the BQTE-function is estimated
q.type = 1
niter = 1000
ylim = c(-11,2.0)
plot(NULL, ylim = ylim, xlim = c(0,20),
     xlab = "Days in controls", ylab = "BQTE (days)",
     main = "Bootstrap samples",xaxs = "i",
     cex.lab = cexlab, cex.axis = cexaxis, cex.main = cexmain)
grid()
xs = 1:max(Control)
at = (4:15)
boot = matrix(0, nrow = niter, ncol = length(at))
y.grid = seq(-11, 2, 0.5)
counts = matrix(0, nrow = length(y.grid), ncol = length(at))
set.seed(1823)
for(ii in 1:niter){
  Tr = sample(Treatment, size = Nt, replace = TRUE)
  Co = sample(Control, size = Nc, replace = TRUE)
  C.q = as.numeric(quantile(Co, prob = probs, type = q.type))
  T.q = as.numeric(quantile(Tr, prob = probs, type = q.type))
  boot[ii,] = as.numeric(approx(C.q, T.q - C.q, xout = at, rule = 2, ties = mean)$y)
  ind = (at-at[1])*nrow(counts) + sapply(boot[ii,],function(z){which.min(abs(z-y.grid))})
  for(jj in ind){counts[ind] = counts[ind] + 1 }
}

cexs = counts[counts > 0]
cexs = 0.25 + (1.5)*(cexs-1)/max(cexs-1)
kk = which(counts > 0, arr.ind = TRUE)
points(at[kk[,"col"]],y.grid[kk[,"row"]], cex = cexs, col = "black", pch = 19)
text(-2, 6.5, "C", pos = 1, cex = 2, xpd = NA)



set.seed(1) #Set random seed for reproducible results
res.bqte = bqte(Treatment = Treatment, Control = Control,
                at = (4:15), bagging = TRUE) #Estimate BQTE

ate = mean(Treatment) - mean(Control)
#Plot results on outcome scale ('plot.rte = F')
plot_bqte(res.bqte, plot.rbqte = F,
          xlim = c(0,20), ylim = ylim,
          xaxs = "i", col = "red",
          col.ci = "black",
          ylab = "BQTE (days)",
          xlab = "Days in controls", cex.lab = cexlab, cex.axis = cexaxis)
abline(ate, 0, col = "black", lt = 1,  lw = 1) #Add ATE line
text(-2, 6.5, "D", pos = 1, cex = 2, xpd = NA)
title(main = "Estimates from bqte( ) function", cex = cexmain)

if(plot.pdf) dev.off()

#
##
### Figures 2,3 and 4.
##
#

### Choose data set
# 1 = Mossad (Figure 2)
# 2 = Zinc acetat (Figure 3)
# 3 = Carrageenan (Figure 4)
# Figure size 900(width) x 636(height)
data.set = 3
plot.pdf = TRUE
set.seed (18)

# Data set 1, Mossad data
if(data.set == 1){
  filename = "https://raw.githubusercontent.com/mjpirinen/bqte/main/data/Mossad.csv"
  yrange.d = c(-12, 2.25)
  yrange.r = c(-80, 0)
  known.ate = c(-4.0, -43)
  at.value = 4:15
  lt.coord = c(7, 1.2)    # write info texts
  ut.coord = c(12.5, -11)  # of UTBQTE and LTBQTE
  label.coord = matrix(c(0.5, 4, 0.5, 15),byrow = T, ncol = 2)
  dat = read.csv(filename)
  duration = dat$Duration
  status = dat$Zinc
  output.file = "Fig2.pdf" #used only if plot.pdf == TRUE
}

# Data set 2, Zinc acetat data
if(data.set == 2){
  filename = "https://raw.githubusercontent.com/mjpirinen/bqte/main/data/ZnAcet.csv"
  yrange.d = c(-8, 0.5)
  yrange.r = c(-60, 0)
  known.ate = c(-2.7, -36)
  label.coord = matrix(c(1.5, 1, 1.5, 13), byrow = T, ncol = 2)
  at.value = 4:12
  dat = read.csv(filename)
  duration = dat$Duration
  status = dat$Zinc
  output.file = "Fig3.pdf" #used only if plot.pdf == TRUE
}

# Data set 3, carrageenan data
if(data.set == 3){
  filename = "https://raw.githubusercontent.com/mjpirinen/bqte/main/data/Carrageenan.csv"
  yrange.d = c(-10, 1.5)
  yrange.r = c(-70, 30)
  known.ate = c(-2.7, -23)
  at.value = 4:19
  label.coord = matrix(c(-0.5, 2.5, -0.5, 45), byrow = TRUE, ncol = 2)
  dat = read.csv(filename)
  duration = dat$Duration
  status = dat$Carr
  output.file = "Fig4.pdf" #used only if plot.pdf == TRUE
}

if(plot.pdf) pdf(output.file, width = 6.6, height = 4.6, pointsize = 9)


Treatment = duration[status == 1]
Control = duration[status == 0]

plot.range = FALSE
res = bqte(Treatment, Control, at = at.value,
           bqte.conf = 0.95, B = 2000,
           K = length(Control), bagging = TRUE,
           tails = TRUE, discrete.range = plot.range,
           interpolation.method = "spline",
           verbose = TRUE)

xlab = "Days in controls"
pch.sty = 19
par(mfrow = c(2,2))
par(mar = c(5,5,2,2))
labs = c("A","B","C","D")
lab.ii = 1
for(rel.scale in c(FALSE,TRUE)){

  #ate = mean(Treatment) - mean(Control)
  yrange = yrange.d
  ylab.a = "BQTE (Days)"
  ylab.b = "L/UTBQTE (Days)"
  if(rel.scale) {
    #ate = 100*ate/mean(Control)
    yrange = yrange.r
    ylab.a = "Relative BQTE (%)"
    ylab.b = "Relative L/UTBQTE (%)"
  }
  plot_bqte(res, plot.ci = TRUE, plot.rbqte = rel.scale,
            pch = pch.sty, col = "red", col.ci = "red",
            xlim = NULL, ylim = yrange,
            xlab = xlab, ylab = ylab.a,
            xaxs = "r", yaxs = "r",
            xaxp = NULL, yaxp = NULL,
            main = "", cex.axis = 1.3, cex.lab = 1.2,
            plot.range = plot.range,
            col.range = "darkblue")
  if(!is.na(known.ate[1 + as.numeric(rel.scale)])){
  abline(h = known.ate[1 + as.numeric(rel.scale)],
         lty = 1, lwd = 1, col = "black")}
  text(label.coord[1+rel.scale,1], label.coord[1+rel.scale,2], labs[lab.ii], cex = 2, xpd = NA)
  lab.ii = lab.ii + 1

  col.utbqte = "blue"
  col.ltbqte = "magenta4"
  plot_tbqte_both(
    x = res,
    x.utbqte = 0, x.ltbqte = 50,
    plot.ci = TRUE, plot.relative = rel.scale,
    pch.utbqte = 25, pch.ltbqte = 24,
    col.utbqte = col.utbqte,
    col.ltbqte = col.ltbqte,
    col.utbqte.ci = col.utbqte,
    col.ltbqte.ci = col.ltbqte,
    jitter = c(0.2, -0.2),
    xlim = NULL, ylim = yrange,
    xlab = xlab, ylab = ylab.b,
    xaxs = "r", yaxs = "r",
    xaxp = NULL, yaxp = NULL,
    cex = 1, cex.axis = 1.3, cex.lab = 1.2,
    main ="", lwd = 1.2)
  text(label.coord[1+rel.scale,1], label.coord[1+rel.scale,2], labs[lab.ii], cex = 2, xpd = NA)
  lab.ii = lab.ii + 1

  if(data.set == 1 & !rel.scale){ #Add info texts
    text(lt.coord[1], lt.coord[2], col = col.ltbqte,
         "LTBQTE is a lower bound \nfor ATE in the lower tail")
    text(ut.coord[1], ut.coord[2], col = col.utbqte,
         "UTBQTE is an upper bound \nfor ATE in the upper tail")}
  if(!is.na(known.ate[1 + as.numeric(rel.scale)])){
    abline(h = known.ate[1 + as.numeric(rel.scale)],
           lty = 1, lwd = 1, col = "black")}
}


if(plot.pdf) dev.off()

#
##
### 4. Compare BQTE and QTE
##
#


# Make Supplementary Figure S5.



filename.1 = "https://raw.githubusercontent.com/mjpirinen/bqte/main/data/Mossad.csv"
Mossad = read.csv(filename.1)
str(Mossad) # Check data looks OK.
table(Mossad$Duration, Mossad$Zinc)

Treatment = subset(Mossad, Zinc == 1)$Duration
table(Treatment)
Control = subset(Mossad, Zinc == 0)$Duration
table(Control)

tt = table(Treatment)
tc = table(Control)
ran = range(c(as.numeric(names(tt)),as.numeric(names(tc))))
pr = matrix(0, ncol = ran[2]-ran[1]+1, nrow = 3)
pr[1,] = ran[1]:ran[2]
pr[2,as.numeric(names(tc))] = as.numeric(tc)/sum(tc)
pr[3,as.numeric(names(tt))] = as.numeric(tt)/sum(tt)
pr
rownames(pr) = c("days","control","treat")
pr

chosen.i = 8 #to highlight for BQTE
q = 0.5 #to mark for QTE
exact_bqte_discrete(at = chosen.i, x = pr[1,], pr.con = pr[2,], pr.trt = pr[3,], verbose = TRUE)
#Quantile levels from controls: (0.48,0.58].
#Treatment values in the range: 5,...,6.
#Find the probability mass of 5 and 6 from treatment group that corresponds to these quantile levels
cdf.t = cumsum(pr[3,])
pr.t.chosen = rep(0,length(pr[3,])) #This is probability in treated corresponding to the chosen value in controls
pr.t.chosen[5:6] = c(cdf.t[5]-0.48, 0.58-cdf.t[5]) #MUST BE SET MANUALLY
pr.c.chosen = rep(0, length(pr[2,]))
pr.c.chosen[chosen.i] = pr[2,chosen.i]
#sum(pr.t.chosen*1:19)/sum(pr.t.chosen) #average treatment value to be used in BQTE(8) estimate

yran = c(0,max(as.vector(pr[2:3,])))
par(mfcol = c(2,2))
bplot = barplot(rbind(pr[3,]-pr.t.chosen, pr.t.chosen), names.arg = as.integer(pr[1,]),
                main = "Treated", ylim = yran, beside = FALSE, col= c("gray","gold"))
x = bplot[as.numeric(quantile(Treatment,q, type = 1))]
arrows(x , 0.2,x , 0.1, code = 2, col = "red", length = 0.1)
text(-2, 0.2, "A", pos = 1, cex = 3, xpd = NA)

bplot = barplot(rbind(pr[2,]-pr.c.chosen, pr.c.chosen), names.arg = as.integer(pr[1,]),
                main = "Control", ylim = yran , beside = FALSE, col= c("gray","gold"), xlab = "Days")
x = bplot[as.numeric(quantile(Control,q, type = 1))]
arrows(x , 0.2,x , 0.1, code = 2, col = "red", length = 0.1)
text(-2, 0.2, "B", pos = 1, cex = 3, xpd = NA)


set.seed(1) #Set random seed for reproducible results
res.bqte = bqte(Treatment = Treatment, Control = Control, at = (4:15), bagging =FALSE) #Estimate BQTE

res.qte = qte(Treatment = Treatment, Control = Control,
              at = c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), bagging = FALSE) #Estimate QTE


ate = mean(Treatment) - mean(Control)
ylim = c(-10.5,0.5)
#Plot results on outcome scale ('plot.rte = F')
plot_bqte(res.bqte, plot.rbqte = F,
          xlim = c(0,20), ylim = ylim,
          xaxs = "i", col = "red",
          col.ci = "black",
          ylab = "BQTE (Days)",
          xlab = "Control group (Days)")
abline(ate, 0, col = "black", lt = 3,  lw = 2) #Add ATE line
text(-1, 3.5, "C", pos = 1, cex = 3, xpd = NA)

plot_qte(res.qte, ylim = ylim, ylab = "QTE (Days)",
         plot.rqte = FALSE, col = "blue", col.ci = "black")
abline(ate, 0, col = "black", lt = 3,  lw = 2)
text(-0.1, 3.5, "D", pos = 1, cex = 3, xpd = NA)



##
## Check how variance of estimators of BQTE and QTE compare for continuous data
## Supplementary Table S1.
##



Nt = 200
Nc = 100
C.a = 3; C.b = 5
T.a = 1.2; T.b = 2
ps = c(0.2, 0.5, 0.8)
xs = as.numeric(qweibull(ps, C.a, C.b))
ys = as.numeric(qweibull(ps, T.a, T.b))
qte.theor.var = ps*(1-ps)/Nt/(dweibull(ys, T.a, T.b))^2 +
  ps*(1-ps)/Nc/(dweibull(xs, C.a, C.b))^2
bqte.theor.var = ps*(1-ps)/(dweibull(ys, T.a, T.b))^2 * (1/Nc + 1/Nt)
set.seed(87)
niter = 100000
res.bqte = matrix(NA, nrow = niter, ncol = length(ps))
res.qte = matrix(NA, nrow = niter, ncol = length(ps))
for(ii in 1:niter){
  Tr = rweibull(Nt, T.a, T.b)
  Co = rweibull(Nc, C.a, C.b)
  Fc = sapply(xs,function(z){mean(Co <= z)})
  res.bqte[ii,] = as.numeric(quantile(Tr, Fc, type = 1)) - xs
  res.qte[ii,] = as.numeric(quantile(Tr, ps, type = 1)) -
    as.numeric(quantile(Co, ps, type = 1))
}
df = data.frame(
  bqte_var_thr = bqte.theor.var,
  bqte_var_emp = apply(res.bqte,2,var),
  qte_var_thr = qte.theor.var,
  qte_var_emp = apply(res.qte,2,var))
rownames(df) = ps
df


#
## Test spline version on bqte()
#

#Use Mossad data
# (The parts commented out could be used for continuous data.)

filename.1 = "https://raw.githubusercontent.com/mjpirinen/bqte/main/data/Mossad.csv"
Mossad = read.csv(filename.1)
Treatment = subset(Mossad, Zinc == 1)$Duration
Control = subset(Mossad, Zinc == 0)$Duration

tt = table(Treatment)
tc = table(Control)
ran = range(c(as.numeric(names(tt)),as.numeric(names(tc))))
pr = matrix(0, ncol = ran[2]-ran[1]+1, nrow = 3)
pr[1,] = ran[1]:ran[2]
pr[2,as.numeric(names(tc))] = as.numeric(tc)/sum(tc)
pr[3,as.numeric(names(tt))] = as.numeric(tt)/sum(tt)
rownames(pr) = c("days","control","treat")
pr

at = c(4:15)
true.val = sapply(at, function(z){
  exact_bqte_discrete(at = z, x = pr[1,], pr.con = pr[2,], pr.trt = pr[3,], verbose = FALSE)})
Nt = length(Treatment)
Nc = length(Control)

#Nt = 200
#Nc = 100
#C.a = 3; C.b = 5
#T.a = 1.2; T.b = 2
#ps = c(0.2, 0.5, 0.8)
#at = as.numeric(qweibull(ps, C.a, C.b))

set.seed(87)
bqte.lin = matrix(NA, nrow = niter, ncol = length(at))
bqte.spline = matrix(NA, nrow = niter, ncol = length(at))

K = min(c(Nt,Nc))
probs = (1:K)/(K+1) # quantile levels at which the BQTE-function is estimated
q.type = 1
niter = 10000
bqte.lin = matrix(NA, nrow = niter, ncol = length(at))
bqte.spline = matrix(NA, nrow = niter, ncol = length(at))
plot(NULL, ylim = c(-10,0), xlim = c(1,20))

#plot(NULL, xlim = c(0,10), ylim = c(-4,1))
for(ii in 1:niter){
  #Tr = rweibull(Nt, T.a, T.b)
  #Co = rweibull(Nc, C.a, C.b)
  Tr = sample(pr[1,], prob = pr[3,], size = Nt, replace = TRUE)
  Co = sample(pr[1,], prob = pr[2,], size = Nc, replace = TRUE)

  C.q = as.numeric(quantile(Co, prob = probs, type = q.type))
  T.q = as.numeric(quantile(Tr, prob = probs, type = q.type))
  bqte.lin[ii,] = as.numeric(approx(C.q, T.q - C.q, xout = at, rule = 2, ties = mean)$y)
  if(ii < 10) lines(C.q, T.q-C.q, col = "black", lty = 1)
  y = T.q - C.q
  sp.mod = lm(y ~ splines::ns(C.q, df = 4))
  bqte.spline[ii,] = as.numeric(predict(sp.mod,newdata = data.frame(C.q = at)))
  if(ii < 10) lines(C.q, predict(sp.mod,newdata = data.frame(C.q = C.q)), col = "red")
}

#xs = as.numeric(qweibull(probs, C.a, C.b))
#ys = as.numeric(qweibull(probs, T.a, T.b))
#true.val = as.numeric(qweibull(ps, T.a, T.b)) - as.numeric(qweibull(ps, C.a, C.b))

lines(at, true.val, col = "limegreen", lwd = 2) #show the true bqte function
#Compute MSE
cbind(rowSums((true.val - t(bqte.lin))^2)/niter,
      rowSums((true.val - t(bqte.spline))^2)/niter)



lin.bqte = bqte(Treatment, Control, at = NULL,
     bqte.conf = 0.95, B = 2000,
     K = length(Control), bagging = TRUE,
     tails = FALSE,
     discrete.range = FALSE,
     exact.sample.quantiles = TRUE,
     interpolation.method = "linear",
     spline.df = 5,
     verbose = TRUE)

spline.bqte = bqte(Treatment, Control, at = NULL,
                bqte.conf = 0.95, B = 2000,
                K = length(Control), bagging = TRUE,
                tails = FALSE,
                discrete.range = FALSE,
                exact.sample.quantiles = TRUE,
                interpolation.method = "spline",
                spline.df = 5,
                verbose = TRUE)

cbind(
  lin.bqte$res[,2],
  spline.bqte$res[,2])

