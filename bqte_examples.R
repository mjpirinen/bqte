# Example analyses using bqte(), qte() and plotting the results using 
# plot_bqte(), plot_qte(), plot_tbqte() plot_tbqte_both() functions.
# Harri Hemila, Matti Pirinen
# June 17, 2024.

# Contents:
# 1. Read in the functions
# 2. Examples of the functions 
# 3. Figures 2,3 and 4 of the main paper


#
##
### 1. Read in the functions
##
#

# Read in bqte functions by source(filename) command
# where filename is path to the file 'bqte_functions.R'
# in your computer.
source("bqte_functions.R")
# Alternatively, you can read the file directly from GitHub by using
filename = "https://raw.githubusercontent.com/mjpirinen/bqte/main/bqte_functions.R" 
source(filename)

# Data files can be found from GitHub 
# https://github.com/mjpirinen/bqte
# Below direct web links are provided.



#
##
### 2. Examples of the functions 
##
#



# Analysis of Mossad data on effect of zinc gluconate lozenges 
# on duration of common cold.
# Data taken from pages 2-4 of the additional file 2 of
# https://doi.org/10.1186/s12874-017-0356-y

filename.1 = "https://raw.githubusercontent.com/mjpirinen/bqte/main/Mossad.csv"
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

#Plot results on outcome scale ('plot.rte = F')
plot_bqte(res.bqte, plot.rbqte = F, 
          xlim = c(0,20), xaxs = "i", col = "red",
          col.ci = "black",
          ylab = "BQTE (Days)", 
          xlab = "Duration in the control group (Days)")
abline(ate, 0, col = "blue", lt = 3,  lw = 2) #Add ATE line
text(18, 0, "A", pos = 1, cex = 3) #Label this one as panel "A"

#Plot results on relative scale ('plot.rte = T')
plot_bqte(res.bqte, plot.rbqte = T, 
          xlim = c(0,20), xaxs = "i", col = "red",
          col.ci = "black",
          xlab="Duration in the control group (Days)")
abline(rate, 0, col = "blue", lt = 3, lw = 2) #Add relative ATE  line
text(18, 0, "B", pos = 1, cex = 3) #Label this one as panel "B"

#
# Note: Same figures but with confidence bands can be drawn by setting 'pch = NULL'
#

plot_bqte(res.bqte, plot.rbqte = F, pch = NULL,
          xlim = c(0,20), xaxs = "i", col = "red",
          col.ci = "gray",
          ylab = "BQTE (Days)", 
          xlab = "Duration in the control group (Days)")
abline(ate, 0, col = "blue", lt = 3,  lw = 2) #Add ATE line
text(18, 0, "A", pos = 1, cex = 3) #Label this one as panel "A"

#Plot results on relative scale ('plot.rte = T')
plot_bqte(res.bqte, plot.rbqte = T, pch = NULL,
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
plot_tbqte(res.bqte, utbqte = T, plot.relative = F, pch = 25,
          xlim = c(0,20), xaxs = "i", col = "blue",
          col.ci = "blue")
abline(ate, 0, col = "black", lt = 3,  lw = 2) #Add known mean effect line
plot_tbqte(res.bqte, utbqte = T, plot.relative = T, pch = 25,
           xlim = c(0,20), xaxs = "i", col = "blue",
           col.ci = "blue")
abline(rate, 0, col = "black", lt = 3, lw = 2) #Add known mean effect line

#PLOT LTBQTEs
plot_tbqte(res.bqte, utbqte = F, plot.relative = F, pch = 24,
           xlim = c(0,20), xaxs = "i", col = "purple",
           col.ci = "purple")
abline(ate, 0, col = "black", lt = 3,  lw = 2) #Add known mean effect line
plot_tbqte(res.bqte, utbqte = F, plot.relative = T, pch = 24,
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

#Check output of qte() function
res.qte

par(mfrow = c(1,2))
plot_qte(res.qte, plot.rqte = F, col = "limegreen", col.ci = "limegreen")
abline(ate, 0, col = "black", lt = 3,  lw = 2) #Add known mean effect line

plot_qte(res.qte, plot.rqte = T, col = "limegreen", col.ci = "limegreen")
abline(rate, 0, col = "black", lt = 3,  lw = 2) #Add known mean effect line


####

# Analysis of data of the three randomized trials 
# on zinc acetate lozenges on duration of common cold.
# The data are from page 8 of supplementary file 2 of
# https://doi.org/10.1093/ofid/ofx059

filename.2 = "https://raw.githubusercontent.com/mjpirinen/bqte/main/ZnAcet.csv"
ZnAcet = read.csv(filename.2)
str(ZnAcet)
table(ZnAcet$Duration, ZnAcet$Zinc)

Treatment = subset(ZnAcet, Zinc == 1)$Duration
Control = subset(ZnAcet, Zinc == 0)$Duration

ate = mean(Treatment) - mean(Control) #average treatment effect
rate = 100*ate/mean(Control)
c(ate = ate, rel.ate = rate)

set.seed(2)
res.bqte = bqte(Treatment = Treatment, Control = Control, tails = T)

res.bqte

par(mfrow = c(1,2))
plot_bqte(res.bqte, plot.rbqte = F,
          xaxs = "i", xlim = c(0,17), col = "red", 
          col.ci = "black",
          ylab = "BQTE (Days)", 
          xlab="Duration in the control group (Days)")
abline(ate, 0, col = "blue", lt = 3,  lw = 2)
text(15, 0, "A", pos = 1, cex = 3)

plot_bqte(res.bqte, plot.rbqte = T,
          xaxs = "i", xlim = c(0,17), col = "red",
          col.ci = "black",
          xlab="Duration in the control group (Days)")
abline(rate, 0, col = "blue", lt = 3, lw = 2)
text(15, 0, "B", pos = 1, cex = 3)




####

# Analysis of data of the two nasal iota-carrageenan trials 
# on duration of common cold.
# The data are from page 3 of the supplementary file of
# https://doi.org/10.1002/prp2.810
# sensored data on day 20 are imputed with day 20 
# (6 in carrageenan, 21 in placebo).


filename.3 = "https://raw.githubusercontent.com/mjpirinen/bqte/main/Carrageenan.csv"
Carrageenan = read.csv(filename.3)
str(Carrageenan)
table(Carrageenan$Duration,Carrageenan$Carr)

Treatment = subset(Carrageenan, Carr == 1)$Duration
Control = subset(Carrageenan, Carr == 0)$Duration

ate = mean(Treatment) - mean(Control) #average treatment effect
rate = 100*ate/mean(Control)
c(ate = ate, rel.ate = rate)


set.seed(3)
res.bqte = bqte(Treatment = Treatment, Control = Control, tails =T,
                at = c(3:19)) #at specifies the estimation points
res.bqte

par(mfrow = c(1,2))
plot_bqte(res.bqte, plot.rbqte = F, 
          xlim = c(0,20), ylim = c(-9,1),
          yaxp = c(-8, 0, 4),
          col = "red", xaxs = "i",
          col.ci = "black",
          ylab = "BQTE (Days)", 
          xlab="Duration in the control group (Days)")
abline(ate, 0, col = "blue", lt = 3, lw = 2)
text(18, 0, "A", pos = 1, cex = 3)

plot_bqte(res.bqte, plot.rbqte = T, 
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


filename.4 = "https://raw.githubusercontent.com/mjpirinen/bqte/main/Carrageenan30.csv"
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
plot (res.bqte$at, res.bqte$bqte, type ="l", xlim = c(0,20), 
      yaxp = c(-10, 0, 5), xaxs = "i",
      ylab = "Mean treatment effect (Days)", 
      xlab="Duration in the control group (Days)")
lines(res30.bqte$at, res30.bqte$bqte, lty = 2, col = "red")
legend("bottomleft", col = c("black","red"), lty = c(1,2), 
       leg = c("cens. to 20", "cens. to 30"))
# Conclusion: 
# No difference observed due to imputation in the range considered here.



#
##
###  3. Figures 2,3 and 4 of the main paper
##
#



### Choose data set
# 1 = Mossad (Figure 2)
# 2 = Zinc acetat (Figure 3)
# 3 = Carrageenan (Figure 4)
# Figure size 900(width) x 636(height)
data.set = 3
set.seed (18)

# Data set 1, Mossad data
if(data.set == 1){
  filename = "https://raw.githubusercontent.com/mjpirinen/bqte/main/Mossad.csv"
  yrange.d = c(-12, 1)
  yrange.r = c(-80, 0)
  known.ate = c(-4.0, -43)
  at.value = NULL
  lt.coord = c(14, 0.15)    # write info texts
  ut.coord = c(6, -10)  # of UTBQTE and LTBQTE
  label.coord = matrix(c(-1,2,-1, 15),byrow = T, ncol = 2)
  dat = read.csv(filename)
  duration = dat$Duration
  status = dat$Zinc
}

# Data set 2, Zinc acetat data
if(data.set == 2){
  filename = "https://raw.githubusercontent.com/mjpirinen/bqte/main/ZnAcet.csv"
  yrange.d = c(-8, 0.5)
  yrange.r = c(-60, 0)
  known.ate = c(-2.7, -36)
  label.coord = matrix(c(-0.25,1.7,-0.25, 15),byrow = T, ncol = 2)
  at.value = NULL
  dat = read.csv(filename)
  duration = dat$Duration
  status = dat$Zinc
}

# Data set 3, carrageenan data
if(data.set == 3){
  filename = "https://raw.githubusercontent.com/mjpirinen/bqte/main/Carrageenan.csv"
  yrange.d = c(-10, 1.5)
  yrange.r = c(-70, 30)
  known.ate = c(-2.7, -23)
  at.value = 3:19
  label.coord = matrix(c(-1.7,3,-1.7, 50),byrow = T, ncol = 2)
  
  dat = read.csv(filename)
  duration = dat$Duration
  status = dat$Carr
}

Treatment = duration[status == 1]
Control = duration[status == 0]

res = bqte(Treatment, Control, at = at.value, 
           bqte.conf = 0.95, B = 2000, 
           K = length(Control), bagging = TRUE,
           tails = TRUE,
           verbose = FALSE)

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
            main = "", cex.axis = 1.3, cex.lab = 1.2)
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

