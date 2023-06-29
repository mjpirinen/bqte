# Example analyses using bqte(), qte() and plotting the results using 
# plot_bqte(), plot_qte() and plot_tbqte functions.
# Harri Hemila, Matti Pirinen
# June 29, 2023.

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
plot_tbqte(res.bqte, utbqte = T, plot.relative = F, pch = 19,
          xlim = c(0,20), xaxs = "i", col = "red",
          col.ci = "black")
abline(ate, 0, col = "blue", lt = 3,  lw = 2) #Add known mean effect line
plot_tbqte(res.bqte, utbqte = T, plot.relative = T, pch = 19,
           xlim = c(0,20), xaxs = "i", col = "red",
           col.ci = "black")
abline(rate, 0, col = "blue", lt = 3, lw = 2) #Add known mean effect line

#PLOT LTBQTEs
plot_tbqte(res.bqte, utbqte = F, plot.relative = F, pch = 19,
           xlim = c(0,20), xaxs = "i", col = "red",
           col.ci = "black")
abline(ate, 0, col = "blue", lt = 3,  lw = 2) #Add known mean effect line
plot_tbqte(res.bqte, utbqte = F, plot.relative = T, pch = 19,
           xlim = c(0,20), xaxs = "i", col = "red",
           col.ci = "black")
abline(rate, 0, col = "blue", lt = 3, lw = 2) #Add known mean effect line


# We can also estimate quantile treatment effects (QTEs) and plot them as follows

res.qte = qte(Treatment = Treatment, Control = Control) #Estimate QTE

#Check output of qte() function
res.qte

par(mfrow = c(1,2))
plot_qte(res.qte, plot.rqte = F)
plot_qte(res.qte, plot.rqte = T)


#
##
###
##
#


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


#
##
###
##
#


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

