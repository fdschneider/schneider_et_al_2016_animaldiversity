######################################################################
#
#  Project: Predator Diversity and Ecosystem Functioning
#                            ------
#                       Analytical Code
#                        Version 1.0.0
#          last edit: 30.09.2014 by Florian Schneider
#
# Copyright (C) 2014 Christian Guill & Florian D. Schneider
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact the authors: [Christian Guill](C.P.Guill@uva.nl), 
#                      [Florian D. Schneider](florian.schneider@univ-montp2.fr)
#
######################################################################


rm(list=ls())

##################################
### read in simulation results ###
##################################

# list of filenames
items <- paste("data\\pdef_2.0\\pdef_2.0_", c(paste("0", as.character(1:9), sep ="")  ,as.character(10:55)), sep ="")

# loop: open files and append to object webstats
webstats <- read.table(items[1], header = T, skip = 0, sep = "")
for(i in 2:55) webstats <- rbind(webstats, read.table(items[i], skip =1, header = T, sep = ""))

##################################
###### derived parameters ########
##################################

S_b_ini <- 30 # initial number of basal species
webstats$S_ini <- webstats$S_c+S_b_ini        # initial number of species
webstats$S_c_fin <- round(webstats$S_c*webstats$P_c, digits = 0)  # final number of consumers, creates NaNs !!!!
webstats$S_c_fin[webstats$S_c_fin == "NaN"] <- 0                     # sets NaN = 0
webstats$S_b_fin <- round(S_b_ini*webstats$P_b, digits = 0)		# final number of basal species
webstats$S_fin <- webstats$S_c_fin+webstats$S_b_fin				# final number of species
webstats$FP <- webstats$FD_fin/webstats$FD_ini				# functional persistence
webstats$iniC <- webstats$iniL/webstats$S_ini^2     # initial Connectance
webstats$finC <- webstats$finL/webstats$S_fin^2    # final Connectance

webstats$Mav <- webstats$meanB/(webstats$findens_c+webstats$findens_b) #  average bodymass
webstats$Mav_c <- webstats$meanB_c/webstats$findens_c   # average bodymass of predators
webstats$Mav_b <- webstats$meanB_b/webstats$findens_b  # average bodymass of basalspecies

webstats$meanTLav_ini[which(webstats$meanTLav_ini == "NaN")] <- 0   # replaces NaN with 0
webstats$meanTLav_fin[which(webstats$meanTLav_fin == "NaN")] <- 0   # replaces NaN with 0
webstats$meanTLmin_ini[which(webstats$meanTLmin_ini == "NaN")] <- 0   # replaces NaN with 0
webstats$meanTLmin_fin[which(webstats$meanTLmin_fin == "NaN")] <- 0   # replaces NaN with 0 
webstats$meanTLeff[which(webstats$meanTLeff == "NaN")] <- 0    # replaces NaN with 0 


##################################
###### correct FD for maxFD ######
##################################

# To correct FD for the actual niche space. The niche space is the maximal area that can be covered by the predators' feeding likelihoods (denominator of *FD*; dark area in Fig. 3*b*). It is defined as the whole range of the bodymass axis ($\log_{10}m_\mathit{max} - \log_{10}m_\mathit{min} = 12$) times 1 (i.e., the maximum of the likelihood curve), except for the body-mass range that cannot be accessed even by the potentially largest predator (with $\log_{10}m_i=12$) with the feeding likelihood $\mathcal{L}_\textrm{max}$. 


success <- function(mi, mj, parms.fr) with(parms.fr,  ((mi/(mj*Ropt)) * exp(1-(mi/(mj*Ropt))))^gam )  # feeding likelihood function
parms.fr <- list(Ropt = 100, gam =  2, aj = 0.25, ai = 0.25, a0 = 200, h0= 0.4, hi = -0.75, hj = 0.25)  # parameters of functional response

mids <- seq(0,12,length = 1200) # steps of prey body mass for which to calculate the feeding likelihood
values <- success(10^(12), 10^(mids), parms.fr) # calculatefeeding likelihood for all 'mids'
values[1:which(values == max(values))] <- 1 # set all values left of the max value to 1
#plot(mids, values, type = "h")

maxFD <- sum(values*12/1200)/12  # summed area of maximally possible FD

webstats$FD_ini <- webstats$FD_ini/maxFD # correct initial FD 
webstats$FD_fin <- webstats$FD_fin/maxFD # correct final FD


##################################
##### drawing valid datasets  ####
##################################

# from the 55.550 generated webs, draw 40.000 valid ones where at least one predator survived.  
set.seed(11) 
webstats <- webstats[sort(sample(webstats$WebID[webstats$fintop_b == 0 ], 40000, replace = FALSE)),]   # selects 40.000 random food webs


##################################
#### analyse B-EF relationship ###
##################################

# generate list of linear models of stocks and rates as function of FD

output <- list()

output[[1]] <- lm( log(meanB_c) ~ log(FD_fin), data = webstats )
output[[2]] <- lm( log(meanB_b) ~ log(FD_fin), data = webstats )
output[[3]] <- lm( log(cons_igp) ~ log(FD_fin), data = subset(webstats, cons_igp != 0) )
output[[4]] <- lm( log(bas_cons) ~ log(FD_fin), data = webstats )
output[[5]] <- lm( log(cons_resp) ~ log(FD_fin), data = webstats )
output[[6]] <- lm( log(bas_resp) ~ log(FD_fin), data = webstats )

lapply(output, summary)




# function to plot the food-webs and the linear model

analyse <- function(x,y, xlab = "", ylab = "",  yrange = "auto", ylim = "auto", xlim = (c(0.25,1))+c(-0.01, 0.01), xval = seq(0.3, 0.95, length = 100), ...) {
    
  model <- lm(log(y) ~ log(x)) # fit log-linear model
  log_Y <- predict(model, newdata = list(x = xval), type = "response", se.fit = TRUE) # predict values for range in 'xval' 
  
  # this part calculates the range of the y-axis to correspond to the range of the plotted model line. 
  if(yrange == "auto") {
    mod_range <- predict(model, newdata = data.frame(x = range(xval))) # get range of Y values for given X values
    mod_span <- max(mod_range) - min(mod_range) # calculate the span of the model on the Y-axis
    explained <- summary(model)$r.squared  # get R-squared of linear model
    y_span <- mod_span / explained  # calculate Y-axis length relative to the model span
    mrgn <- (y_span - mod_span)/2  # calculate upper and lower margins to add 
    yrange <- exp(c( min(mod_range) - mrgn , max(mod_range) + mrgn)) # calculate Y-axis range
  }
  if(ylim == "auto") ylim <- yrange
    
  plot(x,y, log = "xy", col = "#00000010", pch = 20, cex = 0.3, xaxt = "n", yaxt = "n", bty = "n", xlab = xlab, ylab = ylab, ylim = ylim, xlim = xlim , ...) # plot Y and X values of all food webs
  
  abline(v = 10^par()$usr[1], lwd = 1 ) # add a solid line along the y-axis
  
  lines(xval, exp(log_Y$fit), lwd = 2, col = "red") # add linear model
  
  mtext(side = 3, line = -3,
       substitute(
         italic(Y) == a %.% italic(X) ^ b, 
         list(a = round(summary(model)$coefficients[1,1], digits = 2), 
              b = round(summary(model)$coefficients[2,1], digits = 2)
              )
         ) 
       )  # add model equation
}


# plotting stocks and rates as function of consumer Functional Diversity to default screen plotting device (re-scaling might be recommended)

par( mar = c(5,3,0.5,1)+.1, oma = c(4,4,2,0), las = 1, lend = 1, bty = "n")
layout(matrix(c(1:6), ncol = 3, byrow = F))

with(webstats, {

  longtick = -0.06
  shorttick = -0.03
  par(cex.lab = 1.5)
  
  
  # a) Predator biomass
  analyse(FD_fin, meanB_c, ylab = "")
  
  axis(2, at = 10^c(-3,0,3), tck = longtick, las = 1)
  axis(2, at = (rep(10^seq(-4, 4, 1), each = 10) * seq(1,10, length = 10)), labels = NA, tck = shorttick)
  mtext("predator biomass", las = 0, padj = -4,side = 2)
 
  # b) Basal biomass
  analyse(FD_fin, meanB_b, ylab = "")
    
  axis(1, at = c(0.3,0.5,1), tck = longtick)
  axis(1, at = log10(seq(0.3,1, 0.1)), tck = shorttick, labels = NA)
  
  axis(2, at = 10^(0:3), tck = longtick, las = 1)
  axis(2, at = (c(seq(1,10, length = 10), seq(10,100, length = 10), seq(100,1000, length = 10))), tck = shorttick , labels = NA)
  mtext("basal biomass", las = 0, padj = -4, side = 2)
  
  
  # c) intra guild predation
  analyse(FD_fin[cons_igp != 0], cons_igp[cons_igp != 0], ylab = "" )
  
  axis(2, at = 10^(-3:0), tck = longtick, las = 1)
  axis(2, at = (c(seq(0.001,.01, length = 10), seq(.01,.1, length = 10), seq(.1,1, length = 10))), tck = shorttick , labels = NA)
  mtext("intraguild predation", las = 0, padj = -4, side = 2)
  
  # d) total consumption on basal
  analyse(FD_fin, bas_cons, ylab = "")
  
  axis(1, at = c(0.3,0.5,1), tck = longtick)
  axis(1, at = log10(seq(0.3,1, 0.1)), tck = shorttick, labels = NA)
  
  axis(2, at = 10^c(-1:0), tck = longtick, las = 1)
  axis(2, at = (c(seq(0.1,1, length = 10) )), tck =shorttick, labels = NA)
  mtext("consumption on basal", las = 0, padj = -4, side = 2)
  
  
  # e) Predator respiration
  analyse(FD_fin, cons_resp, ylab = "")
 
  axis(2, at =  10^c(-2:1), tck = longtick, las = 1)
  axis(2, at = (c(seq(.1,1, length = 10))), tck = shorttick , labels = NA)
  mtext("predator respiration", las = 0, padj = -4, side = 2)
  
  # f) Basal respiration
  analyse(FD_fin, bas_resp, ylab = "") # ylim = exp(c(-0.51,1)),
  
  axis(1, at = c(0.3,0.5,1), tck = longtick)
  axis(1, at = log10(seq(0.3,1, 0.1)), tck = shorttick, labels = NA)
  
  axis(2, at = 10^c(0:1), tck = longtick, las = 1)
  axis(2, at = c( seq(1,10, length = 10)), tck = shorttick , labels = NA)
  mtext("basal respiration", las = 0, padj = -4, side = 2)
  
})

mtext("Functional Diversity", side = 1, outer = TRUE)

