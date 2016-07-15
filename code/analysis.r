######################################################################
#
#   'Animal diversity and ecosystem functioning in dynamic food webs'
#                            ------
#                       Analytical Code
#                        Version 2.1.0
#   last edit: 15.07.2016 by Christian Guill & Florian D. Schneider
#
# Website: https://github.com/fdschneider/schneider_et_al_2016_animaldiversity
#
# Contact the authors: [Christian Guill](C.P.Guill@uva.nl), 
#                      [Florian D. Schneider](florian.schneider@univ-montp2.fr)
#
#
#
# MIT License
# 
# Copyright (C) 2016 Christian Guill & Florian D. Schneider
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy 
# of this software and associated documentation files (the "Software"), to deal 
# in the Software without restriction, including without limitation the rights 
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
# copies of the Software, and to permit persons to whom the Software is 
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
# SOFTWARE.
#
#
######################################################################


rm(list=ls())

##################################
### read in simulation results ###
##################################

source("code/data.r")
source("code/analyse.r")

##################################
###     fit linear models      ###
##################################


output <- list()
output[[1]] <- lm( log(meanB_c) ~ log(S_c_fin), data = webstats )
output[[2]] <- lm( log(meanB_b) ~ log(S_c_fin), data = webstats )
output[[3]] <- lm( log(cons_igp) ~ log(S_c_fin), data = subset(webstats, cons_igp != 0) )
output[[4]] <- lm( log(bas_cons) ~ log(S_c_fin), data = webstats )
output[[5]] <- lm( log(cons_resp) ~ log(S_c_fin), data = webstats )
output[[6]] <- lm( log(bas_resp) ~ log(S_c_fin), data = webstats )


lapply(output, summary)


# get compact output table: 

out_table <- sapply(output, function(x) summary(x)$coefficients[2,])

out_table <- rbind( 
  sapply(output, function(x) summary(x)$coefficients[1,1]),
  out_table, 
  sapply(output, function(x) summary(x)[]$r.squared) )

t(round(out_table,digits = 4))

##################################
###   result figure (Fig.4)    ###
##################################


par( mar = c(5,3,0.5,1)+.1, oma = c(4,4,2,0), las = 1, lend = 1, bty = "n")
layout(matrix(c(1:6), ncol = 3, byrow = F))

longtick = -0.06
shorttick = -0.03
par(cex.lab = 1.5)

x <- seq(1,90,length = 100)

# a) Predator biomass
analyse(webstats$S_c_fin, webstats$meanB_c, equation = TRUE, ylab = "")

axis(2, at = 10^c(0,1,2), labels = c( 1, 10,100), tck = longtick, las = 1)
axis(2, at = (rep(10^seq(0, 1, 1), each = 10) * seq(1,10, 1)), labels = NA, tck = shorttick)
mtext("animal biomass", las = 0, padj = -3,side = 2)

# b) Basal biomass
analyse(webstats$S_c_fin, webstats$meanB_b, equation = TRUE, ylab = "")

axis(1, at = c(0,50,100), tck = longtick)
axis(1, at = (seq(0,100, 10)), tck = shorttick, labels = NA)

axis(2, at = 10^(1:2), labels = c( 10,100), tck = longtick, las = 1)
axis(2, at = (rep(10^seq(1, 1, 1), each = 10) * seq(1,10, 1)), labels = NA, tck = shorttick)
mtext("plant biomass", las = 0, padj = -3, side = 2)


# c) intra guild predation
analyse(webstats$S_c_fin[webstats$cons_igp != 0], webstats$cons_igp[webstats$cons_igp != 0], 
        ylab = "", equation = TRUE)

axis(2, at = 10^seq(-4,0,2), 
     labels = c(expression(10^{-4}),expression(10^{-2}),expression(10^0)), 
     tck = longtick, las = 1)
axis(2, at = (rep(10^seq(-4, -1, 1), each = 10) * seq(1,10, 1)), labels = NA, tck = shorttick)
mtext("intraguild predation", las = 0, padj = -3, side = 2)

# d) total consumption on basal

analyse(webstats$S_c_fin, webstats$bas_cons, equation = TRUE, ylab = "")

axis(1, at = c(0,50,100), tck = longtick)
axis(1, at = (seq(0,100, 10)), tck = shorttick, labels = NA)

axis(2, at = 10^c(-1:0), labels = c(0.1, 1), tck = longtick, las = 1)
axis(2, at = c(seq(0.1,1, 0.1), seq(1,5, 1) ), tck =shorttick, labels = NA)
mtext("consumption on plants", las = 0, padj = -3, side = 2)


# e) Predator respiration

analyse(webstats$S_c_fin, webstats$cons_resp, equation = TRUE,  ylab = "")

axis(2, at = 10^c(-1:0), labels = c(0.1,1),tck = longtick, las = 1)
axis(2, at = c(seq(.05,0.1, 0.01),seq(.1,1, length = 10)), tck = shorttick , labels = NA)
mtext("animal respiration", las = 0, padj = -3, side = 2)

# f) Basal respiration

analyse(webstats$S_c_fin, webstats$bas_resp, equation = TRUE,  ylab = "")

axis(1, at = c(0,50,100), tck = longtick)
axis(1, at = (seq(0,100, 10)), tck = shorttick, labels = NA)

axis(2, at = 10^c(-1:0), labels = c(0.1, 1), tck = longtick, las = 1)
axis(2, at = c( seq(0.1,1, 0.1), seq(1,5,1)), tck = shorttick , labels = NA)
mtext("plant respiration", las = 0, padj = -3, side = 2)

mtext("animal species richness", side = 1, outer = TRUE)



##################################
###   result figure (Fig.5)    ###
##################################

output_bodysize <- list()
output_bodysize[[1]] <- lm( log(Mav_c) ~ log(S_c_fin), data = webstats )
output_bodysize[[2]] <- lm( log(Mav_b) ~ log(S_c_fin), data = webstats )

lapply(output_bodysize, summary)

out_table_bodysize <- sapply(output_bodysize, function(x) summary(x)$coefficients[2,])

out_table_bodysize <- rbind( 
  sapply(output_bodysize, function(x) summary(x)$coefficients[1,1]),
  out_table_bodysize, 
  sapply(output_bodysize, function(x) summary(x)[]$r.squared) )

t(round(out_table_bodysize,digits = 4))



longtick = -0.06
shorttick = -0.03
par(mfrow = c(2,1))
par( mar = c(3,4,0.5,1)+.1, oma = c(1,0,1,0), las = 1, lend = 1)

# panel a) average animal body mass

analyse(webstats$S_c_fin, webstats$Mav_c, 
        ylab = "animal body mass", modeltype = "log", plotting = TRUE, 
        equation = !raw, ylim = quantile(webstats$Mav_c, c(0.02,0.975)), 
        stats = TRUE )

axis(2, at = 10^c(2,3,4),
     labels = c(expression(10^{2}),expression(10^{3}),expression(10^4)),  
     tck = longtick, las = 1)
axis(2, at = c(seq(10^2,10^3, length = 10), seq(10^3,10^4, length = 10) ), 
     tck =shorttick, labels = NA)

# panel b) average animal body mass

analyse(webstats$S_c_fin, webstats$Mav_b, 
        ylab = "plant body mass", modeltype = "log", 
        equation = !raw, ylim = quantile(webstats$Mav_b, c(0.025,0.975)), 
        stats = TRUE)

axis(1, at = c(0,50,100), tck = longtick)
axis(1, at = seq(0,100, 10), tck = shorttick, labels = NA)

axis(2, at = 10^c(0,1,2),labels = c(1,10,100),  tck = longtick, las = 1)
axis(2, at = rep(seq(1,10, length = 10), times = 2)*rep(10^(0:1), each = 10), 
     tck =shorttick, labels = NA)





##################################
###   sensitivity analysis     ###
##################################


source("code/sensitivity.r")

# Effect of Hill-exponent of feeding rate on relationship between animal species
# richness and ecosystem functions. Influence of Hill-exponent, $q$, on steady
# state (a) animal biomass, (b) plant biomass, (c) intraguild feeding, (d)
# consumption on plants, (e) animal respiration and (d) plant respiration. Solid
# line is linear model of average $q$; dashed lines are covering model
# predictions for 99% of variation in $q$ (99% interquantile range: 0.04 --
# 0.96).

sensitivity(webstats$hill-1,  quantile(webstats$hill-1, c(0.005,0.995)))

# Effect of consumer body mass in attack rate of feeding rate on relationship
# between animal species richness and ecosystem functions.** Influence of
# Scaling exponent, $\\beta_i$, on steady state (a) animal biomass, (b) plant
# biomass, (c) intraguild feeding, (d) consumption on plants, (e) animal
# respiration and (d) plant respiration. Solid line is linear model of average
# $q$; dashed lines are covering model predictions for 99% of variation in
# $\\beta_i$ (99% interquantile range: 0.37 -- 0.57 ).

sensitivity(webstats$a_c, quantile(webstats$a_c, c(0.005, 0.995)))

# Effect of resource body mass in attack rate of feeding rate on relationship
# between animal species richness and ecosystem functions.** Influence of
# Scaling exponent, $\\beta_j$, on steady state (a) animal biomass, (b) plant
# biomass, (c) intraguild feeding, (d) consumption on plants, (e) animal
# respiration and (d) plant respiration. Solid line is linear model of average
# $q$; dashed lines are covering model predictions for 99% of variation in
# $\\beta_j$ (99% interquantile range:  0.07 -- 0.23  ).

sensitivity(webstats$a_p, quantile(webstats$a_p, c(0.005,0.995)))

#Effect of predator interference in attack rate of feeding rate on relationship
#between animal species richness and ecosystem functions.** Influence of
#interference coefficient, $c$, on steady state (a) animal biomass, (b) plant
#biomass, (c) intraguild feeding, (d) consumption on plants, (e) animal
#respiration and (d) plant respiration. Solid line is linear model of average
#$q$; dashed lines are covering model predictions for 99% of variation in $c$
#(99% interquantile range:  0.3 -- 1.3  ).

sensitivity(webstats$pred_int, quantile(webstats$pred_int, c(0.005,0.995)) )

#Effect of consumer body mass scaling in handling time of feeding rate on
#relationship between animal species richness and ecosystem functions.**
#Influence of scaling exponent, $\\eta_i$, on steady state (a) animal biomass,
#(b) plant biomass, (c) intraguild feeding, (d) consumption on plants, (e)
#animal respiration and (d) plant respiration. Solid line is linear model of
#average $q$; dashed lines are covering model predictions for 99% of variation
#in $\\eta_i$ (99% interquantile range:  $-0.56$ -- $-0.41$  ).

sensitivity(webstats$h_c, quantile(webstats$h_c, c(0.005,0.995)))

#Effect of resource body mass scaling in handling time of feeding rate on
#relationship between animal species richness and ecosystem functions.**
#Influence of scaling exponent, $\\eta_j$, on steady state (a) animal biomass,
#(b) plant biomass, (c) intraguild feeding, (d) consumption on plants, (e)
#animal respiration and (d) plant respiration. Solid line is linear model of
#average $q$; dashed lines are covering model predictions for 99% of variation
#in $\\eta_j$ (99% interquantile range:  $-0.71$ -- $-0.61$  ).

sensitivity(webstats$h_p, quantile(webstats$h_p, c(0.005,0.995)))

# Effect of number of plant species richness on relationship between animal
# species richness and ecosystem functions.** Influence of the number of initial
# plant species , $S_P$, on steady state (a) animal biomass, (b) plant biomass,
# (c) intraguild feeding, (d) consumption on plants, (e) animal respiration and
# (d) plant respiration. Yellow and blue lines show alternative simulations with
# low ($S_P = 10$) and high ($S_P = 50$) initial plant species richness,
# respectively.

webstats_Sb10 <- read.table("data/webstats_Sb10.txt", sep = ",")
webstats_Sb50 <- read.table("data/webstats_Sb50.txt", sep = ",")

sensitivity_plus(list(webstats_Sb10, webstats_Sb50), label = c(10,50), 
                 cols = c("#FFAA00","#0060FF") )

# Effect of having a pronounced trophic structure on the relationship between
# animal species richness and ecosystem functions.** Effect of having 50% strict
# herbivores within the animal community on steady state (a) animal biomass, (b)
# plant biomass, (c) intraguild feeding, (d) consumption on plants, (e) animal
# respiration and (d) plant respiration. Blue line show an alternative
# simulation with 50% of animals strictly limited to herbivorous consumption.

webstats_H50 <- read.table("data/webstats_H50.txt", sep = ",")

sensitivity_plus(data = list(webstats_H50), label=  c(50), cols = c("#0060FF") )

# Effect of nutrient supply on relationship between animal species richness and
# ecosystem functions.** Influence of nutrient turnover rate, $D$, on steady
# state (a) animal biomass, (b) plant biomass, (c) intraguild feeding, (d)
# consumption on plants, (e) animal respiration and (d) plant respiration.
# Yellow and blue lines show alternative simulations with low ($D = 0.1$) and
# high ($D = 0.5$) turnover rates, respectively.

webstats_D01 <- read.table("data/webstats_D01.txt", sep = ",")
webstats_D05 <- read.table("data/webstats_D05.txt", sep = ",")

sensitivity_plus(list(webstats_D01, webstats_D05), label = c(0.1,0.5), 
                 cols = c("#FFAA00","#0060FF") )
