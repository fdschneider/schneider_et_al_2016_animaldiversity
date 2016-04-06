######################################################################
#
#   'Animal diversity and ecosystem functioning in dynamic food webs'
#                            ------
#                       Analytical Code
#                        Version 2.0.0
#   last edit: 06.04.2016 by Christian Guill & Florian D. Schneider
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
analyse(webstats$S_c_fin[webstats$cons_igp != 0], webstats$cons_igp[webstats$cons_igp != 0], ylab = "", equation = TRUE)

axis(2, at = 10^seq(-4,0,2), labels = c(expression(10^{-4}),expression(10^{-2}),expression(10^0)), tck = longtick, las = 1)
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




