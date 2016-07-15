

sensitivity <- function(x, values, labels = values, data = webstats, modeltype = "log") {

  data$variable <- x
  
  output <- list()
  if(modeltype == "log") {
    output[[1]] <- lm( log(meanB_c) ~ log(S_c_fin)*variable, data = data )
    output[[2]] <- lm( log(meanB_b) ~ log(S_c_fin)*variable, data = data )
    output[[3]] <- lm( log(cons_igp) ~ log(S_c_fin)*variable, data = subset(data, cons_igp != 0) )
    output[[4]] <- lm( log(bas_cons) ~ log(S_c_fin)*variable, data = data )
    output[[5]] <- lm( log(cons_resp) ~ log(S_c_fin)*variable, data = data )
    output[[6]] <- lm( log(bas_resp) ~ log(S_c_fin)*variable, data = data )
  }
  
  if(modeltype == "logit") {
    output[[1]] <- lm( log(meanB_c) ~ boot::logit(S_c_fin)*variable, data = data )
    output[[2]] <- lm( log(meanB_b) ~ boot::logit(S_c_fin)*variable, data = data )
    output[[3]] <- lm( log(cons_igp) ~ boot::logit(S_c_fin)*variable, data = subset(data, cons_igp != 0) )
    output[[4]] <- lm( log(bas_cons) ~ boot::logit(S_c_fin)*variable, data = data )
    output[[5]] <- lm( log(cons_resp) ~ boot::logit(S_c_fin)*variable, data = data )
    output[[6]] <- lm( log(bas_resp) ~ boot::logit(S_c_fin)*variable, data = data )
  }
  
  
  if(modeltype == "gam") {
    output[[1]] <- lm( log(meanB_c) ~ boot::logit(S_c_fin)*variable, data = data )
    output[[2]] <- lm( log(meanB_b) ~ boot::logit(S_c_fin)*variable, data = data )
    output[[3]] <- lm( log(cons_igp) ~ boot::logit(S_c_fin)*variable, data = subset(data, cons_igp != 0) )
    output[[4]] <- lm( log(bas_cons) ~ boot::logit(S_c_fin)*variable, data = data )
    output[[5]] <- lm( log(cons_resp) ~ boot::logit(S_c_fin)*variable, data = data )
    output[[6]] <- lm( log(bas_resp) ~ boot::logit(S_c_fin)*variable, data = data )
  }
  
  par( mar = c(5,3,0.5,1)+.1, oma = c(4,4,2,0), las = 1, lend = 1, bty = "n", cex.lab = 1.5)
  layout(matrix(c(1:6), ncol = 3, byrow = F))
  longtick = -0.06
  shorttick = -0.03
  
  x <- seq(1,90,length = 100)
  
  # a) Predator biomass
  
  analyse(webstats$S_c_fin, webstats$meanB_c, 
          equation = FALSE, ylab = "", plotmodel = TRUE, modeltype = modeltype)
  
  for( i in values) {
    y <- exp(predict(output[[1]], newdata = list(S_c_fin = x, variable = rep(i,length(x))  )))
    lines(x, y , col = "red2", lty = 2)
    text(tail(x,1), tail(y,1), round(i, 2), col = "red2", pos = 4)
  }
  
  axis(2, at = 10^c(0,1,2), labels = c( 1, 10,100), tck = longtick, las = 1)
  axis(2, at = (rep(10^seq(0, 1, 1), each = 10) * seq(1,10, 1)), labels = NA, tck = shorttick)
  mtext("animal biomass", las = 0, padj = -3,side = 2)
  mtext("a)", side = 3, line = 0, at =-10)
  
  # b) Basal biomass
  
  analyse(webstats$S_c_fin, webstats$meanB_b, 
          equation = FALSE, ylab = "", plotmodel = TRUE, modeltype = modeltype)
  
  for( i in values) {
    y <- exp(predict(output[[2]], newdata = list(S_c_fin = x, variable = rep(i,length(x))  )))
    lines(x, y , col = "red2", lty = 2)
    text(tail(x,1), tail(y,1), round(i, 2), col = "red2", pos = 4)
  }
  
  axis(1, at = c(0,50,100), tck = longtick)
  axis(1, at = (seq(0,100, 10)), tck = shorttick, labels = NA)
  axis(2, at = 10^(1:2), labels = c( 10,100), tck = longtick, las = 1)
  axis(2, at = (rep(10^seq(1, 1, 1), each = 10) * seq(1,10, 1)), labels = NA, tck = shorttick)
  mtext("plant biomass", las = 0, padj = -3, side = 2)
  mtext("b)", side = 3, line = 0, at =-10)
  
  # c) intra guild predation
  
  analyse(webstats$S_c_fin[webstats$cons_igp != 0], webstats$cons_igp[webstats$cons_igp != 0], 
          ylab = "", equation = FALSE , plotmodel = TRUE, modeltype = modeltype)
  
  for( i in values) {
    y <- exp(predict(output[[3]], newdata = list(S_c_fin = x, variable = rep(i,length(x))  )))
    lines(x, y , col = "red2", lty = 2)
    text(tail(x,1), tail(y,1), round(i, 2), col = "red2", pos = 4)
  }
  
  axis(2, at = 10^seq(-4,0,2), 
       labels = c(expression(10^{-4}),expression(10^{-2}),expression(10^0)), tck = longtick, las = 1)
  axis(2, at = (rep(10^seq(-4, -1, 1), each = 10) * seq(1,10, 1)), labels = NA, tck = shorttick)
  mtext("intraguild predation", las = 0, padj = -3, side = 2)
  mtext("c)", side = 3, line = 0, at =-10)
  
  # d) total consumption on basal
  
  analyse(webstats$S_c_fin, webstats$bas_cons, 
          equation = FALSE, ylab = "", plotmodel = TRUE, modeltype = modeltype)
  
  for( i in values) {
    y <- exp(predict(output[[4]], newdata = list(S_c_fin = x, variable = rep(i,length(x))  )))
    lines(x, y , col = "red2", lty = 2)
    text(tail(x,1), tail(y,1), round(i, 2), col = "red2", pos = 4)
  }
  
  axis(1, at = c(0,50,100), tck = longtick)
  axis(1, at = (seq(0,100, 10)), tck = shorttick, labels = NA)
  axis(2, at = 10^c(-1:0), labels = c(0.1, 1), tck = longtick, las = 1)
  axis(2, at = c(seq(0.1,1, 0.1), seq(1,5, 1) ), tck =shorttick, labels = NA)
  mtext("consumption on plants", las = 0, padj = -3, side = 2)
  mtext("d)", side = 3, line = 0, at =-10)
  
  # e) Predator respiration
  
  analyse(webstats$S_c_fin, webstats$cons_resp,
          equation = FALSE,  ylab = "", plotmodel = TRUE, modeltype = modeltype)
  
  for( i in values) {
    y <- exp(predict(output[[5]], newdata = list(S_c_fin = x, variable = rep(i,length(x))  )))
    lines(x, y , col = "red2", lty = 2)
    text(tail(x,1), tail(y,1), round(i, 2), col = "red2", pos = 4)
  }
  
  axis(2, at = 10^c(-1:0), labels = c(0.1,1),tck = longtick, las = 1)
  axis(2, at = c(seq(.05,0.1, 0.01),seq(.1,1, length = 10)), tck = shorttick , labels = NA)
  mtext("animal respiration", las = 0, padj = -3, side = 2)
  mtext("e)", side = 3, line = 0, at =-10)
  
  # f) Basal respiration
  
  analyse(webstats$S_c_fin, webstats$bas_resp,
          equation = FALSE,  ylab = "", plotmodel = TRUE, modeltype = modeltype)
  
  for( i in values) {
    y <- exp(predict(output[[6]], newdata = list(S_c_fin = x, variable = rep(i,length(x))  )))
    lines(x, y , col = "red2", lty = 2)
    text(tail(x,1), tail(y,1), round(i, 2), col = "red2", pos = 4)
  }
  
  axis(1, at = c(0,50,100), tck = longtick)
  axis(1, at = (seq(0,100, 10)), tck = shorttick, labels = NA)
  axis(2, at = 10^c(-1:0), labels = c(0.1, 1), tck = longtick, las = 1)
  axis(2, at = c( seq(0.1,1, 0.1), seq(1,5,1)), tck = shorttick , labels = NA)
  mtext("plant respiration", las = 0, padj = -3, side = 2)
  mtext("f)", side = 3, line = 0, at =-10)
  
  mtext("animal species richness", side = 1, outer = TRUE)
  
}



sensitivity_add <- function(response, data, label, 
                            col = "#00FFFF",  
                            x = seq(1,90,length = 100), 
                            modeltype = "log", 
                            addpoints = TRUE, 
                            ...) {
  
  if("log" %in% modeltype) {
    model <- lm( log(data[,response]) ~ log(S_c_fin), data = data )
  }
  
  if("logit" %in% modeltype) {
    model <- lm( log(data[,response]) ~ boot::logit(S_c_fin)*variable, data = data )
  }
  
  if("gam" %in% modeltype) {
    model <- lm( log(data[,response]) ~ boot::logit(S_c_fin)*variable, data = data )
  }
  
  if(addpoints) {
    points(data[,response] ~ data$S_c_fin, pch = 20, cex = 0.5, col = paste0(col,"10") )
  }
  
  y <- exp(predict(model, newdata = list(S_c_fin = x)))
  lines(x, y , col = col, lty = 1, lwd =2, ...)
  text(tail(x,1), tail(y,1), round(label, 2), col = col, pos = 4)
}

sensitivity_plus <- function(data, label, cols, modeltype = "log") {
  
  par(mar = c(5,3,0.5,1)+.1, oma = c(4,4,2,0), las = 1, lend = 1, bty = "n", cex.lab = 1.5)
  layout(matrix(c(1:6), ncol = 3, byrow = F))
  longtick = -0.06
  shorttick = -0.03
  
  x <- seq(1,90,length = 100)
  
  webstats <- subset(webstats, WebID %in% sample(webstats$WebID, 6000) )
  
  # a) Predator biomass
  
  analyse(webstats$S_c_fin, webstats$meanB_c, 
          equation = FALSE, ylab = "", plotmodel = TRUE, modeltype = modeltype)
  
  lapply(1:length(data), function(i)  sensitivity_add("meanB_c", data[[i]], label[i], cols[i]))
  
  axis(2, at = 10^c(0,1,2), labels = c( 1, 10,100), tck = longtick, las = 1)
  axis(2, at = (rep(10^seq(0, 1, 1), each = 10) * seq(1,10, 1)), labels = NA, tck = shorttick)
  mtext("animal biomass", las = 0, padj = -3,side = 2)
  mtext("a)", side = 3, line = 0, at =-10)
  
  # b) Basal biomass
  
  analyse(webstats$S_c_fin, webstats$meanB_b, 
          equation = FALSE, ylab = "", plotmodel = TRUE, modeltype = modeltype)
  
  lapply(1:length(data), function(i) sensitivity_add("meanB_b", data[[i]], label[i], cols[i]))
  
  axis(1, at = c(0,50,100), tck = longtick)
  axis(1, at = (seq(0,100, 10)), tck = shorttick, labels = NA)
  
  axis(2, at = 10^(1:2), labels = c( 10,100), tck = longtick, las = 1)
  axis(2, at = (rep(10^seq(1, 1, 1), each = 10) * seq(1,10, 1)), labels = NA, tck = shorttick)
  mtext("plant biomass", las = 0, padj = -3, side = 2)
  mtext("b)", side = 3, line = 0, at =-10)
  
  
  # c) intra guild predation
  
  analyse(webstats$S_c_fin[webstats$cons_igp != 0], webstats$cons_igp[webstats$cons_igp != 0], 
          ylab = "", equation = FALSE , plotmodel = TRUE, modeltype = modeltype)
  
  lapply(1:length(data), function(i) sensitivity_add("cons_igp", subset(data[[i]], cons_igp != 0), label[i], cols[i]))
  
  axis(2, at = 10^seq(-4,0,2), 
       labels = c(expression(10^{-4}),expression(10^{-2}),expression(10^0)), tck = longtick, las = 1)
  axis(2, at = (rep(10^seq(-4, -1, 1), each = 10) * seq(1,10, 1)), labels = NA, tck = shorttick)
  mtext("intraguild predation", las = 0, padj = -3, side = 2)
  mtext("c)", side = 3, line = 0, at =-10)
  
  # d) total consumption on basal
  
  analyse(webstats$S_c_fin, webstats$bas_cons, 
          equation = FALSE, ylab = "", plotmodel = TRUE, modeltype = modeltype)
  
  lapply(1:length(data), function(i)  sensitivity_add("bas_cons", data[[i]], label[i], cols[i]))
  
  axis(1, at = c(0,50,100), tck = longtick)
  axis(1, at = (seq(0,100, 10)), tck = shorttick, labels = NA)
  axis(2, at = 10^c(-1:0), labels = c(0.1, 1), tck = longtick, las = 1)
  axis(2, at = c(seq(0.1,1, 0.1), seq(1,5, 1) ), tck =shorttick, labels = NA)
  mtext("consumption on plants", las = 0, padj = -3, side = 2)
  mtext("d)", side = 3, line = 0, at =-10)
  
  # e) Predator respiration
  
  analyse(webstats$S_c_fin, webstats$cons_resp,
          equation = FALSE,  ylab = "", plotmodel = TRUE, modeltype = modeltype)
  
  lapply(1:length(data), function(i)  sensitivity_add("cons_resp", data[[i]], label[i], cols[i]))
  
  axis(2, at = 10^c(-1:0), labels = c(0.1,1),tck = longtick, las = 1)
  axis(2, at = c(seq(.05,0.1, 0.01),seq(.1,1, length = 10)), tck = shorttick , labels = NA)
  mtext("animal respiration", las = 0, padj = -3, side = 2)
  mtext("e)", side = 3, line = 0, at =-10)
  
  # f) Basal respiration
  
  analyse(webstats$S_c_fin, webstats$bas_resp,
          equation = FALSE,  ylab = "", plotmodel = TRUE, modeltype = modeltype) # ylim = exp(c(-0.51,1)),
  
  lapply(1:length(data), function(i)  sensitivity_add("bas_resp", data[[i]], label[i], cols[i]))
  
  axis(1, at = c(0,50,100), tck = longtick)
  axis(1, at = (seq(0,100, 10)), tck = shorttick, labels = NA)
  axis(2, at = 10^c(-1:0), labels = c(0.1, 1), tck = longtick, las = 1)
  axis(2, at = c( seq(0.1,1, 0.1), seq(1,5,1)), tck = shorttick , labels = NA)
  mtext("plant respiration", las = 0, padj = -3, side = 2)
  mtext("f)", side = 3, line = 0, at =-10)
  
  mtext("animal species richness", side = 1, outer = TRUE)
  
}
