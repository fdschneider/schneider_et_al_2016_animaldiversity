###############################################################################
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
###############################################################################


# analyse(): function for full plot and analysis of the relationship between
# animal species richness and ecosystem function
# 
# Usage:
#       analyse( webstats$explanatory, webstats$response)
#
#  The function returns a simple x-y scatterplot and adds a linear model 


analyse <- function(x,y, 
                    xlab = "", ylab = "", 
                    yrange = "data", #"auto", 
                    ylim = "data", 
                    xlim = round(range(x),digits = 1), #c(0,1), # (c(0.25,1))+c(-0.01, 0.01), 
                    xval = seq(min(x), max(x), length = 100), 
                    plotting = TRUE, plotmodel = TRUE, 
                    stats = FALSE, window = 2, 
                    rsquared = FALSE,
                    modeltype = c("log"),
                    equation = FALSE, 
                    axes = FALSE, ...) {
  
  # this part calculates the range of the y-axis to correspond to the range of the plotted model line.
  if(yrange == "auto") {
    mod_range <- predict(model, newdata = data.frame(x = range(xval))) # get range of Y values for given X values
    mod_span <- max(mod_range) - min(mod_range) # calculate the span of the model on the Y-axis
    explained <- summary(model)$r.squared # get R-squared of linear model
    y_span <- mod_span / explained # calculate Y-axis length relative to the model span
    mrgn <- (y_span - mod_span)/2 # calculate upper and lower margins to add
    yrange <- exp(c( min(mod_range) - mrgn , max(mod_range) + mrgn)) # calculate Y-axis range
  }
  
  if(ylim == "auto") ylim <- yrange
  iq95 <- quantile(y, c(0.025,0.975))
  iq99 <- quantile(y, c(0.005,0.995))
  if(ylim == "data") ylim <- 10^(log10(iq99) + c(log10(iq99[2])-log10(iq99[1]) )*c(-.2,0.2))
  
  if(plotting) {
    plot(x,y, log = "y", col = "#00000015",
         pch = 20, cex = 0.6, xaxt = "n", yaxt = "n", bty = "n", 
         xlab = xlab, ylab = ylab, ylim = ylim, xlim = xlim , ...) # plot Y and X values of all food webs
    
  } else {
    plot(NA,NA, log = "y", col = "#00000015", 
         pch = 20, cex = 0.6, xaxt = "n", yaxt = "n", bty = "n", 
         xlab = xlab, ylab = ylab, ylim = ylim, xlim = xlim , ...) # plot Y and X values of all food webs
    
  }
  
  if(stats) {
    
    med <- sapply(1:100, function(i) { median(y[x >= i-window & x <= i+window]) } )
    lower <- sapply(1:100, function(i) { quantile(y[x >= i-window & x <= i+window], c(0.25)) } )
    upper <- sapply(1:100, function(i) { quantile(y[x >= i-window & x <= i+window], c(0.75)) } )
    
    polygon(c(1:100,100:1), c(upper,rev(lower)), col = "#FFFFFF80" , border = NA )
    lines(1:100, upper, lty = 2, col = "#000000")
    lines(1:100, lower, lty = 2, col = "#000000")
    lines(1:100, med, col = "#000000")
    
  }
  
  if(plotmodel) {
    if( "logit" %in% modeltype) {
      model <- lm(I(log(y)) ~ I( boot::logit(x)) ) 
      log_Y <- predict(model, newdata = list(x = xval), type = "response", se.fit = TRUE) # predict values for range in 'xval'
      lines(xval, exp(log_Y$fit), lwd = 2, col = "red") # add linear model
      
    }
    if( "log" %in% modeltype) {
      model <- lm( log(y) ~ log(x) ) 
      log_Y <- predict(model, newdata = list(x = xval), type = "response", se.fit = TRUE) # predict values for range in 'xval'
      lines(xval, exp(log_Y$fit), lwd = 2, col = "red") # add linear model
      
      if(equation) {
        mtext(side = 3, line = -3,
              substitute(
                italic(Y) == a %.% italic(X) ^ b,
                list(a = round( exp(summary(model)$coefficients[1,1]), digits = 2),
                     b = round(summary(model)$coefficients[2,1], digits = 2)
                )
              )
        ) # add model equation
        
      }
      
      if(rsquared) {
        mtext(side = 1, line = -2, at = 80, 
              substitute(
                R ^ 2 == rs,
                list(rs = round(summary(model)[]$r.squared, digits = 2))
              )
        )
      } # plot R squared
    }
    if("gam" %in% modeltype) {
      model <- mgcv::gam(I(log(y)) ~ s(x) )
      log_Y <- predict(model, newdata = list(x = xval), type = "response", se.fit = TRUE) # predict values for range in 'xval'
      lines(xval, exp(log_Y$fit), lwd = 1, col = "blue") # add linear model
      lines(xval, exp(log_Y$fit+1.96*log_Y$se.fit), lwd = 1, lty = 2, col = "blue") # add linear model
      lines(xval, exp(log_Y$fit-1.96*log_Y$se.fit), lwd = 1, lty = 2, col = "blue") # add linear model
    }
  } 
  
  if(axes) {
    yrange <- round(log10(range(y)), 0)
    ydiff <- yrange[2]-yrange[1]
    axis(2, at = seq(yrange[1],yrange[2], 3))
    axis(2, rep(seq(1,10, length = 10), times = ydiff)*rep(10^(1:ydiff), each = 10))
    axis(2)
  }  
  
}
