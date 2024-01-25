# Figure 2 - conditional quantile example

library(latex2exp)

# data
x <- c(0.5, 0.7, 0.86, 1, 1.5, 1.8, 2)
y <- c(0.2, 0.35, 0.53, 0.63, 0.72, 0.87, 1)

ncat <- length(x)
weight <- c(0, (y[-ncat]-y[1]) / (y[ncat-1]-y[1]), 1)
weighted_quantile <- (1-weight)*c(x[1], x) + weight*c(x, x[ncat])

# plot
plot(x, y, xlim=c(0.2, 2.2), ylim=c(0,1), cex = 0.001, 
     xlab = 'y', ylab = 'Probability')
lines(x, y, type="s", col='gray')
lines(weighted_quantile, c(0, y), col="black", lty=1, lwd=2)
lines(c(x[1], x), c(0, y), col="gray35", lty=2)
lines(c(x,x[length(x)]), c(0, y), col="gray35", lty=3)
legend(x = 1.5, y = 0.2, legend=c(TeX('$\\hat{Q}_1(p)$'), TeX('$\\hat{Q}_2(p)$'), TeX('$\\hat{Q}(p)$')),
       col=c("gray35","gray35", "black"), lty=c(2,3,1), cex=1, lwd=c(1,1,2), pt.cex = 2, bty='n')


