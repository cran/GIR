### R code from vignette source 'GIR-guide.Rnw'

###################################################
### code chunk number 1: quick-fit
###################################################
library(GIR)
d <- data.frame(y = rep(c(1, 2, 2, 4, 5), 8),
                group = factor(rep(c("a", "b"), each = 20)))
fit <- gir(y ~ 0 + group, data = d, discrete = "group")
coef(fit)


###################################################
### code chunk number 2: quick-inference
###################################################
bb <- bboot.gir(fit, R = 499L, seed = 42L)
confint(bb)
predict(fit, data.frame(group = c("a", "b")),
        interval = "confidence", bootstrap = bb)
