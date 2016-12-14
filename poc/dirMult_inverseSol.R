a <- runif(50, 0, 10)
x <- runif(50, 0, 10)
n <- 50
y <- n * (a*x) / sum(a*x)

B <- matrix(rep(y/a, length(a)), ncol = length(a)) %*% diag(a)

sol <- eigen(B)

z1 <- Re(sol$vectors[, 1])
z2 <- x
lm(z2 ~ z1)$coeff[2]
sum(x) / sum(z1)

z1 * sum(x) / sum(z1) - x
