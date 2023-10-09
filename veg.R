data (vegetables, package = "psych")
veg <- abs(qnorm(as.matrix(veg)))
w <- matrix(1, 9, 9)- diag(9)
veg <- 2 * veg / sqrt (sum (w * veg * veg))