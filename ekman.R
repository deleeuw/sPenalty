data (ekman, package = "smacof")
ekman <- as.matrix(1 - ekman)
w <- matrix(1, 14, 14)- diag(14)
ekman <- 2 * ekman / sqrt (sum (w * ekman * ekman))
ekmon <- ekman ^ 3
ekmon <- 2 * ekmon / sqrt (sum (w * ekmon * ekmon))