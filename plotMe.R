plotMe2 <- function(hList, labels, s = 1, t = 2) {
  n <- nrow(hList[[1]]$x)
  m <- length (hList)
  par(pty = "s")
  hMatch <- matchMe (lapply (hList, function(r)
    r$x))
  hMat <- matrix (0, 0, 2)
  for (j in 1:m) {
    hMat <- rbind(hMat, hMatch[[j]][, c(s, t)])
  }
  plot(hMat,
       xlab = "dim 1",
       ylab = "dim 2",
       col = c(rep("RED", n*(m-1)), rep("BLUE", n)),
       cex = c(rep(1, n*(m-1)), rep(2, n)))
  for (i in 1:n) {
    hLine <- matrix (0, 0, 2)
    for (j in 1:m) {
      hLine <- rbind (hLine, hMatch[[j]][i, c(s, t)])
    }
    lines(hLine)
  }
  text(hMatch[[m]], labels, cex = .75)
}

plotMe1 <- function(hList, labels) {
  n <- length (hList[[1]]$x)
  m <- length (hList)
  blow <- function (x) {
    n <- length (x)
    return (matrix (c(1:n, x), n, 2))
  }
  hMat <- matrix (0, 0, 2)
  for (j in 1:m) {
    hMat <- rbind(hMat, blow(hList[[j]]$x))
  }
  plot(hMat,
       xlab = "index",
       ylab = "x",
       col = c(rep("RED", n*(m-1)), rep("BLUE", n)),
       cex = c(rep(1, n*(m-1)), rep(2, n)))
  for (i in 1:n) {
    hLine <- matrix (0, 0, 2)
    for (j in 1:m) {
      hLine <- rbind (hLine, blow(hList[[j]]$x)[i, ])
      lines(hLine)
    }
  }
  text(blow(hList[[m]]$x), labels, cex = 1.00)
  for (i in 1:n) {
    abline(h = hList[[m]]$x[i])
  }
}

