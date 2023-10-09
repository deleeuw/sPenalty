

smacofLoss <- function(d, w, delta) {
  return(sum(w * (delta - d) ^ 2) / 4)
}

smacofBmat <- function(d, w, delta) {
  dd <- ifelse(d == 0, 0, 1 / d)
  b <- -dd * w * delta
  diag(b) <- -rowSums(b)
  return(b)
}

smacofVmat <- function(w) {
  v <- -w
  diag(v) <- -rowSums(v)
  return (v)
}

smacofGuttman <- function (x, b, vinv) {
  return (vinv %*% b %*% x)
}

columnCenter <- function (x) {
  return (apply (x, 2, function (z)
    z - mean (z)))
}

smacofComplement <- function(y, v) {
  return (sum (v * tcrossprod(y)) / 4)
}

smacofPenalty <-
  function(w,
           delta,
           p = 2,
           lbd = 0,
           zold = columnCenter(diag(nrow(delta))),
           itmax = 10000,
           eps = 1e-10,
           verbose = FALSE) {
    itel <- 1
    n <- nrow(zold)
    vmat <- smacofVmat(w)
    vinv <- solve(vmat + (1 / n)) - (1 / n)
    dold <- as.matrix(dist(zold))
    mold <- sum(w * delta * dold) / sum(w * dold * dold)
    zold <- zold * mold
    dold <- dold * mold
    yold <- zold[, (p + 1):n]
    sold <- smacofLoss(dold, w, delta)
    bold <- smacofBmat(dold, w, delta)
    told <- smacofComplement(yold, vmat)
    uold <- sold + lbd * told
    repeat {
      znew <- smacofGuttman(zold, bold, vinv)
      ynew <- znew[, (p + 1):n] / (1 + lbd)
      znew[, (p + 1):n] <- ynew
      xnew <- znew[, 1:p]
      dnew <- as.matrix(dist(znew))
      bnew <- smacofBmat(dnew, w, delta)
      tnew <- smacofComplement(ynew, vmat)
      snew <- smacofLoss(dnew, w, delta)
      unew <- snew + lbd * tnew
      if (verbose) {
        cat(
          "itel ",
          formatC(itel, width = 4, format = "d"),
          "sold ",
          formatC(
            sold,
            width = 10,
            digits = 6,
            format = "f"
          ),
          "snew ",
          formatC(
            snew,
            width = 10,
            digits = 6,
            format = "f"
          ),
          "told ",
          formatC(
            told,
            width = 10,
            digits = 6,
            format = "f"
          ),
          "tnew ",
          formatC(
            tnew,
            width = 10,
            digits = 6,
            format = "f"
          ),
          "uold ",
          formatC(
            uold,
            width = 10,
            digits = 6,
            format = "f"
          ),
          "unew ",
          formatC(
            unew,
            width = 10,
            digits = 6,
            format = "f"
          ),
          "\n"
        )
      }
      if (((uold - unew) < eps) || (itel == itmax)) {
        break
      }
      itel <- itel + 1
      zold <- znew
      bold <- bnew
      sold <- snew
      told <- tnew
      uold <- unew
    }
    zpri <- znew %*% svd(znew)$v
    xpri <- zpri[, 1:p]
    return(list(
      x = xpri,
      z = zpri,
      b = bnew,
      l = lbd,
      s = snew,
      t = tnew,
      itel = itel
    ))
  }
