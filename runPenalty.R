


runPenalty <-
  function (w,
            delta,
            lbd,
            p = 2,
            itmax = 10000,
            eps = 1e-10,
            cut = 1e-6,
            write = TRUE,
            verbose = FALSE) {
    m <- length (lbd)
    hList <- as.list (1:m)
    hList[[1]] <-
      smacofPenalty(
        w,
        delta,
        p,
        lbd = lbd[1],
        itmax = itmax,
        eps = eps,
        verbose = verbose
      )
    for (j in 2:m) {
      hList[[j]] <-
        smacofPenalty(
          w,
          delta,
          p,
          zold = hList[[j - 1]]$z,
          lbd = lbd[j],
          itmax = itmax,
          eps = eps,
          verbose = verbose
        )
    }
    mm <- m
    for (i in 1:m) {
      if (write) {
        cat(
          "itel",
          formatC(hList[[i]]$itel, width = 4, format = "d"),
          "lambda",
          formatC(
            hList[[i]]$l,
            width = 10,
            digits = 6,
            format = "f"
          ),
          "stress",
          formatC(
            hList[[i]]$s,
            width = 8,
            digits = 6,
            format = "f"
          ),
          "penalty",
          formatC(
            hList[[i]]$t,
            width = 8,
            digits = 6,
            format = "f"
          ),
          "\n"
        )
      }
      if (hList[[i]]$t < cut) {
        mm <- i
        break
      }
    }
    return(hList[1:mm])
  }

writeSelected <- function(hList, ind) {
  m <- length(hList)
  n <- length(ind)
  mn <- sort(union(union(1:3, ind), m - (2:0)))
  for (i in mn) {
    if (i > m) {
      next
    }
    cat(
      "itel",
      formatC(hList[[i]]$itel, width = 4, format = "d"),
      "lambda",
      formatC(
        hList[[i]]$l,
        width = 10,
        digits = 6,
        format = "f"
      ),
      "stress",
      formatC(
        hList[[i]]$s,
        width = 8,
        digits = 6,
        format = "f"
      ),
      "penalty",
      formatC(
        hList[[i]]$t,
        width = 8,
        digits = 6,
        format = "f"
      ),
      "\n"
    )
  }
}

