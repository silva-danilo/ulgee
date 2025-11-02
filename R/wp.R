
# modificate gamlss wp
wp <- function(object=NULL, xvar=NULL, resid=NULL, n.inter=4, xcut.points=NULL, 
               overlap=0, xlim.all=4, xlim.worm=3.5, show.given=TRUE, line=TRUE,
               ylim.all=12*sqrt(1/length(resid)), cex=1, cex.lab=1, pch=16,
               ylim.worm=12*sqrt(n.inter/length(resid)), bg="wheat", 
               col=viridis::viridis(1), bar.bg=c(num="light blue"), ...) 
{
  panel.fun <- function(x, y, col=par("col"), bg=NA, pch=par("pch"), 
                        cex=par("cex"), col.smooth="red", span=2/3, iter=3, ...){
    qq <- as.data.frame(qqnorm(y, plot = FALSE))
    if (any(is.infinite(qq$y))) 
      line <- FALSE
    qq$y <- qq$y - qq$x
    grid(nx = NA, ny = NA, lwd = 2)
    points(qq$x, qq$y, pch = pch, col = col, bg = bg, cex = cex)
    abline(0, 0, lty = 2, col = col)
    abline(0, 1e+05, lty = 2, col = col)
    yuplim <- 10 * sqrt(1/length(qq$y))
    level <- 0.95
    lz <- -xlim.worm
    hz <- xlim.worm
    dz <- 0.25
    z <- seq(lz, hz, dz)
    p <- pnorm(z)
    se <- (1/dnorm(z)) * (sqrt(p * (1 - p)/length(qq$y)))
    low <- qnorm((1 - level)/2) * se
    high <- -low
    {
      no.points <- length(qq$y)
      total.points <<- total.points + no.points
      no.mis <- sum(abs(qq$y) > ylim.worm)
      cat("number of missing points from plot=", no.mis, 
          " out of ", no.points, "\n")
      if (any(abs(qq$y) > ylim.worm)) 
        warning("Some points are missed out ", "\n", 
                "increase the y limits using ylim.worm")
    }
    if (any(abs(qq$x) > xlim.worm)) {
      warning("Some points are missed out ", "\n", 
              "increase the x limits using xlim.worm")
    }
    lines(z, low, lty = 2, lwd = 0.5)
    lines(z, high, lty = 2, lwd = 0.5)
    if (line == TRUE) {
      fit <- lm(y ~ x + I(x^2) + I(x^3), data = qq)
      s <- spline(qq$x, fitted(fit))
      flags <- s$x > -2.5 & s$x < 2.5
      lines(list(x = s$x[flags], y = s$y[flags]), col = col, 
            lwd = 0.5)
      assign("coef1", coef(fit), envir = parent.frame(n = 3))
      assign("coef", c(coef, coef1), envir = parent.frame(n = 3))
    }
  }
  check.overlap <- function(interval) {
    if (!is.matrix(interval)) {
      stop(paste("The interval specified is not a matrix."))
    }
    if (dim(interval)[2] != 2) {
      stop(paste("The interval specified is not a valid matrix.\nThe number 
                 of columns should be equal to 2."))
    }
    crows = dim(interval)[1]
    for (i in 1:(crows - 1)) {
      if (!(abs(interval[i, 2] - interval[i + 1, 1]) < 
            1e-04)) {
        interval[i + 1, 1] = interval[i, 2]
      }
    }
    return(interval)
  }
  get.intervals <- function(xvar, xcut.points) {
    if (!is.vector(xcut.points)) {
      stop(paste("The interval is not a vector."))
    }
    if (any((xcut.points < min(xvar)) | any(xcut.points > 
                                            max(xvar)))) {
      stop(paste("The specified `xcut.points' are not within the range of the 
                 x: (", min(xvar), " , ", max(xvar), ")"))
    }
    extra <- (max(xvar) - min(xvar))/1e+05
    int <- c(min(xvar), xcut.points, (max(xvar) + 2 * extra))
    ii <- 1:(length(int) - 1)
    r <- 2:length(int)
    x1 <- int[ii]
    xr <- int[r] - extra
    if (any(x1 > xr)) {
      stop(paste("The interval is are not in a increasing order."))
    }
    cbind(x1, xr)
  }
  deparen <- function(expr) {
    while (is.language(expr) && !is.name(expr) && deparse(expr[[1L]])[1L] == 
           "(") expr <- expr[[2L]]
    expr
  }
  if (is.null(object) && is.null(resid)) 
    stop(paste("A fitted object with resid() method or the argument resid should
               be used ", "\n", ""))
  resid <- if (is.null(object)) 
    resid
  else resid(object)
  DataExist <- FALSE
  if (!is.null(object) && any(grepl("data", names(object$call))) && 
      !(object$call["data"] == "sys.parent()()")) {
    DaTa <- eval(object$call[["data"]])
    DataExist <- TRUE
  }
  if (!grepl("$", deparse(substitute(xvar)), fixed = T) && 
      !grepl("~", deparse(substitute(xvar)), fixed = T) && 
      DataExist) {
    xvar <- eval(substitute(xvar), envir = as.environment(DaTa))
  }
  if (is.null(xvar)) {
    qq <- as.data.frame(qqnorm(resid, plot = FALSE))
    qq$y <- qq$y - qq$x
    level <- 0.95
    lz <- -xlim.all
    hz <- xlim.all
    dz <- 0.25
    z <- seq(lz, hz, dz)
    p <- pnorm(z)
    se <- (1/dnorm(z)) * (sqrt(p * (1 - p)/length(qq$y)))
    low <- qnorm((1 - level)/2) * se
    high <- -low
    if (any(abs(qq$y) > ylim.all)) {
      warning("Some points are missed out ", "\n", "increase the y limits using 
              ylim.all")
    }
    if (any(abs(qq$x) > xlim.all)) {
      warning("Some points are missed out ", "\n", "increase the x limits using 
              xlim.all")
    }
    plot(qq$x, qq$y, ylab = "Desvio", xlab = "Quantis teóricos", 
         xlim = c(-xlim.all, xlim.all), ylim = c(-ylim.all, 
                                                 ylim.all), cex = cex, pch = pch,
         bg = bg, cex.lab = cex.lab, cex.axis = 1.5)
    lines(z, low, lty = 2)
    lines(z, high, lty = 2)
    if (line == TRUE) {
      fit <- lm(qq$y ~ qq$x + I(qq$x^2) + I(qq$x^3))
      s <- spline(qq$x, fitted(fit))
      flags <- s$x > -xlim.all & s$x < xlim.all
      lines(list(x = s$x[flags], y = s$y[flags]), col = col, lwd = 2)
    }
    return(invisible(coef(fit)))
  }
  if (!is(xvar, "formula")) {
    if (is.factor(xvar)) 
      stop("Use formula for factors i.e. xvar=~f")
    w <- if (is.null(object)) 
      rep(1, length(resid))
    else object$weights
    if (all(trunc(w) == w)) 
      xvar <- rep(xvar, w)
    if (is.null(xcut.points)) {
      given.in <- co.intervals(xvar, number = n.inter, 
                               overlap = overlap)
      if (overlap == 0) 
        given.in <- check.overlap(given.in)
    }
    else {
      given.in <- get.intervals(xvar, xcut.points)
    }
    total.points <- 0
    coef <- coef1 <- NULL
    y.y <- resid
    x.x <- resid
    coplot(y.y ~ x.x | xvar, given.values = given.in, panel = panel.fun, 
           ylim = c(-ylim.worm, ylim.worm), xlim = c(-xlim.worm, 
                                                     xlim.worm), 
           ylab = "Deviation", xlab = "Unit normal quantile", 
           show.given = show.given, bg = bg, pch = pch, cex = cex, 
           bar.bg = bar.bg)
    if (overlap == 0) {
      if (total.points != length(resid)) 
        warning("the total number of points in the plot is not equal \n to the 
                number of observations in y \n")
    }
    mcoef <- matrix(coef, ncol = 4, byrow = TRUE)
    out <- list(classes = given.in, coef = mcoef)
    return(invisible(out))
  }
  if (is(xvar, "formula")) {
    w <- if (is.null(object)) 
      rep(1, length(resid))
    else object$weights
    total.points <- 0
    coef <- coef1 <- NULL
    y.y <- resid
    x.x <- resid
    if (DataExist) {
      coplot(as.formula(paste("y.y~x.x|", as.character(xvar), 
                              sep = "")[2]), data = DaTa, panel = panel.fun, 
             overlap = overlap, number = n.inter, ylim = c(-ylim.worm, 
                                                           ylim.worm), 
             xlim = c(-xlim.worm, xlim.worm), 
             ylab = "Deviation", xlab = "Unit normal quantile", 
             show.given = show.given, bg = bg, pch = pch, 
             cex = cex, bar.bg = bar.bg, ...)
      mcoef <- matrix(coef, ncol = 4, byrow = TRUE)
      if ("given.values" %in% names(vars <- list(...))) 
        given.in <- vars$given.values
      else {
        rhs <- deparen(xvar)[[2L]]
        if (length(rhs) == 3 && rhs[[1]] != "$") {
          a <- eval(deparen(rhs[[2L]]), envir = as.environment(DaTa))
          b <- eval(deparen(rhs[[3L]]), envir = as.environment(DaTa))
          if (length(n.inter) == 1L) 
            n.inter <- rep(n.inter, 2)
          if (length(overlap) == 1L) 
            overlap <- rep(overlap, 2)
          Inter1 <- if (is.factor(a)) 
            levels(a)
          else co.intervals(unclass(a), number = n.inter[1L], 
                            overlap = overlap[1L])
          Inter2 <- if (is.factor(b)) 
            levels(b)
          else co.intervals(unclass(b), number = n.inter[2L], 
                            overlap = overlap[2L])
          given.in <- list(Inter1, Inter2)
        }
        else {
          a <- eval(deparen(rhs), envir = as.environment(DaTa))
          if (length(n.inter) == 1L) 
            n.inter <- rep(n.inter, 2)
          if (length(overlap) == 1L) 
            overlap <- rep(overlap, 2)
          given.in <- if (is.factor(a)) 
            levels(a)
          else co.intervals(unclass(a), number = n.inter[1L], 
                            overlap = overlap[1L])
        }
      }
      out <- list(classes = given.in, coef = mcoef)
      return(invisible(out))
    }
    else {
      coplot(as.formula(paste("y.y~x.x|", as.character(xvar), 
                              sep = "")[2]), panel = panel.fun, overlap=overlap, 
             number = n.inter, ylim = c(-ylim.worm, ylim.worm), 
             xlim = c(-xlim.worm, xlim.worm), ylab = "Deviation", 
             xlab = "Unit normal quantile", show.given = show.given, 
             bg = bg, pch = pch, cex = cex, bar.bg = bar.bg, 
             ...)
      mcoef <- matrix(coef, ncol = 4, byrow = TRUE)
      if ("given.values" %in% names(vars <- list(...))) 
        given.in <- vars$given.values
      else {
        rhs <- deparen(xvar)[[2L]]
        if (length(rhs) == 3 && rhs[[1]] != "$") {
          a <- eval(deparen(rhs[[2L]]))
          b <- eval(deparen(rhs[[3L]]))
          if (length(n.inter) == 1L) 
            n.inter <- rep(n.inter, 2)
          if (length(overlap) == 1L) 
            overlap <- rep(overlap, 2)
          Inter1 <- if (is.factor(a)) 
            levels(a)
          else co.intervals(unclass(a), number = n.inter[1L], 
                            overlap = overlap[1L])
          Inter2 <- if (is.factor(b)) 
            levels(b)
          else co.intervals(unclass(b), number = n.inter[2L], 
                            overlap = overlap[2L])
          given.in <- list(Inter1, Inter2)
        }
        else {
          a <- eval(deparen(rhs))
          if (length(n.inter) == 1L) 
            n.inter <- rep(n.inter, 2)
          if (length(overlap) == 1L) 
            overlap <- rep(overlap, 2)
          given.in <- if (is.factor(a)) 
            levels(a)
          else co.intervals(unclass(a), number = n.inter[1L], 
                            overlap = overlap[1L])
        }
      }
      out <- list(classes = given.in, coef = mcoef)
      return(invisible(out))
    }
  }
}

# plots
par(mar=c(5.5,5.5,2,2), mfrow=c(1,3))
random.pos <- function(fit.model, X){
  pos <- rep(T, nrow(X))
  ids <- unique(fit.model$id)
  random_time <- numeric(nrow(X))
  for(i in ids){
    pos_i <- rep(0, nrow(X))
    time_i <- sample(fit.model$time[fit.model$id==i], 1)
    pos_i[fit.model$time==time_i & fit.model$id==i] <- 1
    random_time <- random_time + pos_i
  }
  return(as.logical(random_time))
}

# plots
set.seed(6666)
set.seed(1245)
set.seed(9999)
set.seed(56840)
wp(resid=fit_1$rq[random.pos(fit_1, X)], cex.lab=1.8, lwd=2, xlim.all=2.5)
wp(resid=fit_1$rq[random.pos(fit_1, X)], cex.lab=1.8, lwd=2, xlim.all=2.5)
wp(resid=fit_1$rq[random.pos(fit_1, X)], cex.lab=1.8, lwd=2, xlim.all=2.5)
